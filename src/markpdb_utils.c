#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
#include <regex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "atypes.h"
#include "pdb.h"


/*
 * intranges
 */

int
intrange_findrange(int val, intrange *R, int n)
{
  int left, right, k;

  left = 0; right = n-1;
  while (left <= right) {
    k = (left + right) / 2;
    if (R[k].lower <= val && val <= R[k].upper)
      return k;
    else if (val < R[k].lower)
      right = k-1;
    else if (val > R[k].upper)
      left = k+1;
    else {
      fprintf(stderr, "intrange_findrange() ERROR\n");
      exit(-1);
    }
  }
  return -1;
}

/*
 * given selections (selecstr) for sorted ranges, intrange_decode()
 * returns the num of pairs [lower,upper] stored on array R.
 * Selections are expressed either in digits or in alphabets.
 * Each minus sign must be escaped with '\'
 *	expamples: 3,5-10,12-  a-g,i-n  \-3-10,12-
 */
int
intrange_decode(char *selecstr, intrange *R)
{
  int nranges, lower, upper, digits, c, interval;
  char *bp, *cp, *bp2, buf[BUFSIZ], buf2[BUFSIZ];

  if (selecstr == 0) {
    R[0].lower = INT_MIN; R[0].upper = INT_MAX;
    return nranges = 1;
  }

  strcpy(buf, selecstr);
  nranges = 0;

  bp = buf;
  for (bp = buf; isspace(*bp); bp++)
    ;

  for (bp2 = bp; *bp2; bp2++)
    if (isdigit(*bp2)) {
      digits = 1; break;
    } else if (isalpha(*bp2)) {
      digits = 0; break;
    } else if (*bp2 == '\\' && *(bp2+1) == '-') {
      digits = 1; break;
    }
  if (!*bp2)
    return -2;
  if (digits)
    while ((bp = strtok_r(bp, ",", &cp)) || nranges == 0) {

      lower = INT_MIN; upper = INT_MAX;         // in case unspecified
      if (bp != NULL) {
        bp2 = strcpy(buf2, bp);
        bp = strchr(bp2, '-');
        if (bp > bp2 && bp[-1] == '\\')         // minus sign
          bp = strchr(bp+1, '-');
        if (bp == NULL) {                 // sigle value
	  if (*bp2 == '\\') bp2++;
	  lower = upper = atoi(bp2);
        } else {
          if (bp2 != bp){                  // no leading '-'
            if (*bp2 == '\\') bp2++;
            lower = atoi(bp2);
          }
          bp++;
          if (*bp == '\\' ) bp++;
          if (strlen(bp) > 0) upper = atoi(bp);
        }
      }

      if (lower == INT_MIN && nranges != 0)
        return -1;

      if (nranges > 0 && R[nranges-1].upper == lower - 1)
	R[nranges-1].upper = upper;
      else {
	R[nranges].lower = lower;
	R[nranges].upper = upper;
	nranges++;
      }
      bp = NULL;
    }
  else {                // if (!digits).
    lower = INT_MIN; upper = INT_MAX;
    interval = 0;
    while (1) {
      c = *bp++;

      if (c == ',') continue;           // commas may be omitted.

      if (c == 0 || c > ((long)upper + 1)) {
        R[nranges].lower = lower;
        if (interval) {
          R[nranges].upper = (c) ? c : INT_MAX;
        } else
          R[nranges].upper = upper;
        nranges++;
        if (c == 0) break;
        lower = c; upper = INT_MAX;
      }

      if (c == '-') {
        // prefix hyphen -> lower = 0
        // postfix hyphen -> upper = INT_MAX
        upper = INT_MAX;
        interval = 1;
        continue;
      }

      if (c == upper + 1)
        upper++;
      else if (interval) {
        upper = c;
        interval = 0;
      } else
        lower = upper = c;
      printf("c = %d\n", c);
    }
  }

  return nranges;
}


/*
 * is simplified version of landmark_getatoms().
 * On exit, L[0..n-1] contains n groups of numbers stored in R[].
 */
int
lists_of_intranges(char *str, landmark *L, intrange *R)
{
  landmark *lp;
  intrange *rp;
  int n;
  char buf[BUFSIZ], buf2[BUFSIZ], *bp, *cp;

  n = 0;
  lp = L; rp = R;
  strcpy(buf, str);
  bp = buf;					// bp = gr;gr;..;gr
  while ((bp = strtok_r(bp, ";", &cp)) != 0) {
    bp = strcpy(buf2, bp);			// bp = gr = ir,ir,..,ir
    lp->range = rp;
    lp->nranges = intrange_decode(bp, rp);
    lp->name[0] = 0;
    rp += lp->nranges;
    lp++;
    n++; 
    bp = NULL;
  }
  return n;
}



/*
 * landmark atoms
 */

/*
 * decode landmarks given in the
 * selecstr: semicolon seperated list of landmarks. example: "O@2-103,105;CA@2-"
 */
int
landmark_getatoms(char *selecstr, landmark *L, intrange *R)
{
  landmark *lp;
  intrange *rp;
  int n;
  char buf[BUFSIZ], buf2[BUFSIZ], *bp, *bp2, *cp, *cp2;

  n = 0;
  lp = L; rp = R;
  strcpy(buf, selecstr);
  bp = buf;
  while ((bp = strtok_r(bp, ";", &cp)) != 0) {
    bp2 = strcpy(buf2, bp);                     // bp2 -> atomname@intranges
    bp = strtok_r(bp2, "@", &cp2);              // bp -> atomname
    strcpy(lp->name, bp);
    bp = strtok_r(NULL, "@", &cp2);             // bp -> intranges
    lp->range = rp;
    lp->nranges = intrange_decode(bp, rp);
    rp += lp->nranges;
    lp++;
    n++; 
    bp = NULL;
  }
  lp->name[0] = 0;              // sentinel?  使用してないと思う
  return n;
}

#define ATOMNAMELEN PDB_ATOM_name_LEN

/*
 * atom type: string -> unsigned int
 */
unsigned int
atomtype_name2uint(char *name)
{
  int i, j;
  unsigned int k;
  char c;

  for (i = j = 0; i < ATOMNAMELEN; i++, j++)
    if (!isspace(name[i])) break;
  c = ' '; k = 0; 
  for ( ; i < ATOMNAMELEN; i++) {
    k = k << 8; 
    if (c != 0) {
      c = name[i];
      if (c != ' ') k = k | c;
    }
  }
  while (j--)
    k <<= 8;
  return k;
}
// to confirm an unsigned int is ATOMNAMELEN bits size or longer
enum {
  is_uint_size_ok = 1 / (sizeof(unsigned int) >= ATOMNAMELEN)
};

#undef ATOMNAMELEN


/*
 * atom type hash table
 */

#define AT_HASHTABSIZ (1<<4)
#define AT_HASHMASK ~((~0U)<<4)
#define AT_HASH_MOD_PRIME 47
static struct atomtype_hashent *at_hashtab[AT_HASHTABSIZ];
#define AT_HASHENTPOOL_SIZ 256
static struct atomtype_hashent at_hashentpool[AT_HASHENTPOOL_SIZ];
static int at_nhashentsused = 0;

static int
atomtype_hashindex(unsigned int key)
{
  return (key % AT_HASH_MOD_PRIME) & AT_HASHMASK;
}

static int
atomtype_hashsearch(char *name)
{
  struct atomtype_hashent *p;
  unsigned int k;

  if (at_nhashentsused == 0)
    return -1;
  k = atomtype_name2uint(name);
  p = at_hashtab[atomtype_hashindex(k)];
  while (p) {
    if (p->atomtype == k)
      return p->index;
    else
      p = p->next;
  }
  return -1;
}

static int
atomtype_hashinsert(char *name)
{
  struct atomtype_hashent *p;
  int i;
  unsigned int k;

  if (at_nhashentsused == 0)            // initialize
    memset(at_hashtab, 0, sizeof(struct atomtype_hashent *) * AT_HASHTABSIZ);
  if ((i = atomtype_hashsearch(name)) >= 0)
    return i;
  if (at_nhashentsused >= AT_HASHENTPOOL_SIZ)
    return -1;
  k = atomtype_name2uint(name);
  i = atomtype_hashindex(k);
  p = &at_hashentpool[at_nhashentsused];
  p->next = at_hashtab[i];
  at_hashtab[i] = p;
  p->atomtype = k;
  p->index = at_nhashentsused;
  return at_nhashentsused++;
}

extern int *LUTbuf;

/*
 * allocate an area for the LUT of landmark ids (serial numbers).
 * atom types (uint)(rows) by residue ids (int)(column).
 * returns the ptr to the area.
 */
int **
landmark_initlut(int resid_min, int resid_max, int nlandmarks)
{ 
  int **T, i, n;

  // LUT index begins from lower_min (rather than resid_min)
  n = resid_max - resid_min + 1;
  if ((T = (int **)malloc(sizeof(int *) * nlandmarks)) == 0) return 0;
  T[0] = malloc(sizeof(int) * nlandmarks * n);
  if (T[0] == 0) return 0;
  LUTbuf = T[0];
  memset(T[0], 0, sizeof(int) * nlandmarks * n);
  T[0] = T[0] - resid_min;
  for (i = 1; i < nlandmarks; i++)
    T[i] = T[i-1] + n;

  return T;
}

/*
 * builds an LUT of landmark ids (serial numbers).
 * atom types (uint)(rows) by residue ids (int)(column)
 * returns the range of residues that is going to be marked.
 */
intrange
landmark_preplut(int resid_max, int nlandmarks, landmark *L, int **T)
{
  int baseresid, i, j, k, l, id, lower, upper;
  intrange range;

  baseresid = 1;
  for (i = 0; i < nlandmarks; i++)
    for (j = 0; j < L[i].nranges; j++)
      if (L[i].range[j].lower > INT_MIN && L[i].range[j].lower < baseresid)
	baseresid = L[i].range[j].lower;

  for (i = 0; i < nlandmarks; i++) {
    k = atomtype_hashinsert(L[i].name);
    k = atomtype_hashsearch(L[i].name);
    for (j = 0; j < L[i].nranges; j++) {
      lower = L[i].range[j].lower;
      if (lower == INT_MIN) lower = baseresid; // == 1
      upper = L[i].range[j].upper;
      if (upper > resid_max) upper = resid_max;
      for (l = lower; l <= upper; l++)
        T[k][l] = -1;
    }
  }
  
  id = 0;
  //j = (baseresid < resid_min) ? baseresid : resid_min;
  j = baseresid;
  for ( ; j <= resid_max; j++)
    for (i = 0; i < nlandmarks; i++) {
      // if (j < baseresid) continue;
      id++;
      if (T[i][j] < 0) T[i][j] = id;
    }

  range.lower = lower;
  range.upper = upper;
  return range;
}


/*
 * structure selection
 */

/*
 * selec examples: "1abc.pdb" "2nmr.pdb:1-3,5-" "3xry.pdb:A,C-D"
 *     NMR: model ID.  CRYST: chain ID
 */
int
pdb_filenameandchains(char selec[], intrange *R, char *filename)
{
  int n;
  char buf[BUFSIZ], *bp, *cp;

  bp = strcpy(buf, selec);
  bp = strtok_r(bp, ":", &cp);
  strcpy(filename, bp);
  bp = strtok_r(NULL, ":", &cp);
  n = intrange_decode(bp, R);
  return n;
}

static int pdbrec(char *keyword, char *rec)
/* returns 1 if `rec` is the type of `keyword`, otherwise 0. */
{
  int i;
  char kbuf[7], rbuf[7];

  for (i = 0; i < 6; i++) {
    if (keyword[i] == 0 || isspace(keyword[i]))
      break;
    kbuf[i] = keyword[i];
  }
  for ( ; i < 6; i++)
    kbuf[i] = ' ';

  for (i = 0; i < 6; i++) {
    if (rec[i] == 0 || isspace(rec[i]))
      break;
    rbuf[i] = rec[i];
  }
  for ( ; i < 6; i++)
    rbuf[i] = ' ';

  for (i = 0; i < 6; i++)
    if (kbuf[i] != rbuf[i])
      return 0;
  return 1;
}

static char prefetchbuf[BUFSIZ], *prefetched = 0;
// scan for a fragment of single chain/model
// *modelid (OUT) == value given if contained a MODEL record, == -1 otherwise.
int
pdb_scanmark(FILE *fp, char *altloc, int *chainid, int *modelid, int *ishetatm,
	     int *residmin, int *residmax, pdbEntry *pe, int **LUT, int *offset,
	     int *breaks)
{
  int cid, mid, i, k, resseq_min, resseq_max, resid, nres, eof, hetatm, mark;
  int print, id;
  static int ichain;
  char rec[BUFSIZ], buf[BUFSIZ+16], c, *p;

  print = (LUT != 0);
  mark = print;
  eof = 0;
  nres = 0;
  cid = 0;
  *altloc = 0;
  hetatm = 0;
  mid = -1;
  resseq_max = -999;		// rewritten if this segment contains landmarks
  resseq_min = 999;		// rewritten if this segment contains landmarks
  while (1) {
    if (prefetched) {
      strcpy(rec, prefetched);
      prefetched = 0;
    } else if (!fgets(rec, BUFSIZ, fp)) {
      eof = 1;  break;
    }

    if (pdbrec("END", rec) || pdbrec("ENDMDL", rec)) {
      if (print) printf("%s", rec);
      break;
    } else if (pdbrec("TER", rec)) {
      if (print) printf("%s", rec);
      break;
    } else if (pdbrec("MODEL", rec)) {
	strncpy(buf, &rec[PDB_MODEL_serial_LOC], PDB_MODEL_serial_LEN);
	buf[PDB_MODEL_serial_LEN] = 0;
	sscanf(buf, "%d", &mid);
	if (print) printf("%s", rec);
	ichain = 1;
    } else if (pdbrec("ATOM", rec) || pdbrec("HETATM", rec)) {
      if (breaks && breaks[ichain]
	  < atoi(p = strndup(&rec[PDB_ATOM_resSeq_LOC], PDB_ATOM_resSeq_LEN))) {
	// if (p) free(p);
	strcpy(prefetchbuf, rec);
	prefetched = prefetchbuf;
	break;
      }
      // if (p) free(p);
      if (cid > 0) {
        if (cid != rec[PDB_ATOM_chainID_LOC]) { // next chain, push back
          strcpy(prefetchbuf, rec);
          prefetched = prefetchbuf;
          break;
        } else if ((rec[0] == 'H') != hetatm) {	// ATOM -> HETATM or vice versa
          strcpy(prefetchbuf, rec);
          prefetched = prefetchbuf;
          break;
	}
      } else {		// chainID == 0. start a new chain 
        cid = rec[PDB_ATOM_chainID_LOC];
	hetatm = rec[0] == 'H';
	resseq_max = 0;
	resseq_min = 1;
	resid = INT_MIN;
	nres = 0;
	if (mark) {
          if (mid >= 0)
            id = mid;
          else if (*modelid >= 0)
            id = *modelid;
          else
            id = cid;
	  for (i = 0; i < pe->nranges; i++)
	    if (pe->range[i].lower <= id && id <= pe->range[i].upper)
	      break;
	  if (i >= pe->nranges)
	    mark = 0;
	}
      }

      strncpy(buf, &rec[PDB_ATOM_resSeq_LOC], PDB_ATOM_resSeq_LEN);
      buf[PDB_ATOM_resSeq_LEN] = 0;
      sscanf(buf, "%d", &k);
      if (breaks)
	k -= breaks[ichain-1];
      if (k != resid) {         // new residue
        nres++;
        if (k > resseq_max)
          resseq_max = k;
	if (k < resseq_min)
          resseq_min = k;
	resid = k;
      }

      c = rec[PDB_ATOM_altLoc_LOC];
      if (c != ' ') {
        for (p = altloc; *p; p++)
          if (*p == c)
            break;
        if (*p == 0) {
	  *p++ = c;
	  *p = 0;
	}
      }

      if (mark) {
	if (rec[0] == 'A') {
	  strncpy(buf, &rec[PDB_ATOM_name_LOC], PDB_ATOM_name_LEN);
	  buf[PDB_ATOM_name_LEN] = 0;
	  if ((i = atomtype_hashsearch(buf)) >= 0 && LUT[i][k] > 0) {
	    if (strlen(rec) < PDBM_EXTRA_x_LOC+1) {	// +1 for trailing nl
	      rec[strlen(rec)-1] = 0;
	      sprintf(buf, "%-*s\n", PDBM_EXTRA_x_LOC, rec);
	      strcpy(rec, buf);
	    }
	    if (offset)
	      sprintf(buf, "@%0*d", PDBM_EXTRA_ID_LEN, LUT[i][k]+offset[cid]);
	    else
	      sprintf(buf, "@%0*d", PDBM_EXTRA_ID_LEN, LUT[i][k]);
	    memcpy(rec+PDBM_EXTRA_FLAG_LOC, buf, 1+PDBM_EXTRA_ID_LEN);
	  }
	}
      }
      if (print) {
	if (offset) rec[PDB_ATOM_chainID_LOC] = ' ';
        if (breaks) rec[PDB_ATOM_chainID_LOC] = 'A' + ichain-1;
	printf("%s", rec);
      }

    }	// if (pdbrec("ATOM", rec) || pdbrec("HETATM", rec))
  }
  if (breaks) {
    *chainid = ichain;
  } else {
    *chainid = cid;
  }
  ichain++;
  *residmin = resseq_min;
  *residmax = resseq_max;
  *modelid = mid;
  // if (!print) *ishetatm = hetatm;
  *ishetatm = hetatm;
  if (nres == 0 && eof)
    return -1;

  return nres;
}
