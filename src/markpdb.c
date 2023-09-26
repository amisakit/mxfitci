#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "atypes.h"
#include "pdb.h"


/*
 * A tool to mark (backbone and Cb) atoms for superposition using mmfitc.
 */

char *merge(char *a, char *b)
{
  char buf[BUFSIZ], *p, *a0 = a;
  strcpy(buf, a);
  p = a;
  a = buf;
  while (*a || *b)
    if (!*a)
      *p++ = *b++;
    else if (!*b)
      *p++ = *a++;
    else if (*a == *b) {
      *p++ = *a++; b++;
    } else if (*a < *b)
      *p++ = *a++;
    else
      *p++ = *b++;
  *p = 0;
  return a0;
}


#define NRES_LIM 10000
/*
#define CHAINPOOL_SIZ (1<<20)
pdbChain pdbchainpool[CHAINPOOL_SIZ];
int nchainsinpool = 0;
*/
pdbChain *pdbchainpool = 0;	/* storage for pdbChains */

#define MCHAINS_MAX 256		/* chains per a model-record or an entry */
#define NATOMTYPES_MAX 32

#define LM_POOLSIZ NATOMTYPES_MAX
#define IR_POOLSIZ 1024

intrange intrangepool[IR_POOLSIZ];
landmark lm_pool[LM_POOLSIZ];

extern int landmark_getatoms(char *selecstr, landmark *L, intrange *R);
extern int pdb_filenameandchains(char selec[], intrange *R, char *filename);
extern int pdb_scanmark(FILE *fp, char *altloc, int *chainid, int *modelid, int *ishetatm, int *residmin, int *residmax, pdbEntry *pe, int **LUT, int *offset, int *breaks);
extern int **
landmark_initlut(int resid_min, int resid_max, int nlandmarks);
extern intrange
landmark_preplut(int resid_max, int nlandmarks, landmark *L, int **T);
extern int intrange_decode(char *selecstr, intrange *R);


int *LUTbuf = 0;

char usagemessage[] =
  "markpdb [-j] [-s BREAKS] [-f fasta offset] \"PATTERNS\" infile[:RANGE]\n"
  "  Option:\n"
  "    -j  join all chains in the file (or in each model block\n"
  "        if given) into a single chain\n"
  "    -s  chains are splited at BREAKS to form new chains\n"
  "        chains 1-100, 201-200, & 300- are generated by 100,200,300\n"
  "    -f  landmark id's are modified according to the sequences in fasta\n"
  "  \"PATTERNS\" (double quotes are necessary)\n"
  "    a list of PATTERNs separated with semicolons\n"
  "    PATTERN is an atom name followed by the ranges of residues,\n"
  "      with '@' between them\n"
  "  RANGE subsequent to infile\n"
  "    ranges of chains, or ranges of model ids\n"
  "  Example:\n"
  "  ./markpdb \"CA@2,4-\" nmr.pdb:1,3-10\n"
  "  ./markpdb \"CA@2,4-127;N@-7;O@1,2-5,6-\" xray.pdb:a,c-\n"
  "  ./markpdb -s 100,200,300 \"CA\" input.pdb\n";

void usage()
{
  printf("%s", usagemessage);
  exit(-1);
}

static FILE *fp = NULL;
static FILE *fasta = NULL;
static char *progname = 0;

void abnormalexit(char *message, int code)
{
  if (message && code)
    fprintf(stderr, "(%d) %s\n", code, message);
  else if (message)
    fprintf(stderr, "%s\n", message);
  if (fp)
    fclose(fp);
  if (LUTbuf) free(LUTbuf);
  exit(code);
}

struct arguments {
  int nbreaks;			/* < 0 : join, > 0 : split into nbreaks chains*/
  int breakpoints[MCHAINS_MAX];
  int debugprint;
  char *pattern;		/* landmark atoms */
  char *in;			/* a pdb file and the range of chains/blocks */
  char *fasta;			/* a fasta file */
  int fasta_off;		/* offset */
} arguments;

int
parse_args(int argc, char *argv[], struct arguments *arguments)
{
  int opt, i;
  intrange irbuf[MCHAINS_MAX];

  arguments->nbreaks = 0;
  arguments->debugprint = 0;
  arguments->fasta = 0;
  while ((opt = getopt(argc, argv, "js:f:D")) != -1)
    switch (opt) {
    case 'j':
      arguments->nbreaks = -1;
      break;
    case 's':
      arguments->nbreaks = intrange_decode(optarg, irbuf);
      for (i = 0; i < arguments->nbreaks; i++)
	arguments->breakpoints[i+1] = irbuf[i].lower;
      arguments->breakpoints[0] = 0;
      arguments->breakpoints[arguments->nbreaks+1]=pow(10,PDB_ATOM_resSeq_LEN);
      break;
    case 'f':
      arguments->fasta = strdup(optarg);
      arguments->fasta_off = atoi(argv[optind++]);
      break;
    case 'D':
      arguments->debugprint = 1;
      break;
    default:
      usage();
    }
  if (argc - optind != 2)
    usage();
  arguments->pattern = argv[optind];
  arguments->in = argv[optind+1];
  return 0;
}
  

/*
 * read the fasta file fp for the sequence of chain given in the file pdb.
 * offset is the position of resid #1 relative to the head of seq given in fp.
 * e.g., if head of the seq is resid #4 then offset = -3.
 * lim is the limit for length of seq[]
 * OUTPUT:
 *  seq[] : sequence given in the fasta records.
 *  raise[] : raise[j] is the raise amount of LUT for the resid j
 * NOTICE:
 *  resID@pdb > 0.  does not work correctly for negative IDs.
 * If lim < 0, scan_fasta will return the num of residues (including gaps)
 * of the first entry in fp (fasta), and never write on seq[] and raise[].
 */
int
scan_fasta(FILE *fp, int offset, char *pdb, char seq[], int raise[], int lim)
{
  char buf[BUFSIZ], *p, *q;
  int found, ires, nres;

  if (lim < 0) {
    found = 0;
    while (fgets(buf, BUFSIZ, fp) && !found) {
      if (buf[0] != '>')
	continue;
      found = 1;
      nres = 0;
      while (fgets(buf, BUFSIZ, fp)) {
	if (buf[0] == '>') break;
	for (q = buf; *q; q++)
	  if (isalpha(*q) || *q == '-')
	    ++nres;
      }
    }
    return nres;
  }

  for (ires = 1; ires <= -offset; ires++)
    raise[ires] = offset;

  nres = -offset;
    
  found = 0;
  while (fgets(buf, BUFSIZ, fp) && !found) {
    if (buf[0] != '>')
      continue;
    p = buf + 1; 
    while ((q = strtok(p, " \t|"))) {
      while (*q && isspace(*q)) q++;
      p = q + strlen(q) - 1;
      while (p > q && isspace(*p)) *p-- = 0;
      if (strcmp(q, pdb) == 0) {
	found = 1; break;
      }
      p = 0;
    }
    if (!found)
      continue;
    seq[0] = ' ';
    p = seq + 1;
    while (fgets(buf, BUFSIZ, fp)) {
      if (buf[0] == '>') break;
      q = buf;
      while (*q) {
	if (isalpha(*q)) {
	  if (++nres > 0)
	    raise[nres] = offset;
	} else if (*q == '-')
	  offset++;
	else
	  break;

	if (--lim < 0)
	  return -1;
	*p = toupper(*q);
	
	p++; q++;
      }
    }
    *p = 0;
  }
  if (!found) {
    fprintf(stderr, "the pdb not found in the fasta.\n");
    exit(-1);
  }
  return nres;
}

int main(int argc, char *argv[])
{
  int i, j, k, nirs, nres, resid_min, resid_max, modelid, ishetatm, nlandmarks;
  int **LUT;
  int topofblock, newchain, mid, debugprint;
  int nblocks, nchains_block, chain_block[256], offset[256];
  // int nres_block[256];
  int *op, nchains, *fasta_raise;
  intrange *irp = intrangepool, rangemarked;
  char buf[BUFSIZ], fasta_seq[NRES_LIM];
  pdbEntry pe;

  progname = argv[0];
  (void) parse_args(argc, argv, &arguments);
  debugprint = arguments.debugprint;

  nlandmarks = landmark_getatoms(arguments.pattern, lm_pool, irp);
  for (i = 0; i < nlandmarks; i++)
    irp += lm_pool[i].nranges;

  if (debugprint) {
    printf("REMARK| Landmark atoms:\n");
    printf("REMARK|");
    for (i = 0; i < nlandmarks; i++) {
      printf("    %s:", lm_pool[i].name);
      for (j = 0; j < lm_pool[i].nranges; j++)
	printf(" %d-%d",
	       lm_pool[i].range[j].lower, lm_pool[i].range[j].upper);
    }
    printf("\n");
  }

  nirs = pdb_filenameandchains(arguments.in, irp, buf);

  if (debugprint) {
    printf("REMARK| the PDB file and chain/model ID selections:\n");
    printf("REMARK|    %s:", buf);
    for (j = 0; j < nirs; j++)
      printf(" %d-%d", irp[j].lower, irp[j].upper);
    printf("\n");
  }

  //strcpy(pe.file, buf);
  pe.file = strdup(buf);
  pe.nranges = nirs;
  pe.range = irp;
  irp += nirs;

  /*
   * scan the pdb file for the max res-id in main chains (ATOM records).
   */
  fp = fopen(pe.file, "r");
  if (fp == NULL)
    fprintf(stderr, "cannot open file %s.\n", pe.file), usage();

  /*
   * estimate the num of chains (equal to or greater than the true number)
   */
  nchains = 0;
  while (1) {
    nres = pdb_scanmark(fp, buf /*altloc*/, &k /*chainid*/, &modelid, &ishetatm,
                        &resid_min, &resid_max, &pe, 0, 0,
                        arguments.nbreaks > 0 ? arguments.breakpoints : 0);
    if (nres < 0) break;
    if (nres == 0) continue;
    nchains++;
  }
  if (arguments.nbreaks > 0) nchains *= (arguments.nbreaks + 1);
  if (debugprint) {
    printf("REMARK| num of fragments = %d ((over)estimate for nchains)\n",
	   nchains);
  }
  if ((pdbchainpool = malloc(sizeof(pdbChain) *	nchains)) == 0)
    fprintf(stderr, "cannot alloc mem for pdbchainpool\n"), exit(-1);


  rewind(fp);
  if (debugprint) {
    printf("REMARK| Chains read:\n");
  }
  /*
   * a block is either
   *  - records enclosed by MODEL and ENDMDL, or
   *  - entire records in a PDB entry that contains no MODEL record
   */
  nblocks = 0;
  nchains_block = 0;		/* num of chains in the top block */

  pe.nchains = 0;
  topofblock = 0;
  pe.chainptr = &pdbchainpool[0];
  mid = modelid = -1;			// model id when unspecified
  while (1) {
    // modelid >= 0 if the prev frag begins with an MODEL record
    nres = pdb_scanmark(fp, buf /*altloc*/, &k /*chainid*/, &modelid, &ishetatm,
			&resid_min, &resid_max, &pe, 0, 0,
			arguments.nbreaks > 0 ? arguments.breakpoints : 0);
    /*
     * on exit:
     * modelid >= 0 if this frag begins with an MODEL record
     * ichain == chainID char (serial num for chains if nbreaks > 0)
     */
    if (nres < 0) break;
    if (nres == 0) continue;

    newchain = 1;
    if (modelid >= 0) {		/* start of a new MODEL block */
      topofblock = pe.nchains;
      i = pe.nchains;
      mid = modelid;
      nblocks++;
    } else {			/* no beginning MODEL record */
      for (i = topofblock; i < pe.nchains; i++)	// don't look back prev models
	if (k == pe.chainptr[i].chainid) {
	  newchain = 0;
	  break;
	}
    }
    if (!newchain) {
      /*
      if (ishetatm)
	;	// ignore
      else {
      */
	if (resid_max > pe.chainptr[i].resseq_max)
	  pe.chainptr[i].resseq_max = resid_max;
	if (resid_min < pe.chainptr[i].resseq_min)
	  pe.chainptr[i].resseq_min = resid_min;
	/*
      }
	*/
    } else {
      pe.nchains++;
      pe.chainptr[i].index = i;
      pe.chainptr[i].chainid = k;
      pe.chainptr[i].model = mid;
      pe.chainptr[i].entry = 0;

      if ((nblocks == 1 || nblocks == 0))
	nchains_block++;

      /*
      if (ishetatm)
	;		// ignore
      else {
      */
	pe.chainptr[i].resseq_max = resid_max;
	pe.chainptr[i].resseq_min = resid_min;
	/*
      }
	*/
    }

    qsort(buf, strlen(buf), sizeof(char),
	  (int (*)(const void *, const void *))strcoll);
    merge(pe.chainptr[i].altlocs, buf);

    if (nblocks == 1 || nblocks == 0) {	/* first MODEL or no MODELs */
      chain_block[pe.nchains - 1] = pe.chainptr[i].chainid;
    }

    if (debugprint) {
      printf("REMARK| %s chain/model-idx=%d, mid = %d,"
	     "cid = %d, resid_max=%d, altlocs=\"%s\"\n",
	     ishetatm ? "[HET]" : "", i, pe.chainptr[i].model,
	     pe.chainptr[i].chainid, resid_max, pe.chainptr[i].altlocs);
    }
  }
  resid_max = pe.chainptr[0].resseq_max;
  resid_min = pe.chainptr[0].resseq_min;
  for (i = 1; i < pe.nchains; i++) {
    if (pe.chainptr[i].resseq_max > resid_max)
      resid_max = pe.chainptr[i].resseq_max;
    if (pe.chainptr[i].resseq_min < resid_min)
      resid_min = pe.chainptr[i].resseq_min;
  }
  if (debugprint) {
    printf("REMARK| nchains = %d\n", pe.nchains);
    printf("REMARK| resid_max = %d\n", resid_max);
    printf("REMARK| resid_min = %d\n", resid_min);
  }
  
  /*
   * preparing a LUT for lardmark ids.
   */
  // LUT = landmark_initlut(resid_min, resid_max, nlandmarks, lm_pool);
  LUT = landmark_initlut(resid_min, resid_max, nlandmarks);
  rangemarked = landmark_preplut(resid_max, nlandmarks, lm_pool, LUT);

  if (debugprint) {
    printf("REMARK| the LUT for lardmark ids.\n");
    printf("REMARK| %6s ", " ");
    for (j = resid_min; j <= resid_max; j++)
      printf(" %3d", j);
    printf("\n");
    for (i = 0; i < nlandmarks; i++) {
      printf("REMARK| %6s ", lm_pool[i].name);
      for (j = resid_min; j <= resid_max; j++)
	printf(" %03d", LUT[i][j]);
      printf("\n");
    }
    printf("REMARK| range marked [%d,%d]\n",
	   rangemarked.lower,rangemarked.upper);
  }

  if (arguments.fasta) {
    fasta = fopen(arguments.fasta, "r");
    i = scan_fasta(fasta, 0, 0, 0, 0, -1);
    fclose(fasta);

    fasta_raise = malloc(sizeof(int) * (resid_max + i + 1));
    if (fasta_raise == 0) {
      fprintf(stderr, "cannot alloc fast_raise[].\n");
      exit(-1);
    }
    fasta = fopen(arguments.fasta, "r");
    i = scan_fasta(fasta, arguments.fasta_off, pe.file, fasta_seq, fasta_raise,
		   NRES_LIM);
    fclose(fasta);
    for (j = i; j <= resid_max; j++)
      fasta_raise[j] = fasta_raise[i];

    for (i = 0; i < nlandmarks; i++)
      for (j = resid_min; j <= resid_max; j++)
	if (LUT[i][j] != 0)
	  LUT[i][j] += fasta_raise[j];

    if (debugprint) {
      printf("REMARK| %6s ", "LUT raised (LM1)");
      for (j = resid_min; j <= resid_max; j++)
	printf(" %03d", LUT[0][j]);
      printf("\n");
    }
  }

  /*
   * rescan the pdb file and print each line appending id # if landmarks.
   * - records to be printed are ATOM, HETATM, MODEL, END, ENDMDL, & TER.
   *   others are ignored.
   */

  op = 0;
  if (arguments.nbreaks < 0) {		/* merge into a single chain */
    op = offset;
    memset(offset, 0, sizeof(int)*256);
    offset[chain_block[0]] = 0;
    for (i = 1; i < nchains_block; i++)
      offset[chain_block[i]] = offset[chain_block[i-1]]
	+ rangemarked.upper - rangemarked.lower + 1;
  }

  rewind(fp);
  modelid = -1;
  while (1) {
    nres = pdb_scanmark(fp, buf /*altloc*/, &k /*chainid*/, &modelid, &ishetatm,
			&resid_min, &resid_max, &pe, LUT, op,
			arguments.nbreaks > 0 ? arguments.breakpoints : 0);
    if (nres < 0) break;
  }
  fclose(fp);

}