#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atypes.h"
#include "pdb.h"

static double read_double(char rec[], int len)
{
  double d;
  char buf[BUFSIZ];
  strncpy(buf, rec, len); buf[len] = 0;
  sscanf(buf, "%lf", &d);
  return d;
}

static int read_int(char rec[], int len)
{
  int i;
  char buf[BUFSIZ];
  strncpy(buf, rec, len); buf[len] = 0;
  sscanf(buf, "%d", &i);
  return i;
}

static void read_string(char *rec, int len, char s[])
{
  for (; len > 0; rec++,len--)
    if (!isspace(*rec)) break;
  strncpy(s, rec, len);
  for (; len > 0; len--)
    if (!isspace(s[len-1])) break;
  s[len] = 0;
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


/*
 * read_pdbm() reads a contiguous segment of a chain
 * from the 'landmarked' PDB file, `fp`.
 * a 'landmarked' record (either ATOM or HETATM) contains an id following @
 * after the 80-th column.
 * On output,
 * *modelp == MODEL serial num (> 0) if specified, or == -1 otherwise.
 * read_pdbm() returns PDB chain ID.
 */
static char prefetchbuf[BUFSIZ], *prefetched = 0;

int
read_pdbm(FILE *fp,
          int id2ind[],		// id2ind[landmark id] == index
          int id_lim,		// upper limit for landmark ids
          pdbAtom atom[], pdbSite site[],
          int *id_minp, int *id_maxp, int *natomsp, int *nsitesp, int *modelp)
{
  int i, k, chain, resid, len, ind, model, maxid, minid, resseq_max;
  int eof, hetatm;
  int natoms, nsites, nres;
  char rec[BUFSIZ], c;
  pdbSite *sp;

  model = 0;

  for (i = 0; i <= id_lim; i++)
    id2ind[i] = -1;

  eof = hetatm = 0;
  natoms = nsites = nres = 0;
  maxid = 0; minid = INT_MAX;
  chain = 0;			/* pdb chainID as specified in *fp */

  while (1) {
    if (prefetched) {
      strcpy(rec, prefetched);
      prefetched = 0;
    } else if (!fgets(rec, BUFSIZ, fp)) {
      eof = 1;  break;
    }

    if (pdbrec("END", rec) || pdbrec("ENDMDL", rec)) {
      break;
    } else if (pdbrec("TER", rec)) {
      break;
    } else if (pdbrec("MODEL", rec)) {
      model = read_int(rec+PDB_MODEL_serial_LOC, PDB_MODEL_serial_LEN);
      model = -model;
    } else if (pdbrec("ATOM", rec) || pdbrec("HETATM", rec)) {
      if (chain == 0 || model < 0) {			// start a new chain
        chain = rec[PDB_ATOM_chainID_LOC];
	if (model < 0) model = -model;
        hetatm = rec[0] == 'H';
        resid = -1;		// residue ID (resSeq)
        resseq_max = 0;
        nres = 0;		// num of residues
      } else if (chain != rec[PDB_ATOM_chainID_LOC]	// next chain/fragm, or
		 || (rec[0] == 'H') != hetatm) {	// ATOM <-> HETATM,
	strcpy(prefetchbuf, rec);		// push bask
	prefetched = prefetchbuf;
	break;
      }

      k = read_int(rec+PDB_ATOM_resSeq_LOC, PDB_ATOM_resSeq_LEN);
      if (k != resid) {		// new residue ?
        nres++;
        resid = k;
        if (resid > resseq_max)
          resseq_max = resid;
      }

      len = strlen(rec);
      rec[--len] = 0;
      if (len < 80) continue;
      for (i = 80; rec[i]; i++)	// skip white spaces
        if (!isspace(rec[i])) break;

      if (i >= 80 && rec[i] == '@') {		// extra field for landmark id

        sscanf(rec+i+1, "%d", &k);		// k is the landmark id
        if (k <= 0) {
          fprintf(stderr, "id = %d < 0.\n", k);
          return -1;
        } else if (k >= id_lim) {
          fprintf(stderr, "id = %d >= IDNUM_LIM(%d).\n", k, id_lim);
          return -1;
        }
        if (k > maxid) maxid = k;
        if (k < minid) minid = k;
        
        c = rec[PDB_ATOM_altLoc_LOC];
        ind = id2ind[k];
        if (ind < 0) {			// 必ず、リストの末尾に加える。【要検討】
          ind = id2ind[k] = natoms++;
          atom[ind].id = k;
          atom[ind].sp = site + nsites++;
          atom[ind].sp->next = 0;
          sp = atom[ind].sp;
          atom[ind].naltlocs = 0;
        } else {
          sp = atom[ind].sp;
          while (sp->next && sp->altloc != c)
            sp = sp->next;
          if (sp->next) {
            fprintf(stdout, "duplicated site %c for atom id = %d\n", c, k);
            return -1;
          }
          sp = sp->next = site + nsites++;
          sp->next = 0;
          atom[ind].naltlocs++;
        }
        /* sp points a new site */
        sp->altloc = c;
        sp->occ = read_double(rec+PDB_ATOM_occ_LOC, PDB_ATOM_occ_LEN);
	sp->Biso = read_double(rec+PDB_ATOM_tfac_LOC, PDB_ATOM_tfac_LEN);
        /* p, name, resname, resseq are unnecessary on prescanning */
	if (strlen(rec) > PDBM_EXTRA_z_LOC) {
	  sp->p.x = read_double(rec+PDBM_EXTRA_x_LOC, PDBM_EXTRA_x_LEN);
	  sp->p.y = read_double(rec+PDBM_EXTRA_y_LOC, PDBM_EXTRA_y_LEN);
	  sp->p.z = read_double(rec+PDBM_EXTRA_z_LOC, PDBM_EXTRA_z_LEN);
	} else {
	  sp->p.x = read_double(rec+PDB_ATOM_x_LOC, PDB_ATOM_x_LEN);
	  sp->p.y = read_double(rec+PDB_ATOM_y_LOC, PDB_ATOM_y_LEN);
	  sp->p.z = read_double(rec+PDB_ATOM_z_LOC, PDB_ATOM_z_LEN);
	}
        read_string(rec+PDB_ATOM_name_LOC, PDB_ATOM_name_LEN, atom[ind].name);
        read_string(rec+PDB_ATOM_resName_LOC, PDB_ATOM_resName_LEN,
                    atom[ind].resname);
        atom[ind].resseq = resid;
      }         /* if (i >= 80 && rec[i] == '@') */
    }   /* if (ATOM or HEATOM) */
  }     /* while (1) */

  *id_maxp = maxid;
  *id_minp = minid;
  *natomsp = natoms;
  *nsitesp = nsites;
  if (model == 0) model = -1;
  *modelp = model;

  if (eof && maxid == 0)
    return EOF;
  else
    return chain;
}


/*
 * notice read_trans_write_pdbm() processes a pdbm file throughtly at once,
 * while read_pdbm() processes a fragment in a pdbm file at once.
 */

int
read_trans_write_pdbm(FILE *out, pdbEntry *pe,
                      Real *Rbuf, Real *tbuf,
                      int coordsonly)
{
  Matrix3 matR;
  Real *vect, x[3], sx[3];
  int chain, model, k, fwidth, fprec; //, cid
  int nchains, top_block;
  char rec[BUFSIZ], bufxyz[BUFSIZ];
  FILE *in;
  extern int global_chain_id(int ient, int chid, int mdlid);

  in = fopen(pe->file, "r");

  fwidth = PDBM_EXTRA_x_LEN - 1;
  fprec = fwidth - 7;

  matR = NULL; vect = NULL;
  model = 0;
  chain = 0;
  nchains = 0;
  top_block = nchains;
  if (pe->nchains_block > 0)		// pe is consisted of MODEL blocks
    top_block = -(pe->nchains_block);
  while (fgets(rec, BUFSIZ, in)) {
    if (pdbrec("END", rec)) {
      // fprintf(out, "%s", rec);	// END should be put at the very end
      break;
    } else if (pdbrec("ATOM", rec) || pdbrec("HETATM", rec)) {
      if (chain != rec[PDB_ATOM_chainID_LOC]	/* new chain or fragment */
	  || model < 0) {    
        chain = rec[PDB_ATOM_chainID_LOC];
        if (model < 0) { model = -model; top_block += pe->nchains_block; }
        // cid = global_chain_id(pe->id, chain, model);
        for (k = top_block; k < pe->nchains; k++)        // search for the chain
          if (chain == pe->chainptr[k].chainid
	      && (pe->chainptr[k].model < 0 || model == pe->chainptr[k].model))
            break;
        if (k >= pe->nchains) {         // unknown chain
	  fprintf(stderr, "chain (id=%d,model=%d) not found\n", chain, model);
	  exit(-1);
        } else {
          /* read_coordinates() では、archep = archepool + pe->chain[k]->index
             としているが、bufR/tについても同様の保証がある？ */
          matR = (Matrix3)(Rbuf + 9 * pe->chainptr[k].index);
          vect = tbuf + 3 * pe->chainptr[k].index;
        }
      }
      x[0] = read_double(rec + PDB_ATOM_x_LOC, PDB_ATOM_x_LEN);
      x[1] = read_double(rec + PDB_ATOM_y_LOC, PDB_ATOM_y_LEN);
      x[2] = read_double(rec + PDB_ATOM_z_LOC, PDB_ATOM_z_LEN);
      sx[0] = x[0] - vect[0];
      sx[1] = x[1] - vect[1];
      sx[2] = x[2] - vect[2];
      x[0] = sx[0]*matR[0][0] + sx[1]*matR[1][0] + sx[2]*matR[2][0];
      x[1] = sx[0]*matR[0][1] + sx[1]*matR[1][1] + sx[2]*matR[2][1];
      x[2] = sx[0]*matR[0][2] + sx[1]*matR[1][2] + sx[2]*matR[2][2];
      sprintf(bufxyz, "%8.3f%8.3f%8.3f", x[0], x[1], x[2]);
      memcpy(rec+PDB_ATOM_x_LOC, bufxyz, 24);
      if (rec[PDBM_EXTRA_FLAG_LOC] == '@' && strlen(rec) > PDBM_EXTRA_z_LOC) {
	sprintf(bufxyz, " %*.*e %*.*e %*.*e", fwidth, fprec, x[0],
		fwidth, fprec, x[1], fwidth, fprec, x[2]);
	memcpy(rec+PDBM_EXTRA_x_LOC, bufxyz, (fwidth+1)*3);
      }
    } else if (pdbrec("MODEL", rec)) {
      model = read_int(rec+PDB_MODEL_serial_LOC, PDB_MODEL_serial_LEN);
      model = -model;
    } else if (coordsonly
               && !pdbrec("TER", rec) && !pdbrec("ENDMDL", rec))
      continue;
    fprintf(out, "%s", rec);
  }
  fclose(in);
  return 0;
}
