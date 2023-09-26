#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
       #include <sys/types.h>
       #include <sys/stat.h>
       #include <fcntl.h>
#include <unistd.h>
#include "atypes.h"
#include "matrix.h"
#include "memory.h"
#include "molec.h"
#include "pdb.h"
#include "mmfit.h"
#include "defs.h"


#ifdef MKLROOT
int UseMKL = 1;
#else
int UseMKL = 0;
#endif

#define PDB_NENT_MAX 10000
#define IDNUM_LIM 10000			/* upper limit of landmark ID nums  */
#define CHAINPOOL_SIZ (1<<18)
#define OCC_EPS (1.0e-6)		/* to cope with rounding errors     */

#define WORKBUFSIZ (1<<24)

int id2ind[IDNUM_LIM+2];		/* 0..IDNUM_LIM and a sentinel      */
int idcnt[IDNUM_LIM+2];
int idcnttmp[IDNUM_LIM+2];
pdbEntry pdbentry[PDB_NENT_MAX];
pdbChain *pdbchainpool = 0;		/* storage for pdbChains */

#define NSITES_MAX 10000		/* max # of sites(atoms) in an entry */
pdbAtom pdb_atom[NSITES_MAX];		/* reuse pdb_atom and pdb_site       */
pdbSite pdb_site[NSITES_MAX];		/*   for all entries                 */

/*
 * prescan a pdbm file for obtaining
 *  - min & max values of landmark ids (to prepare id2ind[])
 */
int
prescan_pdbm(FILE *fp, pdbEntry *pe, int idcnt[], int id_lim,
	     pdbChain *chainbuf, int nchains_lim, int *minidp, int *maxidp,
	     int *work)
{
  int frag_id_min, frag_id_max, minid, maxid;
  int natoms, nsites, nchains, model_id, cid;
  int i, k, ichain, iatom, *id2ind;
  int top_block;
  pdbChain *chainp;

  id2ind = work;
  minid = *minidp; maxid = *maxidp;
  pe->nchains = nchains = 0;
  model_id = 0;
  pe->nchains_block = -1;	// indicating no MODEL record

  top_block = nchains;

#ifdef DEBUG
  printf("prescaning %s...\n", pe->file);
#endif

  while ((k = read_pdbm(fp, id2ind, id_lim, pdb_atom, pdb_site,
			&frag_id_min, &frag_id_max,
			&natoms, &nsites, &model_id)) != EOF) {
    // cid = global_chain_id(pe->id, k, model_id);
    cid = k;		// k is the chain ID field of the ATOM records
    
#ifdef DEBUG
    printf("|   cid=%d, natoms/nsites=%d/%d", cid, natoms, nsites);
    if (frag_id_max > 0)
      printf(", id=%d:%d  ", frag_id_min, frag_id_max);
    else
      printf("\n");
#endif

    if (frag_id_max == 0)		// skip if this frag contains no id
      continue;
    
    if (frag_id_max > maxid) maxid = frag_id_max;
    if (frag_id_min < minid) minid = frag_id_min;
    for (i = frag_id_min; i <= frag_id_max; i++)
      if (id2ind[i] >= 0)
	idcnt[i]++;
    if (model_id > 0) {			// start of a new MODEL block
      top_block = nchains;
      pe->nchains_block = 1;		// indicating any MODEL records
    }
    
    for (k = top_block; k < pe->nchains; k++)		// search for the chain
      if (cid == pe->chainptr[k].chainid) 
	break;
    if (k < pe->nchains) {			// found in the list
      ichain = k;
      chainp = pe->chainptr + ichain;
    } else {					// a new chain
      ichain = pe->nchains;
      if (ichain == 0) pe->chainptr = chainbuf;
      chainp = pe->chainptr + ichain;
      pe->nchains++;
      nchains++;
      if (nchains > nchains_lim) {
	fprintf(stderr, "too many chains.\n"); return -1;
      }
      chainp->chainid = cid;
      chainp->model = model_id;
      chainp->entry = pe->id;
      //chainp->occ_nparts = 0;	/* num of partitionings for occ on an atom */
      //chainp->occ_nfracs[chainp->occ_nparts] = 0;
      chainp->natoms = 0;
    }

    chainp->natoms += natoms;
    
#ifdef DEBUG
    if (natoms > 0) printf("ichn=%d\n", ichain);
#endif

    /*
     *
     */
    for (iatom = 0; iatom < natoms; iatom++) {
      if (pdb_atom[iatom].sp->occ != 1.0)
	fprintf(stderr, "Warning: fractional occ. (%s %d %d %d)\n",
		pe->file, cid, model_id, iatom);
    }
  }		/* while ((k = read_pdbm(...))) */

  if (pe->nchains_block > 0)
    pe->nchains_block = pe->nchains - top_block;
  else
    pe->nchains_block = 0;
  *minidp = minid; *maxidp = maxid;

  return 0;
}

/*
 * sum up anisotropic cov components to obtain the isotropic cov.
 * assuming anisotropic cov == isotropic cov \otimes I3
 */
void gather_aniso(Real *dest, Real *src, int natoms, int diag)
{
  int i, j, n3;
  Real *wkbuf;

  n3 = natoms * 3;
  if (dest == src) {
    wkbuf = get_workarea("aniso", natoms*natoms, sizeof(Real));
    if (wkbuf == 0)
      abnormalexit("gather_aniso() cannot allocate wkbuf", 0);
  } else
    wkbuf = dest;

  if (diag) {
    for (i = 0; i < natoms; i++)
      wkbuf[i] = (src[i*3] + src[i*3+1] + src[i*3+2]);	// /  3.0;
    if (wkbuf != dest)
      for (i = 0; i < natoms; i++)
	dest[i] = wkbuf[i];
  } else {
    for (i = 0; i < natoms; i++)
      for (j = 0; j < natoms; j++)
        wkbuf[i*natoms+j]
          = (src[(i*3)*n3+(j*3)] + src[(i*3+1)*n3+(j*3+1)]
             + src[(i*3+2)*n3+(j*3+2)]); //	/  3.0;
    if (wkbuf != dest)
      for (i = 0; i < natoms * natoms; i++)
	dest[i] = wkbuf[i];
  }
  if (wkbuf != dest)
    remove_workarea("aniso");
}

/*
 * scatter the isotropic cov components over the anisotropic matrix
 * assuming anisotropic cov == isotropic cov \otimes I3
 */
void scatter_aniso(Real *dest, Real *src, int natoms, int diag)
{
  int i, j, n3;
  Real *wkbuf;

  n3 = natoms * 3;
  if (dest == src) {
    wkbuf = get_workarea("aniso", n3*n3, sizeof(Real));
    if (wkbuf == 0)
      abnormalexit("scatter_aniso() cannot allocate wkbuf", 0);
  } else
    wkbuf = dest;
  if (diag) {
    for (i = 0; i < natoms; i++)
      wkbuf[i*3] = wkbuf[i*3+1] = wkbuf[i*3+2] = src[i];
    if (wkbuf != dest)
      for (i = 0; i < n3; i++)
	dest[i] = wkbuf[i];
  } else {
    for (i = 0; i < natoms; i++)
      for (j = 0; j < natoms; j++)
        wkbuf[(i*3)*n3+(j*3)] = wkbuf[(i*3+1)*n3+(j*3+1)]
          = wkbuf[(i*3+2)*n3+(j*3+2)] = src[i*natoms+j];
    if (wkbuf != dest)
      for (i = 0; i < n3 * n3; i++)
	dest[i] = wkbuf[i];
  }
  if (wkbuf != dest)
    remove_workarea("aniso");
}




#define NGROUPS_LIM 256

Atom *Abuf, *chain_template;
Conformer *Cbuf;
int *nconfms;
Ensemble *Bbuf;
struct jvar *Jbuf;
Real *Xbuf;
Real *Ybuf;
int *Obuf;
Real *Mbuf;
Real *tbuf;
Real *Rbuf;
Real *Zbuf;
Real *Ubuf;
Real *Sbuf;
Real *Vbuf;

int
mmfit_init(int natoms, int ngroups, int nconfms_tot, pdbEntry *pe, int method,
	   int *njvars)
{
  Real *xp, *yp, *zp;
  int i, j, m, iconfm, iatom, iens, *op, *mi, *ix, *ixm, nmi;
  Conformer *conform;
  Atom *ap;

  nconfms = (int *) malloc(sizeof(int) * (ngroups * 4 + 1));
  if (!nconfms) {
    fprintf(stderr, "mmfit_init(): cannot alloc memory.\n");
    return -1;
  }
  // U_i and some workarea are shared between ensembles of the same sizes
				// mi[-1] == num of different mi's, nmi
  mi = nconfms + ngroups + 1;	// mi[] <- unique(sort(nconfms[]))
  ix = mi + ngroups;		// ix[] s.t. sort(nconfms)[k] == nconfms[ix[k]]
  ixm = ix + ngroups;		// ixm[] s.t. nconfoms[k] == mi[ixm[k]]

  for (iens = 0; iens < ngroups; iens++) {
    nconfms[iens] = pe[iens].nchains;
    mi[iens] = nconfms[iens];
    ix[iens] = iens;
  }
  for (i = 0; i < ngroups; i++)
    for (j = i+1; j < ngroups; j++)
      if (mi[j] < mi[i]) {
        m = mi[j];
        mi[j] = mi[i];
        mi[i] = m;
        m = ix[j];
        ix[j] = ix[i];
        ix[i] = m;
      }
  j = 0;
  for (i = 0; i < ngroups; i++) {
    if (mi[i] > mi[j])
      mi[++j] = mi[i];
    ixm[ix[i]] = j;
  }
  nmi = j + 1;		// num of different mi's (== # of S+mi*V)
  nconfms[ngroups] = nmi;

  Abuf = (Atom *) malloc(sizeof(Atom) * natoms * (nconfms_tot + 1));
  chain_template = Abuf + natoms * nconfms_tot;
  Bbuf = (Ensemble *) malloc(sizeof(Ensemble) * ngroups);
  Cbuf = (Conformer *) malloc(sizeof(Conformer) * nconfms_tot);
  Mbuf = (Real *) malloc(sizeof(Real) * 3 * natoms * 3); /* curr, prev, excl */
	/* occ + nconfms[] */
  Obuf = (int *) malloc(sizeof(int) * (natoms * nconfms_tot));
	/* twice for current and previous values, and once for exclusive-ave */
  tbuf = (Real *) malloc(sizeof(Real) * 3 * nconfms_tot  * (2+1));
  Rbuf = (Real *) malloc(sizeof(Real) * 3 * 3 * nconfms_tot  * (2+1));
						/* twice for X and Y */
  Xbuf = (Real *) malloc(sizeof(Real) * 3 * natoms * nconfms_tot * 2);
  Ybuf = Xbuf + 3 * natoms * nconfms_tot;

  Ubuf = (Real *) malloc(sizeof(Real) * natoms * natoms * nmi * NUMJVMATRICES );
  Jbuf = (struct jvar*) malloc(sizeof(struct jvar) * nmi);
  Vbuf = (Real *) malloc(sizeof(Real)	/* full curr/prev, diag c/p, aniso */
			 * (natoms * natoms * 11 + natoms * 2));
  Sbuf = (Real *) malloc(sizeof(Real) * (natoms * natoms * 11 + natoms * 2));
  Zbuf = (Real *) malloc(sizeof(Real) * 3 * natoms * ngroups * 2);

  if (!Abuf || !Cbuf || !tbuf || !Rbuf || !Mbuf || !Xbuf
      || !Vbuf || !Sbuf || !Ubuf || !Jbuf || !Zbuf) {
    fprintf(stderr, "mmfit_init(): cannot alloc memory.\n");
    return -1;
  }
  

  // これで Atom.id も Conform.occ も 0 に初期化される
  memset(Abuf, 0, sizeof(Atom) * natoms * (nconfms_tot + 1));
  memset(Bbuf, 0, sizeof(Ensemble) * ngroups);
  memset(Cbuf, 0, sizeof(Conformer) * nconfms_tot);
  memset(Obuf, 0, sizeof(int) * natoms * nconfms_tot);

  *njvars = nmi;
  for (i = 0; i < nmi; i++) {
    Jbuf[i].mi = mi[i];
    xp = Ubuf + i * NUMJVMATRICES * natoms*natoms;
    Jbuf[i].u    = xp;
    Jbuf[i].sv   = xp + 1 * natoms*natoms; 
    Jbuf[i].b = xp + 3 * natoms*natoms; 
    if ((method & METHOD_FIX) == METHOD_FIX_ML
	|| (method & METHOD_VAR) == METHOD_VAR_REML)
      Jbuf[i].svi  = xp + 2 * natoms*natoms; 
    else
      Jbuf[i].svi  = 0;
    if ((method & METHOD_FIX) == METHOD_FIX_ML) {
      Jbuf[i].c    = xp + 4 * natoms*natoms; 
      Jbuf[i].sc   = xp + 5 * natoms*natoms; 
    } else {
      Jbuf[i].c    = 0;
      Jbuf[i].sc   = 0;
    }
    if (((method & METHOD_FIX) == METHOD_FIX_ML && (nmi > 1))
	|| (method & METHOD_VAR) == METHOD_VAR_REML)
      Jbuf[i].h    = xp + 6 * natoms*natoms; 
    else
      Jbuf[i].h    = 0;
    Jbuf[i].nens = 0;
    if (i < nmi-1)
      Jbuf[i].next = &Jbuf[i+1];
    else
      Jbuf[i].next = NULL;
  }

  xp = Xbuf;
  yp = Ybuf;
  zp = Zbuf;
  conform = Cbuf;
  for (iens = 0; iens < ngroups; iens++) {
    Bbuf[iens].id = iens;
    Bbuf[iens].x = xp;
    Bbuf[iens].y = yp;
    Bbuf[iens].z = zp;
    Bbuf[iens].jvp = &Jbuf[ixm[iens]];
    Bbuf[iens].jvp->nens += 1;
    Bbuf[iens].nconfms = nconfms[iens];
    Bbuf[iens].confmbp = conform;
    xp += 3 * natoms * Bbuf[iens].nconfms;
    yp += 3 * natoms * Bbuf[iens].nconfms;
    zp += 3 * natoms;
    conform += Bbuf[iens].nconfms;
  }
    
  ap = Abuf;
  xp = Xbuf;
  yp = Ybuf;
  op = Obuf;
  iconfm = 0;
  conform = Cbuf;
  for (iconfm = 0; iconfm < nconfms_tot; iconfm++) {
    conform[iconfm].index = iconfm;
    conform[iconfm].head = ap;
    conform[iconfm].natoms = natoms;
    conform[iconfm].x = xp;
    conform[iconfm].y = yp;
    conform[iconfm].o = op;
    conform[iconfm].t = tbuf + (3 * iconfm);
    conform[iconfm].r = Rbuf + (9 * iconfm);
    for (iatom = 0; iatom < natoms; iatom++) {
      ap->p = xp;
      ap->occ = op;
      ap++; op++;
      xp += 3;
    }
    yp += 3 * natoms;
  }

  return 0;
}  


static char *progname = 0;
static char usagemessage[] =
  " Usage\n"
  "  mmfit [Options] pdbm_1 pdbm_2 ... pdbm_r\n"
  "  mmfit [Options] -i listfile\n"
  "    listfile: a file containing paths of pdbm files, one per line.\n" 
  "  Option:\n"
  "  -e epsilon:  used in stopping criteria. default epsilon = epsM^(1/2)\n"
  "  -n niters: upper lim of num of iterations\n"
  "  -i file:  read a list of pdbm files from file\n"
  "  -v [vsvs.dat]: output V & S by est1 and by est2 onto vsvs.dat\n"
  "  -c [m,crd]: output estimated M & Yij onto m.pdb & crd.pdb\n"
  "  -f est1: estimator for fix effect: 0= max Q, 1= max marL\n"
  "  -r est2: estimator for random effect: 0= ML, 1= REML\n"
  "  -a align estimates M on its principal axes.\n"
  "  -d ignore off-diagonal V on IGLS\n"
  "  -s spherical (isotropic) Sigma\n"
  "  -I level: impute coordinates for imcomplete cases\n"
  "         if the atoms are prezented in more than 100*level % chains.\n"
  "  -B balanced design (assuming mi == m/r).\n"
;


/*
 * binary I/O:
 * -b i crdext
 *   overwrite Xbuf with the contents of pefile.crdext files,
 *   where pefile is a basename of the pdbm files.
 * -b o binout
 *   output M, Yij, Zi, V, S [, V, S] to
 *   the file 'd'binout and natoms, nmi, mi's to the file 'i'binout.
 */
char *binaryIO_crdext = 0;
char *binaryIO_binout = 0;
char ibinout[BUFSIZ], dbinout[BUFSIZ];


void usage()
{
  printf("%s", usagemessage);
  exit(-1);
}

void abnormalexit(char *message, int code)
{
  if (message && code)
    fprintf(stderr, "(%d) %s\n", code, message);
  else if (message)
    fprintf(stderr, "%s\n", message);
  exit(code);
}

struct options {
  Real tol;
  int niters;
  int alignM;			// removed
  int impute;
  int isoVar;
  Real impute_lowerlim;
  int debug;
  char *covout;
  char *mpdb, *ypdb;
  int diagV;
  int fix;
  int ran;
  int YinM;			// to be removed
  int balanced;
} options;



int
parse_args(int argc, char *argv[], struct options *options, pdbEntry *pe)
{
  int opt, i, ngroups;
  Real machineps;
  FILE *fp;
  char *cp, *pdbmlist, buf[BUFSIZ];
  char default_covout[] = "vsvs.dat";
  char default_mpdb[] = "m.pdb";
  char default_ypdb[] = "crd.pdb";

  options->isoVar = 0;		/* NOT isotropic (spherical) variance */
  options->diagV = 0;
  options->YinM = 0;
  options->debug = 0;
  options->impute = 0;
  options->covout = 0;
  options->mpdb = options->ypdb = 0;
  options->fix = 0;
  options->ran = 0;
  options->balanced = 0;
  pdbmlist = 0;

  machineps = DBL_EPSILON;
  options->tol = pow(machineps, 1.0 / 2.0);
  options->niters = 5000;			  /* limit for niters */

  while ((opt = getopt(argc, argv, "cde:f:i:n:pvr:sDI:Eb:B")) != -1) 
    switch (opt) {
    case 'e':
      options->tol = atof(optarg);
      break;
    case 'i':
      pdbmlist = strdup(optarg);
      break;
    case 'n':
      options->niters = atoi(optarg);
      break;
    case 's':
      options->isoVar = 1;		// spherical (isotropic) variance
      break;
    case 'D':
      options->debug = 1;
      break;
    case 'I':
      options->impute = 1;
      options->impute_lowerlim = atof(optarg);
      break;
    case 'v':
      if (argv[optind][0] == '-') {
	options->covout = strdup(default_covout);
      } else {
	options->covout = strdup(argv[optind++]);
      }
      break;
    case 'c':
      if (argv[optind][0] == '-') {		// filenames omitted
	options->mpdb = strdup(default_mpdb);
	options->ypdb = strdup(default_ypdb);
      } else {
	cp = strtok(argv[optind], ",");
	if (cp == NULL || strlen(cp) == 0) usage();
	options->mpdb = strcat(strcpy(malloc(strlen(cp)+5), cp), ".pdb");
	cp = strtok(NULL, ",");
	if (cp == NULL || strlen(cp) == 0) usage();
	options->ypdb = strcat(strcpy(malloc(strlen(cp)+5), cp), ".pdb");
	cp = strtok(NULL, ".");
	if (cp != NULL) usage();
	optind++;
      }
      break;
    case 'd':
      options->diagV = 1;
      break;
    case 'f':
      options->fix = atoi(optarg);
      break;
    case 'r':
      options->ran = atoi(optarg);
      break;
    case 'B':
      options->balanced = 1;
      break;
    case 'b':
      switch (*optarg) {
      case 'i':
	binaryIO_crdext = strdup(argv[optind++]);
	break;
      case 'o':
	binaryIO_binout = strdup(argv[optind++]);
	sprintf(ibinout, "i%s", binaryIO_binout);
	sprintf(dbinout, "d%s", binaryIO_binout);
	break;
      }
      break;
    default:
      usage();
    }
  if (pdbmlist) {
    if (argc != optind)
      usage();
    ngroups = 0;
    fp = fopen(pdbmlist, "r");
    while (fgets(buf, BUFSIZ, fp)) {
      for (cp = buf; isspace(*cp); cp++)
	;
      cp[strlen(cp)] = '\0';
      pe->file = strdup(cp);
      for (cp = pe->file+strlen(cp)-1; isspace(*cp); cp--)
	;
      *++cp = '\0';
      ngroups++;
      pe++;
    }
    fclose(fp);
  } else {
    ngroups = argc - optind;
    if (ngroups < 2)
      usage();
    for (i = 0; i < ngroups; i++, pe++)
      pe->file = argv[optind + i];
  }

  return ngroups;
}


/*
 * work[0..id_lim] 
 */
int
read_coordinates(FILE *fp, pdbEntry *pe, int id2ind[], int id_lim,
		 int natoms_all, int work[])
{
  int frag_id_min, frag_id_max;
  int natoms, nsites, model_id, cid;
  int k, ichain, iatom, ind, n;
  int nchains, top_block, found;
  Real d;
  pdbSite *sp;
  Conformer *confmp;
  Atom *ap;

  nchains = 0;
  top_block = nchains;
  if (pe->nchains_block > 0)		// pe is consisted of MODEL blocks
    top_block = -(pe->nchains_block);	// of a constant num of chains
  // model_id = 0;		// ???? model_id = -1; ???
  while ((cid = read_pdbm(fp, work, id_lim, pdb_atom, pdb_site,
                          &frag_id_min, &frag_id_max,
                          &natoms, &nsites, &model_id)) != EOF) {
    // cid = global_chain_id(pe->id, cid, model_id);
    if (frag_id_max == 0) continue;	// skip if this frag contains no id
    if (model_id > 0)
      top_block += pe->nchains_block;
    found = 0; // search for the chain
    if (pe->nchains_block > 0) {	// contains MODELs (NMR/MD)
      for (k = top_block; k < top_block + pe->nchains_block; k++) 
	if (cid == pe->chainptr[k].chainid) {
	  found = 1;
	  break;
	}
    } else {
      for (k = 0; k < pe->nchains; k++) {
	if (cid == pe->chainptr[k].chainid) {
	  found = 1;
	  break;
	}
      }
      printf("\n");
    }

    if (!found) {
      fprintf(stderr, "cid %x not found in entry %s.\n", cid, pe->file);
      return -1;
    }

    ichain = k;
    confmp = Cbuf + pe->chainptr[ichain].index;
    for (iatom = 0; iatom < natoms; iatom++) {
      sp = pdb_atom[iatom].sp;
      d = sp->occ;
      ind = id2ind[pdb_atom[iatom].id];
      if (ind < 0) {
	fprintf(stderr, "Warning: negative index %d for atom #%d (id %d) in %s."
		" ignored\n", ind, iatom, pdb_atom[iatom].id, pe->file);
	continue;
      }

      {
	ap = confmp->head + ind;
	if (d > 0.0) {
	  *ap->occ = 1;
	  // ap->disordered = pdb_atom[iatom].naltlocs > 0;
	  ap->p[0] = sp->p.x;
	  ap->p[1] = sp->p.y;
	  ap->p[2] = sp->p.z;
	  ap->biso = sp->Biso;
	  strcpy(ap->name, pdb_atom[iatom].name);
	  strcpy(ap->resname, pdb_atom[iatom].resname);
	  ap->resid = pdb_atom[iatom].resseq;
	  ap->id = pdb_atom[iatom].id;
	} else {
	  *ap->occ = 0;
	}
      }
    }
  }
  
  for (k = 0; k < pe->nchains; k++) {
    confmp = Cbuf + pe->chainptr[k].index;
    {
      n = 0;
      for (iatom = 0; iatom < natoms_all; iatom++) {
	ap = confmp->head + iatom;
	n += *ap->occ;
      }
      confmp->natoms = n;
    }
  }	

  return 0;
}


void
print_conformer(Real *yp, int *occ, Real *tfac, Atom *templ, int natoms,
		int extra, char *remark, FILE *fp)
{
  int iatom, fwidth, fprec;
  Matrix3 Y;
  Atom *ap;
  char rec[BUFSIZ], buf[BUFSIZ];
  Real biso;

  if (extra) {
    fwidth = PDBM_EXTRA_x_LEN - 1;
    fprec = fwidth - 7;
  }

  fprintf(fp, "REMARK   1  %s\n", remark);
  Y = (Matrix3) yp;
  ap = templ;
  for (iatom = 0; iatom < natoms; iatom++, ap++) {
    if (occ == 0 || occ[iatom] != 0) {
      biso = tfac ? tfac[iatom] : 1.0;
      sprintf(rec, "ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
              ap->id, ap->name, ap->resname, ap->resid,
              Y[iatom][0], Y[iatom][1], Y[iatom][2], 1.0, biso);
      if (!extra)
        fprintf(fp, "%s\n", rec);
      else {
        sprintf(buf, " %*.*e %*.*e %*.*e", fwidth, fprec, Y[iatom][0],
                fwidth, fprec, Y[iatom][1], fwidth, fprec, Y[iatom][2]);
        fprintf(fp, "%-*s@%0*d%s\n",
                PDBM_EXTRA_FLAG_LOC, rec, PDBM_EXTRA_ID_LEN, ap->id, buf);
      }
    }
  }
  fprintf(fp, "TER\n");
}

/*
 *  returns observed data log likelihood presuming complete Yij
 *  workarea: natoms * 3 + 16
 */

Real
loglik(Ensemble Bbuf[], Real *Ybuf, Real *Mbuf, Real Zbuf[], Real diagS[],
       int natoms, int nconfms[], int ngroups, int nconfms_tot, Real term[])
{
  Real ll, s, *yp, *dp, d[3], *zp;
  int i, mi, iconfm, iens;
  Matrix3 M, Y;
  Ensemble *ensp;

  dp = term + 16;

  // 3nm log(2pi)
  ll = term[1] = 3 * natoms * nconfms_tot * (log(2.0*M_PI));
  // 3(m-r)log|Sigma|
  s = 0.0;
  for (i = 0; i < natoms; i++)
    s += log(diagS[i]);
  s *= 3 * (nconfms_tot - ngroups);
  ll += s;
  term[2] = s;

  yp = Ybuf;
  M = (Matrix3) Mbuf;
  zp = Zbuf;
  ensp = Bbuf;
  term[3] = term[4] = term[5] = term[6] = 0.0;
  for (iens = 0; iens < ngroups; iens++, ensp++) {
    mi = nconfms[iens];

    // 3log|Sigma + mi*V|
    term[3] += 3.0 * ensp->jvp->logdet;

    zerovectmatrix(dp, 3, natoms);
    // sum_j tr (Yij - M)^T Sigma^{-1} (Yij - M)
    s = 0;
    for (iconfm = 0; iconfm < mi; iconfm++) {
      Y = (Matrix3) yp;
      for (i = 0; i < natoms; i++) {
        d[0] = Y[i][0] - M[i][0];
        d[1] = Y[i][1] - M[i][1];
        d[2] = Y[i][2] - M[i][2];
        dp[i*3+0] += d[0];
        dp[i*3+1] += d[1];
        dp[i*3+2] += d[2];
        s += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / diagS[i];
      }
      yp += natoms * 3;
    }
    term[4] += s;

    // dp == Di = sum_j (Yij - M)
    // tr Zi^T Sigma^{-1} Di
    s = 0;
    for (i = 0; i < natoms; i++)
      s += (zp[i*3+0]*dp[i*3+0] + zp[i*3+1]*dp[i*3+1]
            + zp[i*3+2]*dp[i*3+2]) / diagS[i];
    term[5] += -s;

    zp += natoms * 3;
  }

  ll += term[3] + term[4] + term[5] + term[6];

  return ll;
}


static int
converged(Real *R, Real *t, Real *M, Real *Z, Real *V,
	  int ngroups, int nconfms_tot, int natoms,
	  Real loglik, Real del[], Real tol)
{
  static int firstcall = 1;
  static Real loglikold;
  int i;
  Real tolr, tolt, tolm, tolv, toll, *prevR, *prevt, *prevM, *prevZ, *prevV;
  Real rmax, tmax, mmax, zmax, vmax, ldiff;

  tolr = tol; tolt = tol; tolm = tol; tolv = tol; toll = tol;
  prevR = R + nconfms_tot * 9;
  prevt = t + nconfms_tot * 3;
  prevM = M + natoms * 3;
  prevZ = Z + natoms * 3 * ngroups;
  prevV = V + natoms * natoms;

  if (firstcall) {
    firstcall = 0;
  } else {
    /* R */
    rmax = 0.0;
    for (i = 0; i < nconfms_tot * 9; i++)
      if (fabs(R[i] - prevR[i]) > rmax) rmax = fabs(R[i] - prevR[i]);
    /* t */
    tmax = 0.0;
    for (i = 0; i < nconfms_tot * 3; i++)
      if (fabs(t[i] - prevt[i]) > tmax) tmax = fabs(t[i] - prevt[i]);
    /* M */
    mmax = 0.0;
    for (i = 0; i < natoms * 3; i++)
      if (fabs(M[i] - prevM[i]) > mmax) {
	mmax = fabs(M[i] - prevM[i]);
	if (mmax > 1.0e10) printf("%d %e %e\n", i, M[i], prevM[i]);
      }
    /* Z */
    zmax = 0.0;
    for (i = 0; i < natoms * 3 * ngroups; i++)
      if (fabs(Z[i] - prevZ[i]) > zmax) zmax = fabs(Z[i] - prevZ[i]);
    /* V */
    vmax = 0.0;
    for (i = 0; i < natoms*natoms; i++)
      if (fabs(V[i] - prevV[i]) > vmax) vmax = fabs(V[i] - prevV[i]);
    /* loglik */
    ldiff = (loglik - loglikold) / (fabs(loglik) + 1.0);
    /* */
    del[0] = ldiff;
    del[1] = rmax;
    del[2] = tmax;
    del[3] = mmax;
    del[4] = zmax;
    del[5] = vmax;
  }
  /* rotate */
  memcpy(prevR, R, sizeof(Real) * nconfms_tot * 9);
  memcpy(prevt, t, sizeof(Real) * nconfms_tot * 3);
  memcpy(prevM, M, sizeof(Real) * natoms * 3);
  memcpy(prevZ, Z, sizeof(Real) * natoms * 3 * ngroups);
  memcpy(prevV, V, sizeof(Real) * natoms*natoms);
  loglikold = loglik;
  return rmax < tolr && tmax < tolt && mmax < tolm && zmax < tolv
    && vmax < tolv && ldiff < toll;
}


int main(int argc, char *argv[])
{
  int ngroups, i, j, k, n, minid, maxid, minid0, maxid0, nchains, ient, nentries, natoms, nconfms_tot, iatom, itarget, iconfm, ichain, iter, niters, iens, ncmin, njvars;
  int nmi, YinM, fd, method;
  pdbEntry *pe;
  Conformer *cp;
  Atom *ap;
  FILE *fp, *fpcrd, *fpm;
  Matrix3 M, X;
  Real shift[3], *diagS, *workbuf, *zp, *rp, logL, sig0, sig0sq, *diagV;
  char buf[BUFSIZ], *bp;
  Real *fullV, *fullS, *rpp[2];

  progname = argv[0];

  ngroups = parse_args(argc, argv, &options, pdbentry);

  niters = options.niters;
  YinM = 0;
  method = (options.fix * METHOD_FIX)
    | (options.ran * METHOD_VAR)
    | (options.balanced * METHOD_EST);

  printf("ngroups = %d\n", ngroups);
  printf("niters lim = %d\n", niters);
  printf("tol: %.3e\n", options.tol);
  printf("alignM %d  isoVar %d  diagV %d  method %x\n",
	 options.alignM, options.isoVar, options.diagV, method);
  printf("impute %d", options.impute);
  if (options.impute)
    printf(" (missing data for atoms that present in at least %g %% of chains.)\n", options.impute_lowerlim * 100);
  else
    printf("\n");
  printf("Output files:");
  if (!options.covout && !options.mpdb && !options.ypdb)
    printf(" none\n");
  else {
    if (options.mpdb) printf("  M %s", options.mpdb); 
    if (options.ypdb) printf("  Y %s", options.ypdb); 
    if (options.covout) printf("  V&S %s", options.covout);
    printf("\n");
  }

  if (init_workarea(WORKBUFSIZ) != 0) {
    fprintf(stderr, "cannot allocate workbuf.\n"); return -1;
  }
  
  pdbchainpool = (pdbChain *) malloc(sizeof(pdbChain) * CHAINPOOL_SIZ);
  if (pdbchainpool == NULL) {
    fprintf(stderr, "cannot allocate chainpool.\n"); return -1;
  }

  nentries = ngroups;
  minid = minid0 = INT_MAX; maxid = maxid0 = 0;
  nchains = 0;			// TOTAL num of chains
  memset(idcnt, 0, sizeof(int) * (IDNUM_LIM+1));
  for (ient = 0; ient < nentries; ient++) {
    pe = &pdbentry[ient];
    fp = fopen(pe->file, "r");
    if (fp == NULL) {
      fprintf(stderr, "cannot open file %s.\n", pe->file); return -1;
    }

    prescan_pdbm(fp, pe, idcnt, IDNUM_LIM, pdbchainpool + nchains,
		 CHAINPOOL_SIZ-nchains, &minid, &maxid, id2ind);
    nchains += pe->nchains;
    if (minid < minid0) minid0 = minid;
    if (maxid > maxid0) maxid0 = maxid;
    fclose(fp);
  }
  minid = minid0; maxid = maxid0;

  /*
   * idcnt[i] is now the num of chains/frags (??not confms) for each landmark i.
   *
   * assign indeces to every landmark atoms
   * id2ind[] is a mapping from landmark id to index.
   * if idcnt[i] == 1, remove i from the landmark list.
   */
  // partial occ には対応しないが、missing にはいずれ対応する可能性もある。

  natoms = 0;
  id2ind[0] = -1;
  for (i = minid; i <= maxid; i++)
    if (idcnt[i] > 1) id2ind[i] = natoms++;
    else if (idcnt[i] == 1) {
      printf("Removing landmark #%d...\n", i);
      id2ind[i] = -1;
      idcnt[i] = 0;
    } else
      id2ind[i] = -1;
  /* natoms is the max num of atoms of the chain among the all chains */

  printf("\n");

#ifdef DEBUG
  printf("\nLandmark IDs in %d:%d\t chain segments by num obs:", minid, maxid);
  k = 0;
  idcnt[maxid+1] = 0;
  for (i = minid; i <= maxid; i++)
    if (idcnt[i] == idcnt[i+1]) k++;
    else {
      printf(" %d", idcnt[i]);
      printf("x%d", k+1);
      k = 0;
    }
  printf("\n");
#endif

  ncmin = nchains;
  if (options.impute)
    ncmin = ceil(nchains * options.impute_lowerlim);
#ifdef DEBUG
  printf("\nLandmarks observed in less than %d/%d chains are ignored.\n",
	 ncmin, nchains);
#endif

  k = 0;
  for (i = minid; i <= maxid; i++)
    if (idcnt[i] >= ncmin)
      id2ind[i] = k++;
    else {
      id2ind[i] = -idcnt[i];
      idcnt[i] = 0;
    }
  natoms = k;

#ifdef DEBUG
  printf("\nLandmark IDs in %d:%d\t chain segments by num obs:", minid, maxid);
  k = 0;
  idcnt[maxid+1] = 0;
  for (i = minid; i <= maxid; i++)
    if (idcnt[i] == idcnt[i+1]) k++;
    else {
      printf(" %d", idcnt[i]);
      printf("x%d", k+1);
      k = 0;
    }
  printf("\n");
  printf("natoms = %d\n\n", natoms);
#endif

  /*
   * pdbChain.index is used for addressing archetype, R, and t.
   */
  for (ichain = 0; ichain < nchains; ichain++)
    pdbchainpool[ichain].index = ichain;

  /*
   * 
   */
  nconfms_tot = nchains;

  /*
   * initialize (prepare memory areas)
   */
  i = mmfit_init(natoms, ngroups, nconfms_tot, pdbentry, method, &njvars);
  if (i < 0)
    abnormalexit("mmfit_init() failed", 0);

  nmi = nconfms[ngroups];

  /*
   * get coordinates
   */
  for (ient = 0; ient < nentries; ient++) {
    pe = &pdbentry[ient];
    fp = fopen(pe->file, "r");
    read_coordinates(fp, pe, id2ind, IDNUM_LIM, natoms, idcnt);
    // printf("nchains = %d\n", pe->nchains);
    fclose(fp);
  }

  /*
   * replace Xbuf with the contents of the binary coordinates files
   */
  if (binaryIO_crdext && *binaryIO_crdext) {
    rp = Xbuf;
    for (ient = 0; ient < nentries; ient++) {
      pe = &pdbentry[ient];
      strcpy(buf, pe->file);
      bp = buf + strlen(buf);
      while (bp-- > buf)
	if (*bp == '.' || *bp == '/') break;
      if (*bp != '.') {		// without extention
	bp = buf + strlen(buf);
	*bp++ = '.';
      } else			// with extention (or leading dot)
	bp++;
      strcpy(bp, binaryIO_crdext);
      fd = open(buf, O_RDONLY);
      if (fd < 0) {
	fprintf(stderr, "cannot open crdfile %s\n", buf);
	exit(-1);
      }
      i = sizeof(Real) * natoms * 3 * pe->nchains;
      if (read(fd, rp, i) < i) {
	fprintf(stderr, "cannot read nchains*natoms*3 values from %s\n", buf);
	exit(-1);
      }
      rp += natoms * 3 * pe->nchains;
      close(fd);
    }
  }  
  
  /*
   * collect missing patterns (oty's)
   * not yet implemented (removed) in the current version for cc cases.
   * 完全データのみに対応のため、このコードブロック occupancy_patterns は実装していない。
   */

  /*
   * prepare chain_template. It will be used for printing the mean structure.
   */
  cp = Cbuf;
  for (i = 0; i < nchains; i++) {
    /* copy atom/residue names/ids from the first conformer. */
    k = natoms;
    for (iatom = 0; iatom < natoms; iatom++) {
      ap = cp->head + iatom;
      if (chain_template[iatom].name[0] == '\0') {
        if (ap->name[0] != '\0') {
          strcpy(chain_template[iatom].name, ap->name);
          strcpy(chain_template[iatom].resname, ap->resname);
          chain_template[iatom].id = ap->id;
          chain_template[iatom].resid = ap->resid;
          k--;
        }
      } else
        k--;
    }
    if (k == 0) break;
  }
  if (k != 0)
    fprintf(stderr, "The names of %d atoms are not specified.\n", k);

#ifdef DEBUG
  printf("num observations for each atom\n");
  for (iatom = 0; iatom < natoms; iatom++) {
    printf("%04d: ", iatom);
    cp = Cbuf;
    for (iens = 0; iens < ngroups; iens++)  {
      k = 0;
      for (iconfm = 0; iconfm < nconfms[iens]; iconfm++) {
	k += cp->o[iatom];
	cp++;
      }
      printf("%3d", k);
    }
    printf("\n");
  }
#endif


  /*
   * prepare shared work area on the process of superposition
   */
  i = natoms*natoms*nmi + natoms*3*2;		// ??loglik() 要修正
  if (i < natoms*(3*natoms+6)) i = natoms*(3*natoms+6);
  workbuf = get_workarea("workbuf", i + 16, sizeof(Real));
  if (workbuf == 0) abnormalexit("workbuf alloc error", 0);

  /*
   * copy reftarget to M. use it as the reference target conformer
   * use the first conformer as the reftarget
   */
  itarget = 0;
  cp = &Cbuf[itarget];

  /*
   * missing atoms in the tentative M are imputed as being the centroid
   */
  M = (Matrix3) Mbuf;
  X = (Matrix3) cp->x;
  shift[0] = shift[1] = shift[2] = 0.0;
  n = 0;
  for (i = 0; i < natoms; i++)
    if (cp->o[i]) {
      shift[0] += X[i][0];
      shift[1] += X[i][1];
      shift[2] += X[i][2];
      n++;
    }
  shift[0] /= n; shift[1] /= n; shift[2] /= n;
  for (i = 0; i < natoms; i++)
    if (cp->o[i]) {
      M[i][0] = X[i][0];
      M[i][1] = X[i][1];
      M[i][2] = X[i][2];
    } else {
      M[i][0] = shift[0];
      M[i][1] = shift[1];
      M[i][2] = shift[2];
    }

  fullS = Sbuf;
  diagS = Sbuf + natoms*natoms * 2;
  /* S <- I */
  for (i = 0; i < natoms; i++)
    diagS[i] = 1.0;
  sig0sq = sig0 = 1.0;
  center_M(Mbuf, diagS, natoms, shift);	// with homoscedastic weight

  /* Z <- 0 */
  zerovectmatrix(Zbuf, natoms * ngroups, 3);

  /* Initial tij, Rij, yij */
  /* Note Rij Zi^T S^{-1} 1 == 0, since Zi == 0 */
  estimate_Rt_EM(Bbuf, Mbuf, diagS, ngroups, natoms, nconfms_tot, 0, YinM, workbuf);

  cp = Cbuf;
  zp = Zbuf;
  for (iens = 0; iens < ngroups; iens++) {
    for (iconfm = 0; iconfm < nconfms[iens]; iconfm++, cp++)
      calculate_Y(cp->x, cp->o, cp->r, cp->t, Mbuf, zp, natoms, cp->y);
    zp += natoms * 3;
  }

  /* V <- I, Vprev <- I */
  if (options.diagV)
    diagV = Vbuf + natoms*natoms * 2;
  else
    diagV = NULL;

  fullV = Vbuf;

  identityvectmatrix(fullV, natoms);
  identityvectmatrix(fullV + natoms * natoms, natoms);
  if (diagV)
    for (i = 0; i < natoms; i++)
      diagV[i] = 1.0;
  

  /* Initial guess for M */
  estimate_M_ML(Ybuf, Bbuf, ngroups, natoms, 1, diagV, Mbuf, workbuf);
  // estimate_M_ML(Ybuf, Bbuf, ngroups, natoms, njvars, diagV, Mbuf, workbuf);

  /* なぜ? */
  if ((method & METHOD_EST) == METHOD_EST_BML)
    center_MZ(Mbuf, Zbuf, diagS, nconfms, ngroups, nconfms_tot, natoms, shift, workbuf);
  else
    center_M(Mbuf, diagS, natoms, shift);

  printf("\nAll conformers have been superimposed on the %d-th conformer.\n\n",
	 itarget);
  
  estimate_Uaux(Jbuf, fullV, diagV, diagS, natoms, workbuf);
  estimate_Z(Bbuf, diagV, Mbuf, natoms, ngroups, workbuf);

  workbuf[0] = 0.0;
  logL = loglik(Bbuf, Ybuf, Mbuf, Zbuf, diagS,
		natoms, nconfms, ngroups, nconfms_tot, workbuf);

  rp = workbuf;
  printf("#loglik = %18.12e (%.3e %.3e %.3e %.3e %.3e %.3e %.3e)\n",
	 logL, rp[0], rp[1], rp[2], rp[3], rp[4], rp[5], rp[6]);

  (void)
    converged(Rbuf, tbuf, Mbuf, Zbuf, fullV, ngroups, nconfms_tot, natoms,
	      logL, workbuf, options.tol);

  printf("\niter\t    loglik  & 'del': changes in the estimates\n");


  for (iter = 0; iter < niters; iter++) {

    /*
     * M
     */
    if ((method & METHOD_FIX) == METHOD_FIX_Q) {
      estimate_M_EM(Ybuf, Zbuf, nconfms, ngroups, natoms, Mbuf);
      center_MZ(Mbuf, Zbuf, diagS, nconfms, ngroups, nconfms_tot, natoms,
		shift, workbuf);
    } else {
      estimate_M_ML(Ybuf, Bbuf, ngroups, natoms, njvars, diagV, Mbuf, workbuf);
      center_M(Mbuf, diagS, natoms, shift);
    }

    /*
     * Rij, tij, and Yij
     */
    if ((method & METHOD_FIX) == METHOD_FIX_Q)
      estimate_Rt_EM(Bbuf, Mbuf, diagS, ngroups, natoms, nconfms_tot, 0, YinM, workbuf);
    else
      estimate_Rt_ML(Bbuf, Mbuf, diagV, diagS, ngroups, natoms, 0, workbuf);

    cp = Cbuf;
    zp = Zbuf;
    for (iens = 0; iens < ngroups; iens++) {
      for (iconfm = 0; iconfm < nconfms[iens]; iconfm++) {
	calculate_Y(cp->x, cp->o, cp->r, cp->t, Mbuf, zp, natoms, cp->y);
	cp++;
      }
      zp += natoms * 3;
    }

    /* V and S */
    if ((method & METHOD_EST) != METHOD_EST_BML) {
      estimate_V(Zbuf, Jbuf, nconfms, ngroups, natoms, method, diagV, fullV, /* Vani */ 0, workbuf);
      estimate_S(Ybuf, Mbuf, Zbuf, Jbuf, nconfms_tot, nconfms, ngroups, natoms,
		 method, diagV, diagS, fullS, /* Sani */ 0, workbuf);
    } else {		// if METHOD_EST == BML 未対応
      sig0sq = estimate_SV_BML(Ybuf, Mbuf, natoms, nconfms_tot, nconfms,
			       ngroups, diagS, fullV, workbuf);
    }
    /*
    for (i = 0; i < natoms*natoms; i++)
      fullS[i] = workbuf[i];
    */
    if (options.isoVar) {
      sig0sq = 0;
      for (i = 0; i < natoms; i++)
	sig0sq += diagS[i];
      sig0sq /= natoms;
      for (i = 0; i < natoms; i++)
	diagS[i] = sig0sq;
    }
    
    estimate_Uaux(Jbuf, fullV, diagV, diagS, natoms, workbuf);
    /*
     * MAP estimates for Zi using tentative Ui
     */
    estimate_Z(Bbuf, diagV, Mbuf, natoms, ngroups, workbuf);
    
    rp = workbuf;
    logL = loglik(Bbuf, Ybuf, Mbuf, Zbuf, diagS,
		  natoms, nconfms, ngroups, nconfms_tot, rp);
    rp = workbuf + 16;
    i = converged(Rbuf, tbuf, Mbuf, Zbuf, fullV, ngroups, nconfms_tot, natoms,
		  logL, rp, options.tol);
    printf("#%5d  -2logLo %18.12e (rdel %.3e)\n",
	   iter+1, logL, rp[0]);
    printf(" |rdel| R %.3e  t %.3e  M %.3e  Z %.3e  V %.3e\n\n",
	   rp[1], rp[2], rp[3], rp[4], rp[5]);
    
    if (i != 0) break;
  
  }

  estimate_V(Zbuf, Jbuf, nconfms, ngroups, natoms, method, diagV, fullV, Vbuf+natoms*natoms*2+natoms*2, workbuf);
      estimate_S(Ybuf, Mbuf, Zbuf, Jbuf, nconfms_tot, nconfms, ngroups, natoms,
		 method, diagV, diagS, fullS, Sbuf+natoms*natoms*2+natoms*2, workbuf);


  rp = workbuf;
  printf("loglik = %18.12e\n\t(%.4e %.4e %.4e %.4e %.4e %.4e)\n",
	 logL, rp[0], rp[1], rp[2], rp[3], rp[4], rp[5]);
  //printf("sig0 = %.8e\n", sig0);


  fd = -1;
  if (binaryIO_binout && *binaryIO_binout) {
    fd = open(dbinout, O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
      fprintf(stderr, "cannot open %s for writing\n", dbinout);
      exit(-1);
    }
  }

  /*
   * the coordinates
   */
  rpp[0] = diagV;
  diagV = Vbuf + natoms*natoms * 2 + natoms;
  for (i = 0; i < natoms; i++)
    diagV[i] = sqrt(fullV[i*natoms+i]);
  rpp[1] = diagS;
  diagS = Sbuf + natoms*natoms * 2 + natoms;
  for (i = 0; i < natoms; i++)
    diagS[i] = sqrt(diagS[i]);

  /* M */
  if (options.mpdb) {
    fpm = fopen(options.mpdb, "w");
    strcpy(buf, "=== mean. V^{1/2} on the temp factor column.");
    fprintf(fpm, "MODEL\n");
    print_conformer(Mbuf, 0, diagV, chain_template, natoms, 0, buf, fpm);
    fprintf(fpm, "ENDMDL\n");
    strcpy(buf, "====== mean. S^{1/2} on the temp factor column.");
    fprintf(fpm, "MODEL\n");
    print_conformer(Mbuf, 0, diagS, chain_template, natoms, 0, buf, fpm);
    fprintf(fpm, "ENDMDL\n");
    fclose(fpm);
  }
  //
  if (fd >= 0)
    write(fd, Mbuf, sizeof(Real) * natoms * 3);
  // restore
  diagV = rpp[0];
  diagS = rpp[1];

  /* Yi */
  if (options.ypdb) {
    fpcrd = fopen(options.ypdb, "w");
    for (ient = 0; ient < nentries; ient++) {
      pe = &pdbentry[ient];
      fprintf(fpcrd, "REMARK   1  ====== %s\n", pe->file);
      read_trans_write_pdbm(fpcrd, pe, Rbuf, tbuf, 1);
    }
    fclose(fpcrd);
  }
  //
  if (fd >= 0)
    write(fd, Ybuf, sizeof(Real) * natoms * 3 * nconfms_tot);

  /* Zi */
  if (fd >= 0)
    write(fd, Zbuf, sizeof(Real) * natoms * 3 * ngroups);


  if (options.covout || fd >= 0) {
    if ((method & METHOD_EST) != METHOD_EST_BML) {
      /* EM updates (estimates) for V and S */
      /* write the current estimates */
      if (options.covout) {
	fp = fopen(options.covout, "a");
	rp = Vbuf;
	for (i = 0; i < natoms; i++) {
	  for (j = 0; j < natoms; j++)
	    fprintf(fp, "%.16e ", rp[i*natoms+j]);
	  fprintf(fp, "\n");
	}
      }
      if (fd >= 0)
	write(fd, Vbuf, sizeof(Real) * natoms * natoms);

      if (options.covout) {
	rp = fullS;
	for (i = 0; i < natoms; i++) {
	  for (j = 0; j < natoms; j++)
	    fprintf(fp, "%.16e ", rp[i*natoms+j]);
	  fprintf(fp, "\n");
	}
	fclose(fp);
      }
      if (fd >= 0)
	write(fd, fullS, sizeof(Real) * natoms * natoms);

    }
  }

  if (fd >= 0) {
    write(fd, Vbuf+natoms*natoms*2+natoms*2, sizeof(Real)*natoms*natoms*9);
    write(fd, Sbuf+natoms*natoms*2+natoms*2, sizeof(Real)*natoms*natoms*9);
  }
  if (fd >= 0)
    close(fd);


  if (binaryIO_binout && *binaryIO_binout) {
    fd = open(ibinout, O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
      fprintf(stderr, "cannot open %s for writing\n", ibinout);
      exit(-1);
    }
    write(fd, &natoms, sizeof(int));
    write(fd, &ngroups, sizeof(int));
    write(fd, nconfms, sizeof(int) * ngroups);
    close(fd);
  }

  return 0;
}
