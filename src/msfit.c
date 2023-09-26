#include <limits.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
       #include <sys/types.h>
       #include <sys/stat.h>
       #include <fcntl.h>
       #include <unistd.h>
#include "molec.h"
#include "pdb.h"
#include "matrix.h"
#include "memory.h"
#include "msfit.h"


#ifdef MKLROOT
int UseMKL = 1;
#else
int UseMKL = 0;
#endif
   
/*
 * msfit.c for complete-case using the diagonal & isotropic working covariance.
 * In complete-case, there is no missing/incomplete atoms in any Xj.
 */

char *progname;

#define PDB_NENT_MAX 10000
#define IDNUM_LIM 10000			/* upper limit of landmark ID nums  */
#define CHAINPOOL_SIZ (1<<18)
#define OCC_EPS (1.0e-6)		/* to cope with rounding errors     */

int id2ind[IDNUM_LIM+2];		/* 0..IDNUM_LIM and a sentinel      */
int idcnt[IDNUM_LIM+2];
int idcnttmp[IDNUM_LIM+2];
pdbEntry pdbentry[PDB_NENT_MAX];
pdbChain *pdbchainpool = 0;		/* storage for pdbChains */

#define NSITES_MAX 10000		/* max # of sites(atoms) in an entry */
pdbAtom pdb_atom[NSITES_MAX];		/* reuse pdb_atom and pdb_site       */
pdbSite pdb_site[NSITES_MAX];		/*   for all entries                 */

/* pointers to storages that will be allocated in msfit_init() */
Atom *Abuf = 0;
Conformer *Cbuf = 0;
Real *Mbuf = 0, *Xbuf = 0, *Ybuf = 0;
Real *tbuf = 0, *Rbuf = 0, *gbuf = 0;
int *Obuf = 0;
struct covar *covar = 0;
Real *covarR = 0;
int *covarI = 0;
Atom *chain_template;
Real *Sbuf = 0;
Real *diagS = 0, *isoS = 0, *anisoS = 0;

#define NLOCSMISS_LIM 62		/* must be < sizeof(int64_t)-1 */
#define NUMOTYS_LIM 5000
int locmiss[NLOCSMISS_LIM], *missarrayptr[NUMOTYS_LIM], nobsonpatn[NUMOTYS_LIM];
int64_t misspatn[NUMOTYS_LIM];

int nlocsmiss;


static struct arguments options;
static int verbose;

#define WORKBUFSIZ (1<<30)



int
convergence_attained(Real *R, Real *t, int nconfms, Real *M, int natoms,
                     Real del[], Real tol[]);
int
partition_into_conformers(Real occ[], int noccs[], int nseqs, Real eps,
                          Real gamma[]);
int
prescan_pdbm(FILE *fp, pdbEntry *pe, int idcnt[], int id_lim,
             pdbChain *chainbuf, int nchains_lim, int *minidp, int *maxidp,
             int *work);
int
msfit_init(int natoms, int nchains, int nconforms, int ngammas[], Real *gp);
int
read_coordinates(FILE *fp, pdbEntry *pe, int id2ind[], int id_lim,
                 Real occ_eps, int natoms_all, Conformer *confmpool, int work[]);
int
allocate_covar(int natoms);
void
print_conformer(Real *yp, int *occ, Atom *templ, int natoms, int extra,
                char *remark, FILE *fp);





void msfit_finalize()
{
  // if (progname) free(progname);
  if (gbuf) free(gbuf);
  if (Abuf) free(Abuf);
  if (Cbuf) free(Cbuf);
  if (Obuf) free(Obuf);
  if (tbuf) free(tbuf);
  if (Rbuf) free(Rbuf);
  if (Mbuf) free(Mbuf);
  if (Xbuf) free(Xbuf);
  if (covar) free(covar);
}  

void abnormalexit(char *message, int code)
{
  if (message && code)
    fprintf(stderr, "(%d) %s\n", code, message);
  else if (message)
    fprintf(stderr, "%s\n", message);
  if (progname) free(progname);
  msfit_finalize();
  exit(code);
}


char usagemessage[] =
  "  -h : homogeneous weight == OLS\n"
  "  -B : ten Berge modification (1977), i.e., fit each on ave of the others\n" 
  "  -a : centerize & reorientate M to place it along PAs, before output\n"
  "  -o v wcovfile : dump the final working cov matrix V to the file wcovfile\n"
  "  -o s|S covfile : output sample cov estimates S to covfile\n"
  "  -o c crdfile : output Yj's & M in high precision to the file crdfile\n"
  "  -c : reduce to the complete case.\n"
  "  -V : verbose\n"
  "  -r : molecule-based weighting using rmsd\n"
  "  -v d : diagonal variance-covariance matrix\n"
  "  -e eR et eM eL niters: stopping criteria\n"
  "      R: epsM^(2/3)*10^eR,  et: epsM^(2/3)*10^et,  eM: epsM^(2/3)*10^eM\n"
  "      loglik: epsM^(1/2)*10^eL,  niters: upper lim of num of iterations\n"
  "  -t cnf: initial reference conformer [default: the first conformer]\n"
  "      examples of cnf: 1xray.pdbm:B  2nmr.pdbm:2 (fname cannot be omitted)\n"
  ;


/*
 * binary I/O:
 * - Requirements:
 *    occ == 1. e.g., nconfoms == 1 for all archetypes.
 * -b i crdext
 *   overwrite Xbuf with the contents of pefile.crdext files,
 *   where pefile is a basename of the pdbm files.
 * -b o binout
 *   output S, M, Yi, Saniso, Ri to the file 'd'binout
 *   and natoms, nchains to the file 'i'binout.
 */
char *binaryIO_crdext = 0;
char *binaryIO_binout = 0;
char ibinout[BUFSIZ], dbinout[BUFSIZ];


void usage(char *progname)
{
  printf("usage: %s [OPTIONS] pdbfile_1 ... pdbfile_n\n", progname);
  printf("%s", usagemessage);
  exit(-1);
}

static void
init_options(struct arguments *op)
{
  op->ols = 0;
  op->ave_of_others = 0;
  op->verbose = 0;
  op->alpha = 0;
  op->tau = 0;
  op->covty = COV_DIAGONAL;
  op->lgkmax = 4;
  op->place_alongPA = 0;
  op->covfile[0] = 0;
  op->wcovfile[0] = 0;
  op->reftarget[0] = 0;
  op->r1weight = 0;
  op->impute0 = 0;
  op->impute0samplecovar = 0;
  op->reduceCC = 0;
  op->initweight = 0;
}

int
parse_arguments(int argc, char *argv[],
		pdbEntry pdbentry[], struct arguments *opts)
{
  int covtyinit, c, i, ii, j, nentries;
  pdbEntry *pe;

  init_options(opts);
  covtyinit = opts->covty;
  
  for (i = 1; i < argc; i++) {
    if (*argv[i] != '-')
      break;
    ii = i;
    for (j = 1; argv[i][j]; j++) {
      c = argv[i][j];
      switch (c) {
      case 'V':
	opts->verbose = 1; break;
      case 'h':
	opts->ols = 1; break;
      case 'B':
	opts->ave_of_others = 1; break;
      case 'a':
	opts->place_alongPA = 1; break;
	/*
      case 'r':
	opts->r1weight = 1;
	break;
      case 'c':
	opts->reduceCC = 1;
	break;
	*/
      case 'v':
	opts->covty &= ~covtyinit;
	switch (*argv[++i]) {
	case 'd':
	  opts->covty |= COV_DIAGONAL;
	  break;
	  /*
	case 's':
	  opts->covty |= COV_SHRINK;
	  opts->alpha = atof(argv[++i]);
	  break;
	case 't':
	  opts->covty |= COV_TRUNCATE;
	  if (argv[i][1] == 's')
	    opts->covty |= COV_SD_MEST;
	  else if (argv[i][1] == 'r')
	    opts->covty |= COV_COR_KENDALL;
	  else if (argv[i][1] == 'b')
	    opts->covty |= COV_SD_MEST | COV_COR_KENDALL;
	  opts->tau = atof(argv[++i]);
	  if (isdigit(*argv[i+1]))
	    opts->lgkmax = atof(argv[++i]);
	  break;
	  */
	default:
	  usage(argv[0]);
	}
	break;
      case 'o':
	switch (*argv[++i]) {
	case 'v':
	  i++;
	  strcpy(opts->wcovfile, argv[i]);
	  break;
	case 'c':
	  i++;
	  strcpy(opts->crdfile, argv[i]);
	  break;
	case 'S':
	  i++;
	  strcpy(opts->covfile, argv[i]);
	  break;
	case 's':
	  i++;
	  strcpy(opts->covfile, argv[i]);
	  break;
	default:
	  usage(argv[0]);
	}
	break;
      case 'e':
	opts->tol[0] *= pow(10.0, atof(argv[++i]));
	opts->tol[1] *= pow(10.0, atof(argv[++i]));
	opts->tol[2] *= pow(10.0, atof(argv[++i]));
	opts->tol[3] *= pow(10.0, atof(argv[++i]));
	opts->tol[4] =  atoi(argv[++i]);
	break;
      case 't':
	++i;
	strcpy(opts->reftarget, argv[i]);
	break;
      case 'W':
	opts->initweight = atoi(argv[++i]);
	break;
      case 'b':
	switch (*argv[++i]) {
	case 'i':
	  binaryIO_crdext = strdup(argv[++i]);
	  break;
	case 'o':
	  binaryIO_binout = strdup(argv[++i]);
	  sprintf(ibinout, "i%s", binaryIO_binout);
	  sprintf(dbinout, "d%s", binaryIO_binout);
	  break;
	}
	break;
      default:
	usage(argv[0]);
      }
      if (i > ii)
	break;
    }
  }

  nentries = 0;
  for ( ; i < argc; i++) {
    pe = &pdbentry[nentries];
    pe->id = nentries;
    nentries++;
    pe->file = argv[i];
  }
  if (nentries == 0)
    usage(argv[0]);
  return nentries;
}



/*
 * return the (serial) index of reftarget conformer specied as str
 * given as filename[:chain/model id]
 * Notice: 
 * - if str is empty, this func returns 0 (index of the first conformer) 
 * - chains and conformerss are virtually the same and are given same indeces.
 * - returns the smallest index if the model contains multiple conformers.
 * - returns the smallest index if chain or model id is omitted.
 */
int
indexof_reftarget_conformer(char *str, pdbEntry *pe, int nentries,
			    Conformer *confmpool)
{
  char *bp, *cp;
  int refconfm, k, cid, mdlid, ient;
  Conformer *confmp;

  if (str[0] == 0)
    return 0;		/* the first conformer */

  bp = strtok_r(str, ":", &cp);		/* filename */
  for (ient = 0; ient < nentries; ient++, pe++)
    if (strcmp(bp, pe->file) == 0)
      break;
  if (ient >= nentries)
    abnormalexit("reftarget file not found\n", 0);
  bp = strtok_r(NULL, ":", &cp);	/* chain/model id */
  if (*bp) {
    if (isdigit(*bp)) {			/* model num. chain id unknown */
      mdlid = atoi(bp);
      for (k = 0; k < pe->nchains; k++)
	//if (match_two_ids(pe->chainptr[k].id, pe->id, -1, mdlid) == 0)
	if (pe->chainptr[k].model == mdlid)
	  break;
    } else {
      cid = *bp;
      for (k = 0; k < pe->nchains; k++)
	if (cid == pe->chainptr[k].chainid)
	  break;
    }
    if (k >= pe->nchains)
      abnormalexit("reftarget chain not found\n", 0);
  } else
    k = 0;
  confmp = confmpool + pe->chainptr[k].index;
  refconfm = confmp->index;
  return refconfm;
}


/*
 * work: >= natoms * 6
 */
Real
loglik(Real *Mbuf, Real *Ybuf, Conformer *cp, int nconfms, int natoms,
       Real *diagS, Real *work)
{
  Real *mp, *yp, *dp, *xp, s, s1, logL, chisq, sumsq;
  Matrix3 M, Y, D;
  int iconfm, i, k, iatom, nu, sumnu;
  Real t1, t2, t3;

  xp = work;
  dp = xp + natoms * 3;
  mp = Mbuf;
  M = (Matrix3) mp;
  D = (Matrix3) dp;
  yp = Ybuf; 

  chisq = sumnu = 0;
  sumsq = 0;		/* sum of the unweighted squared residuals */
  t1 = t2 = t3 = 0;
  logL = 0.0;
  for (iconfm = 0; iconfm < nconfms; iconfm++, cp++, yp += natoms*3) {
    Y = (Matrix3) yp;
    /* D <- Yj - M */
    i = 0;
    for (iatom = 0; iatom < natoms; iatom++)
      if (cp->o[iatom]) {
	D[i][0] = Y[iatom][0] - M[iatom][0];
	D[i][1] = Y[iatom][1] - M[iatom][1];
	D[i][2] = Y[iatom][2] - M[iatom][2];
	i++;
      }
    if (i != natoms) 
      fprintf(stderr, "ERROR\n"), exit(-1);

    /* X <- Sj^{-1} D = Sj^{-1} (Yj - M) */
    nu = natoms;
    sumnu += nu;
    /*
     * X is stored on xp[] in column-major order.
     * This is a trace of the code for full matrix covars,
     * in which dpotrs_() or cholesky_solve() were used.
     */
    i = 0;
    for (iatom = 0; iatom < natoms; iatom++)
      if (cp->o[iatom]) {
	xp[i+0*nu] = D[i][0] / diagS[iatom];
	xp[i+1*nu] = D[i][1] / diagS[iatom];
	xp[i+2*nu] = D[i][2] / diagS[iatom];
	i++;
      }

    /* t1 += tr D^T X = tr (Yj-M)^T Sj^{-1} (Yj-M) */
    s = s1 = 0.0;
    for (i = 0; i < 3; i++)
      for (k = 0; k < nu; k++) {
	s += D[k][i] * xp[k+i*nu];
	s1 += D[k][i] * D[k][i];
      }
    cp->rmsd = sqrt(s / nu);
    cp->rmsd_uw = sqrt(s1 / nu);
    t1 += s;
    chisq += s;
    sumsq += s1;

    /* 3 log det Sj */
    for (iatom = 0; iatom < natoms; iatom++)
      if (cp->o[iatom])
	t2 += 3.0 * log(diagS[iatom]);
    /* 3 nuj log 2pi */
    t3 += 3.0 * nu * log(2.0 * M_PI);
  }
  logL = t1 + t2 + t3;
  work[0] = logL;
  work[1] = t1;
  work[2] = t2;
  work[3] = t3;
  work[4] = sqrt(sumsq / sumnu);
  chisq /= ((sumnu - natoms) * 3.0);
  chisq /= natoms / (natoms - 2);
  work[5] = chisq;
  return logL;
}


int
main(int argc, char *argv[])
{
  int nentries, nchains, nconforms, natoms, niters;
  int minid, maxid, ient, ichain, iconfm, iatom, iter, i, j, k, ic;
  int refconfm, found, n, nu, bulkUpdate, fd;
  int *ngammas, itarget;
  Real *gamma, *gp, *rp, *dp, *spwork, *mat0n3, matI3[9], rotbuf[9], *yp;
  Real logL, chisq, uwrmsd, *tol, machineps, shift[3], s;
  Real qual[10], sx[3], *mp;
  Matrix3 M, X, V, R;
  pdbEntry *pe;
  pdbChain *chainp;
  Conformer *cp;
  Atom *ap;
  FILE *fp;
  char reftarget[BUFSIZ], buf[BUFSIZ], *bp;


  progname = strdup(argv[0]);
  machineps = DBL_EPSILON;
  tol = options.tol;
  tol[0] = tol[1] = tol[2] = pow(machineps, 2./3.0);
  /* are tols for R, t & M, and the tol for loglik is */
  tol[3] = pow(machineps, 0.5);
  /* limit for niters is */
  tol[4] = 200;

  nentries = parse_arguments(argc, argv, pdbentry, &options);
#ifdef DEBUG
  verbose = 1;
#endif

  niters = tol[4];

  if (init_workarea(WORKBUFSIZ) != 0) {
    fprintf(stderr, "cannot allocate workbuf.\n"); return -1;
  }

  pdbchainpool = (pdbChain *) malloc(sizeof(pdbChain) * CHAINPOOL_SIZ);
  if (pdbchainpool == NULL) {
    fprintf(stderr, "cannot allocate chainpool.\n"); return -1;
  }

  minid = INT_MAX; maxid = 0;
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
    fclose(fp);
  }
  /*
   * idcnt[i] is now the num of chains/frags (not confms) for each landmark i.
   *
   * assign indeces to every landmark atoms
   * id2ind[] is a mapping from landmark id to index.
   * if idcnt[i] == 1, remove i from the landmark list.
   */

  if (nchains == 1) {
    /* fake two confms to pass the following landmark counting */
    for (i = minid; i <= maxid; i++)
      if (idcnt[i] == 1) idcnt[i] += 1;
  }

  natoms = 0;
  id2ind[0] = -1;
  for (i = minid; i <= maxid; i++)
    if (idcnt[i] > 1) id2ind[i] = natoms++;
    else if (idcnt[i] == 1) {
      printf("Removing landmark #%d...\n", i);
      id2ind[i] = -1;
      idcnt[i] = 0;
      chainp = pdbchainpool;
      k = 0;
      for (ichain = 0; ichain < nchains; ichain++, chainp++) {
	rp = chainp->occ_frac;
	for (j = 0; j < chainp->occ_nparts; j++) {
	  if (chainp->occ_firstid[j] == i) {
	    if (k != 0)
	      fprintf(stderr, "ERROR: multiple chains for %d\n", i), exit(-1);
	    n = 0;
	    for (k = j; k < chainp->occ_nparts-1; k++) {
	      chainp->occ_firstid[k] = chainp->occ_firstid[k+1];
	      chainp->occ_nfracs[k] = chainp->occ_nfracs[k+1];
	      n += chainp->occ_nfracs[k];
	    }
	    memmove(rp, rp + chainp->occ_nfracs[j], sizeof(Real) * n);
	    chainp->occ_nparts--;
	    k = -1;
	    // break;
	  }
	  rp += chainp->occ_nfracs[j];
	}
      }
      // if (k >= 0)
      // fprintf(stderr, "ERROR: id %d not found\n", i), exit(-1);
    } else
      id2ind[i] = -1;
  /* natoms is the max num of atoms of the chain among the all chains */

  if (nchains == 1) {
    for (i = minid; i <= maxid; i++)
      if (idcnt[i] == 2) idcnt[i] = 1;
  }

  printf("\n\n");
  if (verbose) {
  chainp = pdbchainpool;
  k = 0;
  for (ichain = 0; ichain < nchains; ichain++, chainp++) {
    rp = chainp->occ_frac;
    for (j = 0; j < chainp->occ_nparts; j++) {
      printf("%d ", chainp->occ_firstid[j]);
      for (i = 0; i < chainp->occ_nfracs[j]; i++)
	printf("%f ", *rp++);
      printf(";\t");
    }
    printf("\n");
  }
  }

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

  if (options.reduceCC) {	/* reduce to the complete case */
    k = 0;
    for (i = minid; i <= maxid; i++)
      if (idcnt[i] == nchains)
	id2ind[i] = k++;
      else
	id2ind[i] = -idcnt[i];
    natoms = k;
  }

  /*
   * pdbChain.index is used for addressing archetype, R, and t.
   */
  chainp = pdbchainpool;
  for (ichain = 0; ichain < nchains; ichain++, chainp++)
    chainp->index = ichain;
  
  /*
   * determine partition coefficients of conformers of each archetype
   */
  nconforms = 0;		// TOTAL num of conformers
  gamma = get_workarea("gamma", CHAINPOOL_SIZ*10, sizeof(Real));
  ngammas = get_workarea("ngammas", CHAINPOOL_SIZ, sizeof(int));
  dp = get_workarea("prefixsum", 1000, sizeof(Real));
  gp = gamma;
  chainp = pdbchainpool;
  for (ichain = 0; ichain < nchains; ichain++) {
    
    n = 0;
    rp = chainp->occ_frac;
    for (j = 0; j < chainp->occ_nparts; j++) {
      s = 0.0;
      for (i = 0; i < chainp->occ_nfracs[j]; i++) {
	s += rp[n];
	dp[n] = s;
	n++;
      }
    }
    for (i = 0; i < n-1; i++) {
      k = i;
      for (j = i+1; j < n; j++)
	if (dp[j] < dp[k])
	  k = j;
      if (k > i) {
	s = dp[i];
	dp[i] = dp[k];
	dp[k] = s;
      }
    }

    gp[0] = dp[0];
    i = j = 1;
    for (; i < n; i++)
      if (dp[i] != dp[i-1])
	gp[j++] = dp[i] - dp[i-1];
    ngammas[ichain] = j;

    /*
      ngammas[ichain]
	= partition_into_conformers(chainp->occ_frac, chainp->occ_nfracs,
				    chainp->occ_nparts, OCC_EPS, gp);
    */

      if (ngammas[ichain] <= 0)
	abnormalexit("partition_into_conformers() failed", 0);

      nconforms += ngammas[ichain];
      gp += ngammas[ichain];
      chainp++;

      if (exceed_ulim("gamma", gp) 
	  || exceed_ulim("ngammas", ngammas))
	abnormalexit("insufficient workarea (gamma,ngammas).", 0);
  }			/* for (ichain = 0; ichain < pe->nchains; ichain++) */

  /*
  if (nconforms == 2 && nentries != 2)
    abnormalexit("Place the two conformers into seperate files.\n", 0);
  */

  if (verbose) {
    printf("\n\nChains and the partition coefficients.\n");
    gp = gamma;
    pe = pdbentry;
    nu = 0;	// accum num of chains
    for (ient = 0; ient < nentries; ient++, pe++) {
      printf("\n%s\n", pe->file); 
      for (ichain = 0; ichain < pe->nchains; ichain++) {
	chainp = pe->chainptr + ichain;
	printf("  ichain = %d\n|    Occupancy :", ichain);
	k = 0;
	for (i = 0; i < chainp->occ_nparts; i++) {
	  printf(i==0? "%5d %2d":"|\t\t%5d %2d", i, chainp->occ_nfracs[i]);
	  for (j = 0; j < chainp->occ_nfracs[i]; j++)
	    printf(" %f", chainp->occ_frac[k++]);
	  printf("\n");
	}
	printf("     gamma's =");
	for (i = 0; i < ngammas[nu]; i++)
	  printf(" %g", gp[i]);
	printf("\n");
	gp += ngammas[nu];
	nu++;
      }
    }
  }


  /*
   * stop in case of incomplete data
   */
  ic = 0;
  gp = gamma;
  pe = pdbentry;
  nu = 0;	// accum num of chains
  for (ient = 0; ient < nentries; ient++, pe++) {
    for (ichain = 0; ichain < pe->nchains; ichain++) {
	chainp = pe->chainptr + ichain;
	k = 0;
	for (i = 0; i < chainp->occ_nparts; i++) {
	  for (j = 0; j < chainp->occ_nfracs[i]; j++)
	    if (chainp->occ_frac[k++] < 1.0) 
	      ic = 1;
	}
	for (i = 0; i < ngammas[nu]; i++)
	  if (gp[i] < 1.0)
	    ic = 1;
	gp += ngammas[nu];
	nu++;
    }
  }
  if (ic)
    abnormalexit("Incomplete data.", 0);


  /*
   * initialize (prepare memory areas)
   */
  i = msfit_init(natoms, nchains, nconforms, ngammas, gamma);
  if (i < 0)
    abnormalexit("msfit_init() failed", 0);

  remove_workarea("prefixsum");
  remove_workarea("ngammas");
  remove_workarea("gamma");


  /*
   * get coordinates
   */
  for (ient = 0; ient < nentries; ient++) {
    pe = &pdbentry[ient];
    fp = fopen(pe->file, "r");
    read_coordinates(fp, pe, id2ind, IDNUM_LIM, OCC_EPS, natoms, Cbuf, idcnt);
    fclose(fp);
  }

  /*
   * replace Xbuf with the contents of the binary coordinates files.
   * occ must be 1 for all marked atoms.
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

  nlocsmiss = 
    occupancy_patterns(Cbuf, nconforms, natoms, NLOCSMISS_LIM, NUMOTYS_LIM,
		       locmiss, misspatn, missarrayptr, nobsonpatn, &numotys);
  if (nlocsmiss < 0)
    abnormalexit("error occupancy_patterns()", nlocsmiss);
  */

  if (verbose) {
    printf("\n\nConformers\n");
    cp = Cbuf;
    for (i = 0; i < nchains; i++, cp++) {
      printf("Conformer %03d gamma = %g  natoms = %d/%d\n",
	     cp->index, *(cp->g), cp->natoms, natoms);
    }
    printf("\n\n");
  }


  /*
   * allocate memory for covariance matrices
   */
  allocate_covar(natoms);

  /*
   * prepare chain_template. It will be used for printing the mean structure.
   */
  cp = Cbuf;
  /* copy atom/residue names/ids from the first conformer. */
  for (i = 0; i < nchains; i++, cp++) {
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


  refconfm = 0;

  mat0n3 = get_workarea("mat0n3", natoms * 3, sizeof(Real));
  if (mat0n3 == 0) abnormalexit("mat0n3 alloc error", 0);
  zerovectmatrix(mat0n3, natoms, 3);
  identityvectmatrix(matI3, 3);

  /*
   * prepare shared work area on the course of superposition
   * loglik: natoms * 6, average_Y: natoms,
   * estimate_S: natoms^2 + natoms * 3
   */
  i = natoms * 10 + 20;
  if (i < natoms * natoms + natoms * 3)
    i = natoms * natoms + natoms * 3;
  spwork = get_workarea("superposwork", i, sizeof(Real));
  if (spwork == 0) abnormalexit("superposwork alloc error", 0);

  /*
   * copy reftarget to M. use it as the reference target conformer
   */
  refconfm = indexof_reftarget_conformer(reftarget, pdbentry, nentries, Cbuf);
  printf("Using %s (conformer %d) as the tentative reference\n\n",
	 reftarget, refconfm);
  found = 0;
  cp = Cbuf;
  for (iconfm = 0; iconfm < nconforms; iconfm++, cp++)
    if (cp->index == refconfm) {
      found = 1;
      break;
    }
  if (!found)
    abnormalexit("the reference target conformer is not found", 0);
  itarget = iconfm;

  /*
   * missing atoms in the tentative M are imputed by the centroid
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

  centerM(Mbuf, NULL, natoms, shift);

  /*
   * initial conformer-based weights. R1 norms as weights.
   */
  cp = Cbuf;
  for (i = 0; i < nconforms; i++)
    cp[i].rmsd = 1.0;

  /* 
   * initial S
   */
  for (i = 0; i < natoms; i++)
    diagS[i] = 1.0;

  logL = loglik(Mbuf, Xbuf, Cbuf, nconforms, natoms, diagS, spwork);
  rp = spwork;
  uwrmsd = rp[4];
  chisq = rp[5];
  printf("\n");
  printf("iter %03d\tloglik = %e (%.3e %.3e %.3e)\n",
	 iter=-1, rp[0], rp[1], rp[2], rp[3]);
  printf("\t\tchisq = %.3e   uwrmsd = %e\n", chisq, uwrmsd);
  /* set `previous values' in convergence_attained(). */
  spwork[0] = logL;		// actually spwork[0] == logL
  (void) convergence_attained(Rbuf, tbuf, nconforms, Mbuf, natoms, spwork, tol);
  printf("\n");

  /* pre-iteration OLS */
  /*
   * interdependence between R & t:
   * R depends on t.
   * t might depend on R and M but only on incomplete-case.
   * On complete-cases, tj == gj (centroid) by taking Rj M^T S^{-1} 1 == 0.
   * By this reason, the two functions are called in this order
   * at this time when the estimates tj are not determined.
   */
  for (iconfm = 0; iconfm < nconforms; iconfm++, cp++) {
    cp = &Cbuf[iconfm];
    estimate_t(cp->x, cp->o, diagS, natoms, cp->t);
    estimate_R(cp->x, cp->t, cp->o, Mbuf, diagS, natoms, cp->r);
    calculate_Y(cp->x, cp->r, cp->t, natoms, cp->y);
  }

  if (nconforms == 2 || nconforms == 1) {
    average_Y(Cbuf, nconforms, natoms, -1, options.r1weight, Mbuf, spwork);
    centerM(Mbuf, diagS, natoms, shift);
    for (iconfm = 0; iconfm < nconforms; iconfm++) {
      // translate Y
      cp = &Cbuf[iconfm];
      X = (Matrix3) cp->y;
      for (i = 0; i < natoms; i++) {
	X[i][0] -= shift[0];
	X[i][1] -= shift[1];
	X[i][2] -= shift[2];
      }
    }
    niters = 0;		/* never enter into the main loop */
  }

  logL = loglik(Mbuf, Xbuf, Cbuf, nconforms, natoms, diagS, spwork);
  rp = spwork;
  uwrmsd = rp[4];
  chisq = rp[5];
  printf("\n");
  printf("iter %03d\tloglik = %e (%.3e %.3e %.3e)\n",
	 iter=0, rp[0], rp[1], rp[2], rp[3]);
  printf("\t\tchisq = %.3e   uwrmsd = %e\n", chisq, uwrmsd);
  /* set `previous values' in convergence_attained(). */
  spwork[0] = logL;		// actually spwork[0] == logL
  (void) convergence_attained(Rbuf, tbuf, nconforms, Mbuf, natoms, spwork, tol);
  printf("\n");

  for (iter = 0; iter < niters; iter++) {

    if (options.ave_of_others) {
      average_Y(Cbuf, nconforms, natoms, -1, options.r1weight, Mbuf, spwork);
      centerM(Mbuf, diagS, natoms, shift);
      bulkUpdate = 1;
      /*
       * estimate R_j and t_j.
       *
       * bulkUpdate is nescessary when every M are shifted by the same amounts,
       * which appears to produce more accurate results
       * than stepwise update using centerM.
       */
      for (iconfm = 0; iconfm < nconforms; iconfm++) {
	cp = &Cbuf[iconfm];
	mp = Mbuf + 3 * natoms * 2;
	M = (Matrix3) mp;
	average_Y(Cbuf, nconforms, natoms, iconfm, options.r1weight, mp, spwork);
	if (bulkUpdate)		/* translate every M by `shift` */
	  for (iatom = 0; iatom < natoms; iatom++) {
	    M[iatom][0] -= shift[0];
	    M[iatom][1] -= shift[1];
	    M[iatom][2] -= shift[2];
	  }
	else
	  centerM(mp, diagS, natoms, 0);

	if (bulkUpdate) {
	  memcpy(cp->r+nconforms*9*2, cp->r, sizeof(Real)*9);
	  rp = cp->r; cp->r += nconforms * 9 * 2;
	  estimate_R(cp->x, cp->t, cp->o, mp, diagS, natoms, cp->r);
	  cp->r = rp;
	  memcpy(cp->t+nconforms*3*2, cp->t, sizeof(Real)*3);
	  rp = cp->t; cp->t += nconforms * 3 * 2;
	  estimate_t(cp->x, cp->o, diagS, natoms, cp->t);
	  cp->t = rp;
	} else {
	  estimate_R(cp->x, cp->t, cp->o, mp, diagS, natoms, cp->r);
	  estimate_t(cp->x, cp->o, diagS, natoms, cp->t);
	}
      }
      if (bulkUpdate) {
	memcpy(Rbuf, Rbuf+nconforms*9*2, sizeof(Real) * 9 * nconforms);
	memcpy(tbuf, tbuf+nconforms*3*2, sizeof(Real) * 3 * nconforms);
      }
    } else {		// if (!options.ave_of_others)
      average_Y(Cbuf, nconforms, natoms, -1, options.r1weight, Mbuf, spwork);
      centerM(Mbuf, diagS, natoms, shift);
      for (iconfm = 0; iconfm < nconforms; iconfm++) {
	cp = &Cbuf[iconfm];
	estimate_R(cp->x, cp->t, cp->o, Mbuf, diagS, natoms, cp->r);
	estimate_t(cp->x, cp->o, diagS, natoms, cp->t);
      }
    }

    for (iconfm = 0; iconfm < nconforms; iconfm++) {
      cp = &Cbuf[iconfm];
      calculate_Y(cp->x, cp->r, cp->t, natoms, cp->y);
    }

    /* update diagS */
    estimate_S(Cbuf, Mbuf, natoms, nconforms, diagS, anisoS, isoS, spwork);

    if (options.ols)		/* force W=I. i.e., OLS */
      for (i = 0; i < natoms; i++)
	diagS[i] = 1.0;


    logL = loglik(Mbuf, Ybuf, Cbuf, nconforms, natoms, diagS, spwork);

    if (!isfinite(logL)) {
      fflush(stdout);
      fprintf(stderr, "logL is either NaN or infinite. aborting...\n");
      abort();
    }

    rp = spwork;
    chisq = rp[5];
    uwrmsd = rp[4];
    printf("\n");
    printf("iter %03d\tloglik = %e (%.3e %.3e %.3e)\n",
	   iter+1, rp[0], rp[1], rp[2], rp[3]);
    /* CHECK CONVERGENCE */
    spwork[0] = logL;		// actually rp[0] == logL

    rp = spwork;
    i = convergence_attained(Rbuf, tbuf, nconforms, Mbuf, natoms, spwork, tol);
    printf("\t\tchisq = %.3e   rmsd/del = %e / %.3e\n", chisq, uwrmsd, rp[4]);
    printf("\t\tdL = %.3e  dM = %.3e  dR = %.3e  dt = %.3e\n",
	   rp[0], rp[1], rp[2], rp[3]);
    printf("\n");
    if (i != 0) break;
  }

  printf("Iteration finished:  %d / %d\n", iter+1, niters);
  rp = spwork;
  printf("\tdL: %.3e  dM: %.3e  dR: %.3e  dt: %.3e dr: %.3e\n",
	 rp[0], rp[1], rp[2], rp[3], rp[6]);

  qual[0] = iter + 1;
  qual[1] = niters;
  qual[3] = rp[0];
  qual[4] = rp[1];
  qual[5] = rp[4];

  rp = spwork;

  if (options.place_alongPA) {
    /* M <- M - shift, so that the centroid of M is on the origin.   */
    centerM(Mbuf, diagS, natoms, shift);
    /* M <- M V, so that M lays along its principal axes. */
    alignM(Mbuf, natoms, rotbuf);
    V = (Matrix3) rotbuf;
    /*
     * (Yj - shift) V = ((Xj - 1n x tj^T) Rj - 1n x shift^T) V
     * = (Xj - 1n x (tj + Rj shift)^T) Rj V
     */
    cp = Cbuf;
    for (iconfm = 0; iconfm < nconforms; iconfm++, cp++) {
      R = (Matrix3) (cp->r);
      /* tj <- tj + Rj s */ 
      for (i = 0; i < 3; i++) {
	s = 0.0;
	for (k = 0; k < 3; k++)
	  s += R[i][k] * shift[k];
	cp->t[i] += s;
      }
      /* Rj <- Rj V */
      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  s = 0.0;
	  for (k = 0; k < 3; k++)
	    s += R[i][k] * V[k][j];
	  rp[i*3+j] = s;
	}
      }
      for (i = 0; i < 9; i++)
	cp->r[i] = rp[i];
      /* update Yj */
      calculate_Y(cp->x, cp->r, cp->t, natoms, cp->y);
    }
  }

  estimate_S(Cbuf, Mbuf, natoms, nconforms, diagS, anisoS, isoS, spwork);

  logL = loglik(Mbuf, Ybuf, Cbuf, nconforms, natoms, diagS, spwork);
  rp = spwork;
  printf("\tloglik = %e (%.3e %.3e %.3e)\n", rp[0], rp[1], rp[2], rp[3]);
  printf("\tchisq = %.5e   rmsd = %e\n", rp[5], rp[4]);
  printf("\n");

  /* print conformers */
  printf("Conformers\n");
  printf("   id   gamma     nu     rmsd (weighted unweighted)\n");
  pe = pdbentry;
  k = 0;
  for (i = 0; i < nchains; i++) {
    cp = Cbuf + i;
    printf("  %05d   %5g   %5d   %10e  %10e   %s\n",
	   cp->index, *(cp->g), cp->natoms, cp->rmsd, cp->rmsd_uw, pe->file);
    if (++k >= pe->nchains) {
      pe++; k = 0;
    }
  }
  printf("\n\n");

  if (options.covfile[0]) {	/* write the sample covariance matrix */
    fp = fopen(options.covfile, "w");
    if (fp == NULL)
      fprintf(stderr, "cannot open a file %s for writing covar.\n",
	      options.covfile);
    else {
      rp = isoS;
      for (i = 0; i < natoms; i++) {
	for (j = 0; j < natoms; j++)
	  fprintf(fp, "%.16e ", rp[i*natoms + j]);
	fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }


  if (options.crdfile[0]) {
    fp = fopen(options.crdfile, "w");
    if (fp == NULL) {
      fprintf(stderr, "cannot open a file %s for writing the crds.\n",
	      options.crdfile);
      fp = stdout;
    }
    /*
      for (i = 0; i < nchains; i++) {
	arp = Bbuf + i;
	for (iconfm = 0; iconfm < arp->nconfms; iconfm++) {
	  cp = arp->confmbp + iconfm;
	  if (nconforms == 2) {		// pairwise superposition
	    if (iconfm == itarget) {
	      for (i = 0; i < natoms * 3; i++)
		cp->y[i] = cp->x[i];	// target, untransformed
	    } else {
	      calculate_Y(cp->x, arp->r, arp->t, natoms, cp->y);
	      Y = (Matrix3) cp->y;
	      R = (Matrix3) Bbuf[itarget].r;
	      rp = Bbuf[itarget].t;	// apply R^{-1} & -t of target confm
	      for (i = 0; i < natoms; i++) {
		sx[0] = Y[i][0]*R[0][0] + Y[i][1]*R[0][1] + Y[i][2]*R[0][2];
		sx[1] = Y[i][0]*R[1][0] + Y[i][1]*R[1][1] + Y[i][2]*R[1][2];
		sx[2] = Y[i][0]*R[2][0] + Y[i][1]*R[2][1] + Y[i][2]*R[2][2];
		Y[i][0] = sx[0] + rp[0];
		Y[i][1] = sx[1] + rp[1];
		Y[i][2] = sx[2] + rp[2];
	      }
	    }
	  } else			// multiple superposition
	    calculate_Y(cp->x, arp->r, arp->t, natoms, cp->y);
	  sprintf(buf, "conformer %d   gamma = %g", cp->id, *(cp->g));
	  print_conformer(cp->y, cp->o, chain_template, natoms,
			  1 , buf, fp);
	}
      }
    strcpy(buf, "====== mean structure");
    print_conformer(Mbuf, 0, chain_template, natoms, 1 , buf, fp);
    fclose(fp);
    */
  } else {
    fp = stdout;
  }

  fprintf(fp, "REMARK   1  ====== COORDINATES\n");
  if (nconforms > 2) {	/* multiple superposition onto the mean confm */
    for (ient = 0; ient < nentries; ient++) {
      pe = &pdbentry[ient];
      printf("REMARK   1  ====== %s\n", pe->file);
      //fp = fopen(pe->file, "r"); infile は pe で渡す
      read_trans_write_pdbm(fp, pe, Rbuf, tbuf, 1);
      //fclose(fp);
    }
  } else if (nconforms == 2) {		/* pairwise superposition. */
    /* M <- M R^T + t0^T, where R=I */
    dp = Cbuf[itarget].t;
    mp = Mbuf;
    for (i = 0; i < natoms; i++) {
      mp[0] += dp[0];  mp[1] += dp[1];  mp[2] += dp[2];
      mp += 3;
    }
    /* Y <- Y R^T + t0^T, where R=I */
    yp = Ybuf;
    for (iconfm = 0; iconfm < nconforms; iconfm++) {
      for (i = 0; i < natoms; i++) {
	yp[0] += dp[0];  yp[1] += dp[1];  yp[2] += dp[2];
	yp += 3;
      }
    }

    ic = 0;
    if (itarget == 0) ic = 1;

    /* t1 <- t1 - R1 t0 */
    R = (Matrix3) Cbuf[ic].r;
    dp = Cbuf[itarget].t;
    gp = Cbuf[ic].t;
    sx[0] = dp[0]*R[0][0] + dp[1]*R[0][1] + dp[2]*R[0][2];
    sx[1] = dp[0]*R[1][0] + dp[1]*R[1][1] + dp[2]*R[1][2];
    sx[2] = dp[0]*R[2][0] + dp[1]*R[2][1] + dp[2]*R[2][2];
    gp[0] -= sx[0]; gp[1] -= sx[1]; gp[2] -= sx[2];
    rp = spwork;
    memcpy(rp, Cbuf[itarget].t, sizeof(Real) * 3);
    memset(Cbuf[itarget].t, 0, sizeof(Real) * 3);
    memcpy(rp+3, Cbuf[itarget].r, sizeof(Real) * 9);
    identityvectmatrix(Cbuf[itarget].r, 3);

    if (nentries == 2)
      pe = &pdbentry[itarget];
    else
      pe = &pdbentry[0];	/* two confms are in the single file */
    fprintf(fp, "REMARK   1  ====== %s\n", pe->file);
  
    read_trans_write_pdbm(fp, pe, Rbuf, tbuf, 1);

    if (nentries == 2) {
      /* nentries = nchains = nconforms = 2, placed in two files. */
      pe = &pdbentry[ic];
      fprintf(fp, "REMARK   1  ====== %s\n", pe->file);

      read_trans_write_pdbm(fp, pe, Rbuf, tbuf, 1);
    }
    /* restore R & t for the target. */
    rp = spwork;
    memcpy(Cbuf[itarget].t, rp, sizeof(Real) * 3);
    memcpy(Cbuf[itarget].r, rp+3, sizeof(Real) * 9);
  } else {			/* nconforms == 1 */
    /* M <- M R^T + t0^T, where R=I */
    dp = Cbuf[itarget].t;
    mp = Mbuf;
    for (i = 0; i < natoms; i++) {
      mp[0] += dp[0];  mp[1] += dp[1];  mp[2] += dp[2];
      mp += 3;
    }
    /* Y <- Y R^T + t0^T, where R=I */
    yp = Ybuf;
    for (iconfm = 0; iconfm < nconforms; iconfm++) {
      for (i = 0; i < natoms; i++) {
	yp[0] += dp[0];  yp[1] += dp[1];  yp[2] += dp[2];
	yp += 3;
      }
    }

    identityvectmatrix(Cbuf[0].r, 3);
    memset(Cbuf[0].t, 0, sizeof(Real) * 3);
    pe = &pdbentry[0];
    read_trans_write_pdbm(fp, pe, Rbuf, tbuf, 1);
  }


  strcpy(buf, "====== mean structure");
  print_conformer(Mbuf, 0, chain_template, natoms, 0 /* extra */, buf, fp);

  if (fp != stdout)
    fclose(fp);



  if (binaryIO_binout && *binaryIO_binout && nconforms >= 1) {
    fd = open(dbinout, O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
      fprintf(stderr, "cannot open %s for writing\n", dbinout);
      exit(-1);
    }
    /* S */
    write(fd, isoS, sizeof(Real) * natoms * natoms);
    /* M */
    write(fd, Mbuf, sizeof(Real) * natoms * 3);
    /* Yi:  */
    write(fd, Ybuf, sizeof(Real) * natoms * 3 * nchains);
    /* (3n)^2 S */
    write(fd, anisoS, sizeof(Real) * natoms * natoms * 9);
    /* Ri */
    write(fd, Rbuf, sizeof(Real) * 9 * nconforms);
    /* ti */
    write(fd, tbuf, sizeof(Real) * 3 * nconforms);
    close(fd);

    fd = open(ibinout, O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
      fprintf(stderr, "cannot open %s for writing\n", ibinout);
      exit(-1);
    }
    write(fd, &natoms, sizeof(int));
    write(fd, &nconforms, sizeof(int));
    close(fd);
  }



  remove_workarea("superposwork");
  remove_workarea("mat0n3");

  msfit_finalize();

  printf("Normal termination.  %d/%d %.3e %.3e %.3e\n",
	 iter+1, niters, qual[3], qual[4], qual[5]);
  return 0;
}

int
convergence_attained(Real *R, Real *t, int nconfms, Real *M, int natoms,
		     Real del[], Real tol[])
{
  static int firstcall = 1;
  static Real loglikold, rmsdold;
  Real *Rold, *told, *Mold, rmax, tmax, mmax, ldiff, loglik, rmsd;
  Real tolr, tolt, tolm, toll;
  int i;

  loglik = del[0];
  rmsd = del[4];
  tolr = tol[0]; tolt = tol[1]; tolm = tol[2]; toll = tol[3];
  Rold = R + nconfms * 9;
  told = t + nconfms * 3;
  Mold = M + natoms * 3;
  if (firstcall) {
    firstcall = 0;
  } else {
    /* R */
    rmax = 0.0;
    for (i = 0; i < nconfms * 9; i++)
      if (fabs(R[i] - Rold[i]) > rmax) rmax = fabs(R[i] - Rold[i]);
    /* t */
    tmax = 0.0;
    for (i = 0; i < nconfms * 3; i++)
      if (fabs(t[i] - told[i]) > tmax) tmax = fabs(t[i] - told[i]);
    /* M */
    mmax = 0.0;
    for (i = 0; i < natoms * 3; i++)
      if (fabs(M[i] - Mold[i]) > mmax) mmax = fabs(M[i] - Mold[i]);
    /* loglik */
    ldiff = fabs(loglik - loglikold) / (fabs(loglik) + 1.0);
    /* */
    del[0] = ldiff;
    del[1] = mmax;
    del[2] = rmax;
    del[3] = tmax;
    del[4] = fabs(rmsd - rmsdold) / (fabs(rmsd) + 1.0);
  }
  /* rotate */
  memcpy(Rold, R, sizeof(Real) * nconfms * 9);
  memcpy(told, t, sizeof(Real) * nconfms * 3);
  memcpy(Mold, M, sizeof(Real) * natoms * 3);
  loglikold = loglik;
  rmsdold = rmsd;
  return rmax < tolr && tmax < tolt && mmax < tolm && ldiff < toll;
}
  


/*
 * partition_into_conformers()
 * calculates the partition coefficients (gammas) according to the occupancy
 * factors of atoms provided in frac[] and nfracs[].
 * the partition is a sequence of propotions which sum upto 1 or less.
 * nseqs: total num of sequences.
 * occ[]: list of the sequences of the occupancy factors.
 * noccs[k]: length of the sequence wrt the k-th atom in the list.
 * eps: to cope with rounding errors. occ < eps is considered as zero.
 */

int
partition_into_conformers(Real occ[], int noccs[], int nseqs, Real eps,
			  Real gamma[])
{  
  Real min, **prop, *rp;
  int i, n, *nprops, noccs_tot; 
  
  noccs_tot = 0;			/* total num of occupancy factors */
  for (i = 0; i < nseqs; i++) noccs_tot += noccs[i];
  n = 0;
  nprops = get_workarea("nprops", PDB_NPARTNS_MAX, sizeof(int));
  prop = get_workarea("partitioning", PDB_NPARTNS_MAX, sizeof(Real *));
  rp = get_workarea("prop", PDB_OCCFRAC_BUFSIZ, sizeof(Real));
  if (nprops == 0 || prop == 0 || rp == 0) {
    fprintf(stderr, "partition_into_conformers: working area is too small\n");
    return -1;
  }

  prop[0] = rp;
  nprops[0] = noccs[0];
  for (i = 1; i < nseqs; i++) {
    prop[i] = prop[i-1] + nprops[i-1];
    nprops[i] = noccs[i];
  }
  memcpy(prop[0], occ, noccs_tot * sizeof(Real));

  n = 0;
  while (1) {
    min = 2;
    for (i = 0; i < nseqs; i++)
      if (nprops[i] > 0 && *prop[i] < min) min = *prop[i];
    if (min > 1) break;		/* no props remain */
    gamma[n++] = min;		/* min among the heads of remaining seqs */
    for (i = 0; i < nseqs; i++)
      if (nprops[i] > 0) {	/* elements in the seq */
        *prop[i] -= min;	/* subtract min from the top of the seq */
        if (*prop[i] <= eps) {	/* consider it as zero */
          *prop[i] = 0.0;
          nprops[i]--;
          prop[i]++;
        }
      }
  }
  remove_workarea("prop");
  remove_workarea("partitioning");
  remove_workarea("nprops");

  return n;
}

/*
 * prescan a pdbm file for obtaining
 *  - partionings on each chain (to prepare partition coefs for conformers)
 *  - min & max values of landmark ids (to prepare id2ind[])
 */
int
prescan_pdbm(FILE *fp, pdbEntry *pe, int idcnt[], int id_lim,
	     pdbChain *chainbuf, int nchains_lim, int *minidp, int *maxidp,
	     int *work)
{
  int frag_id_min, frag_id_max, minid, maxid;
  int natoms, nsites, nchains, model_id, cid;
  int i, j, k, ichain, iatom, *id2ind;
  int top_block;
  Real *occp;
  pdbChain *chainp;
  pdbSite *sp;

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
      chainp->occ_nparts = 0;	/* num of partitionings for occ on an atom */
      chainp->occ_nfracs[chainp->occ_nparts] = 0;
      chainp->natoms = 0;
    }

    chainp->natoms += natoms;
    
#ifdef DEBUG
    if (natoms > 0) printf("ichn=%d, frac occ ", ichain);
#endif

    /*
     * list of partitions of occupancy on this chain
     * occ_nparts: total num of partitions
     * occ_nfracs[0..occ_nparts-1]: num of fractions (parts) for each partition
     * occ_frac[0..]: joined sequencies of fractions
     */
    for (iatom = 0; iatom < natoms; iatom++) {
    
#ifdef DEBUG
    /* print fractional occupancy factors */
      for (sp = pdb_atom[iatom].sp, j = 0; sp; sp = sp->next, j++)
	if (sp->occ == 1.0) break;
	else {
	  if (j == 0) printf("| %g(%d) ", sp->occ, pdb_atom[iatom].id);
	  else printf("%g(%d) ", sp->occ, pdb_atom[iatom].id);
	}
#endif
    
      sp = pdb_atom[iatom].sp;
      occp = chainp->occ_frac;
      k = 1;			// necessary for the first entry of partition
      for (i = 0; i < chainp->occ_nparts; i++) {
	k = 0;
	for (sp = pdb_atom[iatom].sp, j = 0; sp; sp = sp->next, j++)
	  if (sp->occ != occp[j]) { k = 1; break; }
	if (k == 0) break;
	occp += chainp->occ_nfracs[i];
      }
      /* if k == 0, iatom shares the i-th partition. */
      if (k != 0) {			// new partitioning
	for (sp = pdb_atom[iatom].sp, j = 0; sp; sp = sp->next, j++)
	  occp[j] = sp->occ;
	chainp->occ_nfracs[chainp->occ_nparts] = j;
	chainp->occ_firstid[chainp->occ_nparts] = pdb_atom[iatom].id;
	chainp->occ_nparts++;
      }
    }
#ifdef DEBUG
    printf(chainp->occ_nparts > 1 ? "\n" : "(none)\n");
#endif    
  }		/* while ((k = read_pdbm(...))) */

  if (pe->nchains_block > 0)
    pe->nchains_block = pe->nchains - top_block;
  else
    pe->nchains_block = 0;
  *minidp = minid; *maxidp = maxid;

  return 0;
}



int
msfit_init(int natoms, int nchains, int nconforms, int ngammas[], Real *gp)
{
  Real *rp, *xp, *yp;
  int ich, iatom, *op;
  Conformer *cp;
  Atom *ap;

  gbuf = (Real *) malloc(sizeof(Real) * nconforms);	// gamma
  Abuf = (Atom *) malloc(sizeof(Atom) * natoms * (nconforms + 1));
  chain_template = Abuf + natoms * nconforms;
  Cbuf = (Conformer *) malloc(sizeof(Conformer) * nconforms);
  Obuf = (int *) malloc(sizeof(int) * natoms * nconforms);
	/* twice for current and previous values, and one for exclusive-ave */
  tbuf = (Real *) malloc(sizeof(Real) * 3 * nchains  * (2+1));
  Rbuf = (Real *) malloc(sizeof(Real) * 3 * 3 * nchains  * (2+1));
  Mbuf = (Real *) malloc(sizeof(Real) * 3 * natoms * 3); /* curr, prev, excl */
						/* twice for X and Y */
  Xbuf = (Real *) malloc(sizeof(Real) * 3 * natoms * nconforms * 2);
  Ybuf = Xbuf + 3 * natoms * nconforms;
  if (!gbuf || !Abuf || !Cbuf || !Obuf || !tbuf || !Rbuf
      || !Mbuf || !Xbuf) {
    fprintf(stderr, "msfit_init(): cannot alloc memory.\n");
    return -1;
  }

  // これで Atom.id も Conform.occ も 0 に初期化される
  memset(Abuf, 0, sizeof(Atom) * natoms * (nconforms + 1));
  memset(Cbuf, 0, sizeof(Conformer) * nconforms);
  memset(Obuf, 0, sizeof(int) * natoms * nconforms);
  ap = Abuf;
  xp = rp = Xbuf;
  yp = Ybuf;
  op = Obuf;
  cp = Cbuf;
  memcpy(gbuf, gp, nconforms * sizeof(Real));
  gp = gbuf;
  for (ich = 0; ich < nchains; ich++, cp++) {
    if (ngammas[ich] > 1) {
      fprintf(stderr, "Multiple conformers belong chain %d.\n", ich);
      exit(-1);
    }
    cp->index= ich;
    cp->head = ap;
    cp->natoms = natoms;
    cp->x = xp;
    cp->y = yp;
    cp->r = Rbuf + (9 * ich);
    cp->t = tbuf + (3 * ich);
    cp->o = op;
    cp->g = gp;
    for (iatom = 0; iatom < natoms; iatom++) {
      ap->p = xp;
      ap->occ = op;
      ap++; op++;
      xp += 3;
    }
    *gp++ = 1.0;
    yp += 3 * natoms;
  }

  return 0;
}  


/*
 * work[0..id_lim] 
 */

int
read_coordinates(FILE *fp, pdbEntry *pe, int id2ind[], int id_lim,
		 Real occ_eps, int natoms_all, Conformer *confmpool, int work[])
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
    confmp = confmpool + pe->chainptr[ichain].index;
    for (iatom = 0; iatom < natoms; iatom++) {
      sp = pdb_atom[iatom].sp;
      d = sp->occ;
      ind = id2ind[pdb_atom[iatom].id];
      if (ind < 0) {
	fprintf(stderr, "Warning: negative index %d for atom #%d (id %d) in %s."
		" ignored\n", ind, iatom, pdb_atom[iatom].id, pe->file);
	continue;
      }
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
	d -= *(confmp->g);
	if (d < occ_eps) {
	  d = 0.0;
	  sp = sp->next;
	  if (sp) d = sp->occ;
	}
      } else {
	*ap->occ = 0;
	// ap->disordered = 1;
      }
    }
  }
  
  for (k = 0; k < pe->nchains; k++) {
    confmp = confmpool + pe->chainptr[k].index;
    n = 0;
    for (iatom = 0; iatom < natoms_all; iatom++) {
      ap = confmp->head + iatom;
      // if (*ap->occ == 0 && ap->disordered == 0)
      //   ap->disordered = 1;		/* missing as well on the chain */
      n += *ap->occ;
    }
    confmp->natoms = n;
  }

  return 0;
}


/*
 * occupancy_patterns() returns nlocsmiss.
 * locs[0..nlocsmiss-1] ... the indeces of missing atoms throughout the entire conformers
 * patn[0..notys-1] ... contains notys found missing patterns (bitwise OR of locs[])
 * ptrocc[0..notys-1] ... respective pointers to the arrays of occs
 * nobs[0..notys-1] ... num of observed atoms
 */
int
occupancy_patterns(Conformer confm[], int nconfms, int natoms, int nlocslim,
		   int notyslim, int locs[], int64_t patn[], int *ptrocc[],
		   int nobs[], int *notysp)
{
  int nlocsmiss, iconfm, iatom, found, i, numotys, *npatn;
  int64_t *otyof, bit;

  otyof = get_workarea("otyof", nconfms, sizeof(int64_t));
  npatn = get_workarea("npatn", nconfms, sizeof(int));
  if (otyof == 0 || npatn == 0) {
    fprintf(stderr, "occupancy_patterns: working area is too small\n");
    return -1;
  }

  /* num of locs at which the atom is missing at least in one confm */
  nlocsmiss = 0;
  
  /*
   * locs stores the indeces of missing atoms throughout the all conformers
   * the atom locs[k] is missing in the conformer i if the k-th lsb of otyof[i] is on.
   */

  for (iconfm = 0; iconfm < nconfms; iconfm++)
    otyof[iconfm] = 0;

  for (iatom = 0; iatom < natoms; iatom++) {
    found = 0;
    for (iconfm = 0; iconfm < nconfms; iconfm++)
      if (confm[iconfm].o[iatom] == 0) {
	if (!found) {
	  found = 1;
	  locs[nlocsmiss] = iatom;
	  bit = (int64_t)1 << nlocsmiss;
	  nlocsmiss++;
	  if (nlocsmiss > nlocslim) 
	    return -2;
	}
	otyof[iconfm] |= bit;
      }
  }

  /*
   * copy otys into patn[] removing duplicates between the conformers.
   * and give the pointer to that occ array of each missing pattern.
   */
  patn[0] = 0;			/* the first oty is always that of un-missing */
  nobs[0] = natoms;
  numotys = 1;
  ptrocc[0] = 0;
  for (iconfm = 0; iconfm < nconfms; iconfm++) {
    found = 0;
    for (i = 0; i < numotys; i++)
      if (otyof[iconfm] == patn[i]) {
	found = 1;
	break;
      }
    if (!found) {		// note that i == numotys
      patn[numotys] = otyof[iconfm];
      ptrocc[numotys] = confm[iconfm].o;
      nobs[numotys] = confm[iconfm].natoms;
      npatn[numotys] = 0;
     i = numotys++;
      if (numotys > notyslim)
	return -3;
    }
    confm[iconfm].oty = i;	// index of the oty-list
    npatn[i]++;
    if (i == 0 && ptrocc[0] == 0)
      ptrocc[0] = confm[iconfm].o;
  }
  /*
   * In the current implementation, it is assumed that
   * there exists at least one conformer that completes all landmark atoms.
   */
  if (ptrocc[0] == 0)
    fprintf(stderr, "Consider assigning ptrocc[0]\n"), exit(-1);

  if (verbose) {
    printf("\nMissing atoms and type of Sigma (oty):\n");
    printf("nlocs = %d: ", nlocsmiss);
    for (i = 0; i < nlocsmiss; i++)
      printf(" %d", locmiss[i]);
    printf("\nnum of conformers of each oty: ");
    for (i = 0; i < numotys; i++)
      printf(" (%d) 0x%lx x %d ", i, patn[i], npatn[i]);
    printf("\n");
  }

  remove_workarea("npatn");
  remove_workarea("otyof");

  *notysp = numotys;
      
  return nlocsmiss;
}


  /* diagS, isoS, and anisoS */
int
allocate_covar(int natoms)
{
  Sbuf = (Real *)
    malloc(sizeof(Real) * (natoms + (natoms*natoms) + (3*natoms*3*natoms)));
  if (Sbuf == 0) {
    fprintf(stderr, "Cannot allocate memory for covar.\n");
    return -1;
  }
  diagS = Sbuf;
  isoS = diagS + natoms;
  anisoS = isoS + natoms * natoms;
  return 0;
}


/*
 * OUTPUT:
 * if all covar[].unfactorized == 0
 *	covar[].logdet	log det Sj
 *	covar[].s	cholesky decomped Sj
 *	cover[].inv	Sj^{-1}
 *	return dodiag=0
 * dodiag or any covar[].unfactorized > 0
 *	covar[].logdet	log det diag Sj
 *	covar[].s	unchanged (Sj?)
 */
int
covar_factorize(struct covar covar[], int numotys, int natoms, int dodiag)
{
  int oty, nu, i, j;
  char buf[BUFSIZ];

  /*
  M = (Matrix3) mp;
  */
  if (!dodiag) {
    for (oty = 0; oty < numotys; oty++) {
      nu = covar[oty].nu;
      if (UseMKL)
	dpotrf_ ( "U", &nu, covar[oty].s, &nu, &i );
      else
	i = cholesky_decomp(covar[oty].s, nu);
      covar[oty].unfactorized = i;
      if (covar[oty].unfactorized == 0) {	/* successfully factorized */
	covar[oty].logdet = 0.0;
	for (i = 0; i < nu; i++)
	  covar[oty].logdet += log(covar[oty].s[i*(nu+1)]);
	covar[oty].logdet *= 2.0;		/* twice because cholesky */
	/* S^{-1} */
	if (UseMKL) {
	  memcpy(covar[oty].inv, covar[oty].s, sizeof(Real)*nu*nu);
	  dpotri_( "U", &nu, covar[oty].inv, &nu, &i );
	} else
	  i = cholesky_invert(covar[oty].s, covar[oty].inv, nu);
	if (i != 0)
	  sprintf(buf, "cholesky inv (%d)", i), abnormalexit(buf, 0);
	/* S^{-1} [M 1] */
	/*
	rp = work;
	extract_Q(mp, 3, natoms, 3, covar[oty].q, rp, nu, 1);
	for (i = 0; i < nu; i++)
	  rp[i+3*nu] = 1.0;
	ndim = 4;
	if (UseMKL)
	  dpotrs_("U", &nu, &ndim, covar[oty].s, &nu, rp, &nu, &i);
	else
	  i = cholesky_solve(covar[oty].s, rp, nu, ndim);
	if (i != 0)
	  sprintf(buf, "cholesky solve (%d)", i), abnormalexit(buf, 0);
	scatter_Q(rp, nu, 1, natoms, 4, covar[oty].q, covar[oty].sol, 4);
	*/
      } else {				/* factorization failed */
	covar[oty].logdet = -INFINITY;
	memset(covar[oty].inv, 0, sizeof(Real) * natoms * natoms);
	/* restore S remaining on the upper triangular portion */
	for (i = 0; i < nu; i++)
	  for (j = i+1; j < nu; j++)
	    covar[oty].s[j*nu + i] = covar[oty].s[i*nu + j];
	for (i = 0; i < natoms; i++)
	  if ((j = covar[oty].q[i]) >= 0)
	    covar[oty].s[j*nu + j] = covar[oty].diag[i];
      }
      
      /*
      if (verbose)
	printf("\t\tnondiag oty = %d, unfact = %d, logdet = %e\n",
	       oty, covar[oty].unfactorized, covar[oty].logdet);
      */
    }
    for (oty = 0; !dodiag && oty < numotys; oty++)
      if (covar[oty].unfactorized) {
	for (i = 0; i < numotys; i++)
	  printf("\t\toty = %d, logdet = %e, err = %d\n",
		 i, covar[i].logdet, covar[i].unfactorized);
	printf("\t\tSwitching to the diag method...\n"); 
	dodiag = 1;
      }
  }

  if (dodiag) {
    for (oty = 0; oty < numotys; oty++) {
      covar[oty].logdet = 0.0;
      //      zerovectmatrix(covar[oty].inv, natoms, natoms);
      for (i = 0; i < natoms; i++)
	if (covar[oty].o[i]) {
	  covar[oty].logdet += log(covar[oty].diag[i]);
	  /*
	    s = 1.0 / covar[oty].diag[i];
	  */
	  /* S^{-1} */
	  //	  covar[oty].inv[i*(natoms+1)] = s;
	  /* S^{-1} [M 1] */
	  /*
	  covar[oty].sol[i*4+0] = M[i][0] * s;
	  covar[oty].sol[i*4+1] = M[i][1] * s;
	  covar[oty].sol[i*4+2] = M[i][2] * s;
	  covar[oty].sol[i*4+3] = s;
	  */
	}
      /*
      printf("\t\tdiag oty = %d, unfact = %d, logdet = %e\n",
      oty, covar[oty].unfactorized, covar[oty].logdet);
      */
    }
  }
  return dodiag;
}


/*
int
covar_weighted(struct covar *covp, Real *mp, int natoms, int dodiag, Real *work)
{
  int nu, i, ndim;
  Real *rp, s;
  Matrix3 M;

  rp = work;
  M = (Matrix3) mp;
  nu = covar->nu;
  if (!dodiag) {
    extract_Q(mp, 3, natoms, 3, covp->q, rp, nu, 1);
    for (i = 0; i < nu; i++)
      rp[i+3*nu] = 1.0;
    ndim = 4;
    if (UseMKL)
      dpotrs_("U", &nu, &ndim, covp->s, &nu, rp, &nu, &i);
    else
      i = cholesky_solve(covp->s, rp, nu, ndim);
    if (i != 0)
      return i;
    scatter_Q(rp, nu, 1, natoms, 4, covp->q, covp->sol, 4);
  } else {
    for (i = 0; i < natoms; i++)
      if (covp->o[i]) {
	s = 1.0 / covp->diag[i];
	// S^{-1} [M 1]
	covp->sol[i*4+0] = M[i][0] * s;
	covp->sol[i*4+1] = M[i][1] * s;
	covp->sol[i*4+2] = M[i][2] * s;
	covp->sol[i*4+3] = s;
      } else
	covp->sol[i*4+0] = covp->sol[i*4+1] = covp->sol[i*4+2] = covp->sol[i*4+3] = 0.0;
  }
  return 0;
}
*/

void
print_conformer(Real *yp, int *occ, Atom *templ, int natoms, int extra,
		char *remark, FILE *fp)
{
  int iatom, fwidth, fprec;
  Matrix3 Y;
  Atom *ap;
  char rec[BUFSIZ], buf[BUFSIZ];

  if (extra) {
    fwidth = PDBM_EXTRA_x_LEN - 1;
    fprec = fwidth - 7;
  }

  fprintf(fp, "REMARK   1  %s\n", remark);
  Y = (Matrix3) yp;
  ap = templ;
  for (iatom = 0; iatom < natoms; iatom++, ap++) {
    if (occ == 0 || occ[iatom] != 0) {
      sprintf(rec, "ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f%6.2f",
	      ap->id, ap->name, ap->resname, ap->resid,
	      Y[iatom][0], Y[iatom][1], Y[iatom][2], 1.0);
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
