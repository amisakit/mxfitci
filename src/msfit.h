#ifndef _MSFIT_H
#define _MSFIT_H
#include <stdio.h>
#include "atypes.h"
#include "pdb.h"
#include "molec.h"

#define COV_NORMAL 0
#define COV_DIAGONAL 1
#define COV_SHRINK 2
#define COV_TRUNCATE 4
#define COV_SD_MEST 8
#define COV_COR_KENDALL 16
#define COV_NOTGIVEN (-1)

struct covar {
  Real *s;      // packed as nu x nu when incomplete
  Real *diag;   // unpacked with gaps when incomplete
  Real *inv;    // packed as nu x nu
  // Real *sol;    // S^{-1} [M 1]
  Real logdet;
  //int soldiag;  // the sol above is that of diagonal eqns or not
  int solved;
  int nu;       // natoms of this oty
  int *q;       // sub-permutation
  int *o;       // a pointer to the array of occs of this oty
  int unfactorized;     // s unfactorized and inv not contain S^{-1}
};

struct arguments {
  int verbose;
  int ols;
  int ave_of_others;
  Real alpha, tau;
  Real lgkmax;
  int r1weight;
  int reduceCC;
  int initweight;
  int place_alongPA;
  int covty;
  int impute0;
  Real tol[5];
  int impute0samplecovar;
  char covfile[BUFSIZ];
  char crdfile[BUFSIZ];
  char wcovfile[BUFSIZ];
  char reftarget[BUFSIZ];
};


/* mmfitc.c */
extern void abnormalexit(char *message, int code);

/* args.c */
extern int parse_arguments(int argc, char *argv[],
			   pdbEntry pdbentry[], struct arguments *opts);

/* perm.c */
extern void prepare_Q(int q[], int o[], int ndim);
extern void extract_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb);
extern void scatter_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb);
extern void extract_Q(Real *A, int nca, int n, int m, int q[], Real *B, int ncb, int trb);
extern void scatter_Q(Real *A, int nca, int tra, int n, int m, int q[], Real *B, int ncb);

/* superimposems.c */
extern int
centerM(Real *mp, Real *diagS, int natoms, Real *t);
extern int
estimate_t(Real *xp, int *op, Real *diagS, int natoms, Real *t);
extern int
alignM(Real *mp, int natoms, Real *vout);
extern int
estimate_R(Real *xp, Real *t, int *op, Real *mp, Real *diagS, int natoms, Real *rp);
extern int
average_Y(Conformer *cp, int nconfms, int natoms, int iexcl,
          int r1weight, Real *mp, Real *work);
extern void
calculate_Y(Real *vectX, Real *vectR, Real *t, int natoms, Real *vectY);


/* covarms.c */
extern int
estimate_S(Conformer *cp, Real *mp, int natoms, int nconfms,
	   Real *diagS, Real *saniso, Real *siso, Real *work);
#endif
