#ifndef _MATRIX_H
#define _MATRIX_H

#define det3x3(A) (A[0][0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1] * (A[1][0]*A[2][2] - A[1][2]*A[2][0]) + A[0][2] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]))

/* matrix.c */
extern void zerovectmatrix(Real *a, int m, int n);
extern void identityvectmatrix(Real *a, int n);
extern void printMatrix3(char *name, Matrix3 A, int m);
extern Real * multvectmatrix(Real *a, Real *b, Real *c, int m, int l, int n, int AT, int BT);
extern void prepare_Q(int q[], int o[], int ndim);
extern void extract_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb);
extern void scatter_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb);
extern void extract_Q(Real *A, int nca, int n, int m, int q[], Real *B, int ncb, int trb);
extern void scatter_Q(Real *A, int nca, int tra, int n, int m, int q[], Real *B, int ncb);

/* chol.c */
extern int cholesky_decomp(Real *ap, int n);
extern int cholesky_solve(Real *lp, Real *b, int n, int m);
extern int cholesky_invert(Real *lp, Real *xt, int n);

/* svdgl.c */
extern int svdgl(Real *A, Real *U, Real *V, int m, int n, Real *w);
extern int seiggl(Real *A, Real *V, int n, Real *w);

#endif
