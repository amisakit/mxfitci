#include <math.h>
#include <stdio.h>
#include "atypes.h"


#define CHOLESKY_EPS 0.0

/*
 * Cholesky decomposition, A = L L^T.
 * ap (input):  a symmetric PD matrix A of n x n
 * ap (output): L (lower and diagonal. upper triangular unchainged)
 */
int
cholesky_decomp(Real *ap, int n)
{
  Real d, s;
  int i, j, k;

  for (j = 0; j < n; j++) {
    s = ap[j*n+j];
    for (k = 0; k < j; k++)
      s -= ap[j*n+k] * ap[j*n+k];
    if (s <= CHOLESKY_EPS)
      return j + 1;
    d = sqrt(s);
    ap[j*n+j] = d;
    for (i = j + 1; i < n; i++) {
      s = ap[i*n+j];
      for (k = 0; k < j; k++)
	s -= ap[i*n+k] * ap[j*n+k];
      ap[i*n+j] = s / d;
    }
  }
  return 0;
}


/*
 * solve L (L^T X^T) = B^T
 */
int
cholesky_solve(Real *lp, Real *b, int n, int m)
{
  Real *x, *y, s;
  int i, j, k;

  for (i = 0; i < m; i++) {
    y = b;
    for (j = 0; j < n; j++) {
      y[j] = b[j];
      s = 0.0;
      for (k = 0; k < j; k++)
	s += lp[j*n+k] * y[k];
      y[j] -= s;
      y[j] /= lp[j*n+j];
    }
    x = b;
    for (j = n-1; j >= 0; j--) {
      x[j] = y[j];
      s = 0.0;
      for (k = j+1; k < n; k++)
	s += lp[k*n+j] * x[k];
      x[j] -= s;
      x[j] /= lp[j*n+j];
    }
    b += n;
  }
  return 0;
}


/*
 * X^T = (L L^T)^{-1} == X
 */
int
cholesky_invert(Real *lp, Real *xt, int n)
{
  Real s;
  int i, j, k;

  /* (L^{-1})^T */
  for (i = n - 1; i >= 0; i--) {
    s = 1.0;
    for (j = i; j < n; j++) {
      for (k = i; k < j; k++)
	s -= lp[j*n+k] * xt[i*n+k];
      xt[i*n+j] = s / lp[j*n+j];
      s = 0.0;
    }
  }
  /* (L^{-1})^T L^{-1} */
  for (i = 0; i < n; i++) {
    for (j = n - 1; j >= i; j--) {
      s = 0.0;
      for (k = j; k < n; k++)
	s += xt[k+i*n] * xt[k+j*n];
      xt[i+j*n] = s;
    }
  }
  for (i = 0; i < n; i++)
    for (j = i; j < n; j++)
      xt[i*n+j] = xt[i+j*n];
  return 0;
}



/*
Matrix Inversion Using Cholesky Decomposition
Aravindh Krishnamoorthy, Deepak Menon
int
cholesky_inv(Real *lp, int n, Real *xp)
{
  Real d, s;
  int i, j, k;

  for (j = n-1; j >= 0; j--) {
    d = 1.0 / lp[j*n+j];
    for (i = j; i >= 0; i--) {
      s = 0.0;
      for (k = i+1; k < n; k++)
	s += lp[k*n+i] * xp[k*n+j];
      xp[i*n+j] = (d - s) / lp[i*n+i];
      xp[j*n+i] = xp[i*n+j];
      d = 0.0;
    }
  }
  return 0;
}
*/

  
/*

main()
{
  Real A[] = {
    4.16, -3.12,  0.56, -0.10,
    -3.12, 5.03, -0.83,  1.18,
    0.56, -0.83,  0.76,  0.34,
    -0.10, 1.18,  0.34,  1.18};
  Real b[] = {1,2,3,4,7,8,9,10};
  Real C[100], D[100], s;
  
  int i, j, k, n, m;

  n = 4; m = 2;
  for (i = 0; i < n*n; i++) C[i] = A[i];
  cholesky_decomp(C, n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%e ", C[i*n+j]);
    printf("\n");
  }
  printf("\n");
  (void) cholesky_solve(C, b, n, m);
  for (j = 0; j < n*m; j++)
    printf("%e ", b[j]);
  printf("\n");
  printf("\n");
  for (k = 0; k < m; k++)
  for (i = 0; i < n; i++) {
    s = 0.0;
    for (j = 0; j < n; j++)
      s += A[i*n+j]*b[j+k*n];
    printf("%e ", s);
  }
  printf("\n");
  printf("\n");
    
  (void) cholesky_inv(C, n, D);
  for (i = 0; i < n; i++)
    for (j = i+1; j < n; j++)
      C[j*n+i] = C[i*n+j];
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      s = 0.0;
      for (k = 0; k < n; k++)
	s += A[i*n+k] * D[k*n+j];
      printf("%e ", s);
    }
    printf("\n");
  }
  printf("\n");


  for (i = 0; i < n*n; i++) C[i] = A[i];
  cholesky_decomp(C, n);
  cholesky_invert(C, D, n);
  for (i = 0; i < n; i++)
    for (j = i+1; j < n; j++)
      C[j*n+i] = C[i*n+j];
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      s = 0.0;
      for (k = 0; k < n; k++)
	s += A[i*n+k] * D[k*n+j];
      printf("%e ", s);
    }
    printf("\n");
  }
  printf("\n");


  for (i = 0; i < n*n; i++) C[i] = A[i];
  dpotrf_ ( "U", &n, C, &n, &i );
  dpotri_ ( "U", &n, C, &n, &i );
  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++)
      C[j*n+i] = C[i*n+j];
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      s = 0.0;
      for (k = 0; k < n; k++)
	s += A[i*n+k] * C[k*n+j];
      printf("%e ", s);
    }
    printf("\n");
  }
}  

*/
