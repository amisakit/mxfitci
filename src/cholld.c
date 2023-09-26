#include <quadmath.h>
#include <math.h>
#include <stdio.h>
#include "atypes.h"

long double wk[1000000];

#define CHOLESKY_EPS 0.0

/*
 * Cholesky decomposition, A = L L^T.
 * ap (input):  a symmetric PD matrix A of n x n
 * ap (output): L (lower and diagonal. upper triangular unchainged)
 */
int
cholesky_decomp(Real *arg_ap, int n)
{
  long double d, s, *ap;
  int i, j, k;

  ap = wk;
  for (i = 0; i < n*n; i++)
    ap[i] = arg_ap[i];

  for (j = 0; j < n; j++) {
    s = ap[j*n+j];
    for (k = 0; k < j; k++)
      s -= ap[j*n+k] * ap[j*n+k];
    if (s <= CHOLESKY_EPS)
      return j + 1;
    d = sqrtq(s);
    ap[j*n+j] = d;
    for (i = j + 1; i < n; i++) {
      s = ap[i*n+j];
      for (k = 0; k < j; k++)
	s -= ap[i*n+k] * ap[j*n+k];
      ap[i*n+j] = s / d;
    }
  }

  for (i = 0; i < n*n; i++)
    arg_ap[i] = ap[i];
  return 0;
}


/*
 * solve L (L^T X^T) = B^T
 */
int
cholesky_solve(Real *arg_lp, Real *arg_b, int n, int m)
{
  long double *x, *y, s, *lp, *b;
  int i, j, k;

  lp = wk;
  for (i = 0; i < n*n; i++)
    lp[i] = arg_lp[i];
  b = wk + n * m;
  for (i = 0; i < n*m; i++)
    b[i] = arg_b[i];

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

  b = wk + n * m;
  for (i = 0; i < n*m; i++)
    arg_b[i] = b[i];
  return 0;
}


/*
 * X^T = (L L^T)^{-1} == X
 */
int
cholesky_invert(Real *arg_lp, Real *arg_xt, int n)
{
  long double s, *lp, *xt;
  int i, j, k;

  lp = wk;
  for (i = 0; i < n*n; i++)
    lp[i] = arg_lp[i];
  xt = wk + n*n;
  for (i = 0; i < n*n; i++)
    xt[i] = arg_xt[i];
  
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

  xt = wk + n*n;
  for (i = 0; i < n*n; i++)
    arg_xt[i] = xt[i];
  return 0;
}

