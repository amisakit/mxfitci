#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "atypes.h"
#include "matrix.h"

/*
 * *vectmatrix deals with vect(A) = vec(A^T) = [a11, a12, ..., a1n, a21, ... ]
 */

void zerovectmatrix(Real *a, int m, int n)
{
  memset(a, 0, sizeof(Real) * m * n);
}

void identityvectmatrix(Real *a, int n)
{
  int i;
  memset(a, 0, sizeof(Real) * n * n);
  for (i = 0; i < n-1; i++) {
    *a = 1.0;
    a += n + 1;
  }
  *a = 1.0;
}

void printMatrix3(char *name, Matrix3 A, int m)
{
  int i, j;
  printf("%s:\n", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < 3; j++)
      printf("%8.2f", A[i][j]);
    printf("\n");
  }
  printf("\n");
}

/*
 * returns C = A^(T*(1-AT)) (m x l) * B^(T*(1-BT)) (l x n) 
 */ 
Real *
multvectmatrix(Real *a, Real *b, Real *c, int m, int l, int n, int AT, int BT)
{
  Real *ap, *a0p, *bp, *b0p, *cp, s;
  int i, j, k;

  switch (AT*2+BT) {
  case 0:
    cp = c;
    a0p = a;
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
	s = 0.0;
	ap = a0p;			/* &A[i][0]; */
	bp = b + j;			/* &B[0][j]; */
	for (k = 0; k < l; k++, bp += n)
	  s += *ap++ * *bp;		/* s += A[i][k] * B[k][j]; */
	*cp++ = s;			/* C[i][j] = s; */
      }
      a0p += l;
    }
    break;

  case 1:
    cp = c;
    a0p = a;
    for (i = 0; i < m; i++) {
      b0p = b;
      for (j = 0; j < n; j++) {
	s = 0.0;
	ap = a0p;			/* &A[i][0]; */
	bp = b0p;			/* &B[j][0]; */
	for (k = 0; k < l; k++)
	  s += *ap++ * *bp++;		/* s += A[i][k] * B[j][k]; */
	*cp++ = s;			/* C[i][j] = s; */
	b0p += l;
      }
      a0p += l;
    }
    break;

  case 2:
    cp = c;
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
	s = 0.0;
	ap = a + i;			/* &A[0][i]; */
	bp = b + j;			/* &B[0][j]; */
	for (k = 0; k < l; k++, ap += m, bp += n)
	  s += *ap * *bp;		/* s += A[k][i] * B[k][j]; */
	*cp++ = s;			/* C[i][j] = s; */
      }
    }
    break;

  case 3:
    cp = c;
    for (i = 0; i < m; i++) {
      b0p = b;
      for (j = 0; j < n; j++) {
	s = 0.0;
	ap = a + i;			/* &A[0][i]; */
	bp = b0p;			/* &B[j][0]; */
	for (k = 0; k < l; k++, ap += m)
	  s += *ap * *bp++;		/* s += A[k][i] * B[j][k]; */
	*cp++ = s;			/* C[i][j] = s; */
	b0p += l;
      }
    }
    break;
  }

  return c;
}




/* perm.c */


/* array that works as the sub-permutation Q */
void prepare_Q(int q[], int o[], int ndim)
{
  int *op, i, j;

  op = o;
  j = 0;
  for (i = 0; i < ndim; i++)
    if (*op++)
      q[i] = j++;
    else
      q[i] = -1;
  return;
}
  
/*
 * extract elements of A into B st B = Q A Q^T. 
 * A : an n x n matrix on the 2-D array of nca columns
 * B : the nu x nu matrix on the 2-D array of ncb columns
 * q : an array realizing a sub-permutation Q
 */
void extract_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb)
{
  int i, j;

  for (i = 0; i < n; i++)
    if (q[i] >= 0)
      for (j = 0; j < n; j++)
	if (q[j] >= 0)
	  B[q[i]*ncb + q[j]] = A[i*nca + j];
  return;
}

/*
 * distribute elements of A into B st B = Q^T A Q. 
 */
void scatter_QQ(Real *A, int nca, int n, int q[], Real *B, int ncb)
{
  int i, j;

  for (i = n-1; i >= 0; i--)
    if (q[i] >= 0)
      for (j = n-1; j >= 0; j--)
	if (q[j] >= 0)
	  B[i*ncb + j] = A[q[i]*nca + q[j]];
  return;
}  

/*
 * extract rows of A into B st B = Q A. 
 * A : n x m matrix on the 2-D array of nca columns
 * B : nu x m matrix on the 2-D array of ncb columns
 *		if B should be stored in column-major order then set trb
 * q : an array realizing a sub-permutation Q
 */
void extract_Q(Real *A, int nca, int n, int m, int q[], Real *B, int ncb, int trb)
{
  int i, j;

  if (trb) {
    for (i = 0; i < n; i++)
      if (q[i] >= 0)
	for (j = 0; j < m; j++)
	  B[q[i] + j*ncb] = A[i*nca + j];
  } else {
    for (i = 0; i < n; i++)
      if (q[i] >= 0)
	for (j = 0; j < m; j++)
	  B[q[i]*ncb + j] = A[i*nca + j];
  }
  return;
}

/*
 * distribute rows of A into B st B = Q^T A. 
 * A : nu x m matrix on the 2-D array of nca columns
 *		if A is given in column-major order then set tra
 * B : n x m matrix on the 2-D array of ncb columns
 */
void scatter_Q(Real *A, int nca, int tra, int n, int m, int q[], Real *B, int ncb)
{
  int i, j;

  if (tra) {
    for (i = n-1; i >= 0; i--)
      if (q[i] >= 0)
	for (j = 0; j < m; j++)
	  B[i*ncb + j] = A[q[i] + j*nca];
  } else {
    for (i = n-1; i >= 0; i--)
      if (q[i] >= 0)
	for (j = 0; j < m; j++)
	  B[i*ncb + j] = A[q[i]*nca + j];
  }
  return;
}  

