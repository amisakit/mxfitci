/*
 * svd and symmetric eigen-decomposition
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "atypes.h"
#include "matrix.h"

/*
 * Householder vector.
 * Golub & van Loan: Algorithm 5.1.1
 * x and v can be the same region.
 */
static Real
house(Real x[], int n, Real v[])
{
  Real sig, beta, mu;
  int i;

  // normalize x[] and stores in v[0..n-1]
  sig = 0.0;
  for (i = 0; i < n; i++)
    sig += x[i] * x[i];
  if (sig <= 0.0)
    return 0.0;
  sig = sqrt(sig);
  for (i = 0; i < n; i++)
    v[i] = x[i] / sig;

  if (v[0] == 1.0) {
    beta = 0.0;
  } else {
    mu = 1.0;
    // more precise? than computing as sig = 1-v[0]^2 or sig = (1-v[0])*(1+v[0])
    sig = 0.0;
    for (i = 1; i < n; i++)
      sig += v[i] * v[i];
    if (v[0] <= 0.0)
      v[0] = v[0] - mu;
    else
      v[0] = -sig / (v[0] + mu);
    beta = 2.0 * v[0] * v[0] / (sig + v[0] * v[0]);
    for (i = 1; i < n; i++)
      v[i] /= v[0];
  }
  v[0] = 1.0;
  return beta;
}

/*
 * Bidiagonalization.
 * Golub & van Loan: Sec 5.4.8, Algorithm 5.4.2
 */
static void
bidiag(Real *A, int m, int n, Real beta[], Real v[], Real *H)
{
  int jj, i, j, k;
  Real bet, s, *x = v, *betaU = beta, *betaV = beta + m;
  Real *uv;

  // where to store the essential parts of Householder vectors
  if (H) uv = H; else uv = A;
  for (jj = 0; jj < n; jj++) {
    for (i = jj; i < m; i++) x[i] = A[i*n+jj];
    bet = house(x + jj, m - jj, v + jj);
    betaU[jj] = bet;
    for (j = jj; j < n; j++) {
      s = 0.0;
      for (k = jj; k < m; k++)
        s += v[k] * A[k*n+j];
      s *= bet;
      for (i = jj; i < m; i++)
          A[i*n+j] -= v[i] * s;
    }
    for (i = jj+1; i < m; i++)
      uv[i*n+jj] = v[i];

    if (jj < n-2) {
      for (j = jj+1; j < n; j++) x[j] = A[jj*n+j];
      bet = house(x+jj+1, n - jj - 1, v+jj+1);
      betaV[jj] = bet;
      for (i = jj; i < m; i++) {
        s = 0.0;
        for (k = jj+1; k < n; k++)
          s += A[i*n+k] * v[k];
        s *= bet;
        for (j = jj+1; j < n; j++)
          A[i*n+j] -= s * v[j];
      }
      for (j = jj+2; j < n; j++)
        uv[jj*n+j] = v[j];
    }
  }
}

/*
 * Tridiagonalization of a symmetric matrix A.
 * Golub & van Loan Alg 8.3.1
 * v is a work area of 3n.
 */
static void
tridiag(Real *A, int n, Real beta[], Real v[], Real *H)
{
  int i, j, k;
  Real bet, s, *x, *w, *p;
  Real *q;

  x = v;
  w = v + n;
  p = w + n;
  // where to store the essential parts of Householder vectors
  if (H) q = H; else q = A;
  for (k = 0; k < n-2; k++) {
    for (i = k+1; i < n; i++) x[i] = A[i*n+k];
    beta[k] = bet = house(x + k + 1, n - k - 1, v + k + 1);
    for (i = k+1; i < n; i++) {
      s = 0.0;
      for (j = k+1; j < n; j++)
        s += A[i*n+j] * v[j];
      s *= bet;
      p[i] = s;
    }
    s = 0.0;
    for (j = k+1; j < n; j++)
      s += p[j] * v[j];
    s *= bet / 2;
    for (i = k+1; i < n; i++)
      w[i] = p[i] - s * v[i];
    s = 0.0;
    for (j = k+1; j < n; j++)
      s += A[j*n+k] * A[j*n+k];
    s = sqrt(s);
    A[(k+1)*n+k] = A[k*n+(k+1)] = s;
    for (i = k+1; i < n; i++)
      for (j = k+1; j < n; j++)
        A[i*n+j] -= v[i] * w[j] + w[i] * v[j];
    for (i = k+2; i < n; i++)
      q[k*n+i] = v[i];
  }
  return;
}

/*
 * explicitly form an m x m orthogonal matrix U
 * as a product of Householder matrices U_1 U_2 ... U_r
 * stored the lower triangular of A and beta.
 * Golub & van Loan Sec 5.1.6.
 */
static void
houseMatrixU(Real *A, Real beta[], int m, int n, int r, int init, Real *U, Real v[])
{
  int i, j, k, jj, from;
  Real s;

  if (init) {
    for (i = 0; i < m; i++) {
      for (j = 0; j < m; j++)
        U[i*m+j] = 0;
      U[i*m+i] = 1.0;
    }
  }
  for (jj = r - 1; jj >= 0; jj--) {
    for (i = jj+1; i < m; i++)
      v[i] = A[i*n+jj];
    v[jj] = 1;
    from = 0;
    if (init) from = jj;
    for (j = from; j < m; j++) {
      s = 0.0;
      for (k = jj; k < m; k++)
	s += v[k] * U[k*m+j];
      s *= beta[jj];
      for (i = jj; i < m; i++)
	U[i*m+j] -= v[i] * s;
    }
  }
  return;
}

static void
houseMatrixV(Real *A, Real beta[], int n, int r, int init, Real *V, Real v[])
{
  int i, j, k, jj, from;
  Real s;

  if (init) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
        V[i*n+j] = 0;
      V[i*n+i] = 1.0;
    }
  }
  for (jj = r - 1; jj >= 0; jj--) {
    for (j = jj+2; j < n; j++)
      v[j] = A[jj*n+j];
    v[jj+1] = 1.0;
    from = 0;
    if (init) from = jj + 1;
    for (j = from; j < n; j++) {
      s = 0.0;
      for (k = jj+1; k < n; k++)
	s += v[k] * V[k*n+j];
      s *= beta[jj];
      for (i = jj+1; i < n; i++)
	V[i*n+j] -= v[i] * s;
    }
  }
  return;
}


/*
 * Givens rotation. Golub & van Loan Sec 5.1.8.
 * https://en.wikipedia.org/wiki/Givens_rotation
 * ``we require r to be positive'' (r^2 = a^2 + b^2)
 */
static void givens(Real a, Real b, Real *c, Real *s)
{
  Real tau;

  if (b == 0.0) {
    *c = copysign(1.0, a); *s = 0.0;
  } else if (fabs(b) > fabs(a)) {
    tau = a / b;
    *s = -1.0 / copysign(sqrt(1.0 + tau * tau), b);
    *c = -*s * tau;
  } else {
    tau = b / a;
    *c = 1.0 / copysign(sqrt(1.0 + tau * tau), a);
    *s = -*c * tau;
  }
}

// annihilate a rather than b

static void givens2(Real a, Real b, Real *c, Real *s)
{
  Real tau;

  if (a == 0.0) {
    *c = copysign(1.0, b); *s = 0.0;
  } else if (fabs(b) > fabs(a)) {
    tau = a / b;
    *c = 1.0 / sqrt(1.0 + tau * tau);
    *s = *c * tau;
  } else {
    tau = b / a;
    *s = 1.0 / sqrt(1.0 + tau * tau);
    *c = *s * tau;
  }
}


/*
 * a stack storing Givens rotations
 */
#define GR_COL 0		/* apply from right */
#define GR_ROW 1		/* apply from left  */
#define GR_BUFSIZ 50000

struct givensrot {
  int p;	/* offset */
  int n;	/* size */
  int dir;	/* GR_COL | GR_ROW */
  Real c, s;
  int i, k;
} grbuf[GR_BUFSIZ];

static int ngrots;

static void
grs_init()
{
  ngrots = 0;
}

static void grs_push(int i, int k, Real c, Real s, int p, int n, int dir)
{
  struct givensrot *rot;
  if (ngrots > GR_BUFSIZ) {
    fprintf(stderr, "grbuf exceeded its limit.\n");
    exit(-1);
  }
  rot = grbuf + ngrots;
  rot->i = i;
  rot->k = k;
  rot->p = p;
  rot->n = n;
  rot->c = c;
  rot->s = s;
  rot->dir = dir;
  ngrots++;
}

static struct givensrot *grs_pop()
{
  if (ngrots <= 0)
    return NULL;
  else
    return &(grbuf[--ngrots]);
}


static void
givensMatricesUV(Real *U, Real *V, int m, int n)
{
  Real *A, c, s, tau1, tau2;
  struct givensrot *q;
  int i, j, k, nc;

  while ((q = grs_pop())) {
    if (q->dir == GR_COL) {
      A = V; nc = n;
    } else {
      A = U; nc = m;
    }
    if (A == 0) continue;
    i = q->i;
    k = q->k;
    c = q->c;
    s = q->s;
    for (j = q->p; j < q->p + q->n; j++) {
      tau1 = A[(i+q->p)*nc+j];
      tau2 = A[(k+q->p)*nc+j];
      A[(i+q->p)*nc+j] = c * tau1 + s * tau2;
      A[(k+q->p)*nc+j] = -s * tau1 +  c * tau2;
    }
  }
  return;
}


/*
 * Golub & van Loan Algorithm 8.6.1
 * d[i] and f[i] are the diagonal and superdiagonal of a bidiagonal matrix B,
 * having no zeros.
 * golubkahn() returns B <- U^T B V (overwritten).
 * U and V are pushed onto the stack.
 */
static void
golubkahn(Real d[], Real f[], int offset, int n)
{
  Real a, mu, tau1, tau2, c, s, u, v, y, z;
  int k;

  a = 0.0;
  if (n > 2) a = f[n-3] * f[n-3];
  a = (a + d[n-2] * d[n-2] - d[n-1] * d[n-1] - f[n-2] * f[n-2]) / 2.0;
  mu = d[n-1] * d[n-1] + f[n-2] * f[n-2] + a
    - copysign(1.0, a) * sqrt(a * a + d[n-2]*d[n-2]*f[n-2]*f[n-2]);
  y = d[0] * d[0] - mu;
  z = d[0] * f[0];
  for (k = 0; k < n-1; k++) {
    givens(y, z, &c, &s);
    grs_push(k, k+1, c, s, offset, n, GR_COL);
    if (k > 0)
      f[k-1] = c * f[k-1] - s * u;
    tau1 = d[k]; tau2 = f[k];
    d[k] = c * tau1 - s * tau2;
    f[k] = s * tau1 + c * tau2;
    v = -s * d[k+1];
    d[k+1] = c * d[k+1];
    y = d[k];
    z = v;
    givens(y, z, &c, &s);
    grs_push(k, k+1, c, s, offset, n, GR_ROW);
    d[k] = c * d[k] - s * v;
    tau1 = f[k]; tau2 = d[k+1];
    f[k] = c * tau1 - s * tau2;
    d[k+1] = s * tau1 + c * tau2;
    if (k < n-2) {
      tau1 = 0.0; tau2 = f[k+1];
      u = -s * f[k+1];
      f[k+1] = c * f[k+1];
      y = f[k];
      z = u;
    }
  }
  return;
}

/*
 * Golub & van Loan Algorithm 8.6.2
 * SVD of a bidiagonal matrix
 */
void
bdsvd(Real d[], Real f[], int n)
{
  int i, p, q, k;
  Real y, z, c, s, epsB;
  Real eps = __DBL_EPSILON__;		// too small ??

  epsB = fabs(d[n-1]);
  for (i = 0; i < n-1; i++) {
    z = fabs(d[i]) + fabs(f[i]);
    if (z > epsB) epsB = z;             // consider more preferable norm
  }
  epsB *= eps;

  q = 0;
  while (q < n) {

    for (i = 0; i < n-1; i++)
      if (fabs(f[i]) <= eps * (fabs(d[i]) + fabs(d[i+1])))
        f[i] = 0.0;

    // q is the size of the largest diagonal matrix, B33, at lower right of B.
    // diagonals may contain zeros.
    for (i = n-2; i >= 0 && (f[i] == 0.0); i--)
      ;
    if (i < 0) {
      q = n; p = 0;
    } else {
      q = n - i - 2;
      for (; i >= 0 && (f[i] != 0.0); i--)
        ;
      // (n-p-q) is the size of a unreduced bidiagonal adjacent to B33.
      if (i < 0)
        p = 0;
      else
        p = i + 1;
    }

    if (q < n) {

      for (i = p; i < n-q; i++)
        if (fabs(d[i]) <= epsB)
          break;

      if (i < n-q) {
        d[i] = 0.0;
        if (i < n-q-1) {
          // zero out the i-th row
          y = f[i]; if (i < n-1) f[i] = 0;
          for (k = i+1; k < n-q; k++) {
            z = d[k];
            givens2(y, z, &c, &s);
            grs_push(i-p, k-p, c, s, p, n-p-q, GR_ROW);
            d[k] = s * y + c * z;
	    if (k < n-1) {
	      y = -s * f[k];
	      f[k] = c * f[k];
	    }
          }
        }
        // zero out the i-th column
        z = f[i-1]; f[i-1] = 0;
        for (k = i-1; k >= p; k--) {
          y = d[k];
          givens(y, z, &c, &s);
          grs_push(k-p, i-p, c, s, p, n-p-q, GR_COL);
          d[k] = c * y - s * z;
          if (k > 0) {
            z = s * f[k-1];
            f[k-1] = c * f[k-1];
          }
        }
      } else
        golubkahn(d + p,  f + p,  p,  n - p - q);
    }
  }
  return;
}


/*
 * Implicit symmetric QR step
 * Golub & van Loan Algorithm 8.3.2
 */
static void
isQRstep(Real d[], Real f[], int n, int p) // isQRstep(Real d[], Real f[], int n)
{
  int k;
  Real a, mu, x, z, c, s, dk, dk1, fk;
  int offset;

  offset = p; // offset = 0;
  a = (d[n-2] - d[n-1]) / 2.0;
  mu = d[n-1] - f[n-2]*f[n-2] /
    ( a + copysign(1.0, a) * sqrt( a * a + f[n-2]*f[n-2] ));
  x = d[0] - mu;
  z = f[0];
  for (k = 0; k < n-1; k++) {
    givens(x, z, &c, &s);
    grs_push(k, k+1, c, s, offset, n, GR_COL);
    if (k > 0)
      f[k-1] = c * x - s * z;
    dk = d[k]; dk1 = d[k+1]; fk = f[k];
    d[k] = c*c*dk + s*s*dk1 - 2.0*c*s*fk;
    f[k] = c*s*(dk - dk1) + (c*c-s*s)*fk;
    d[k+1] = s*s*dk + c*c*dk1 + 2.0*c*s*fk;
    if (k < n-2) {
      x = f[k];
      z = -s * f[k+1];
      f[k+1] = c*f[k+1];
    }
  }
  return;
}

/*
 * Eigen decomposition of a symmetric tridiagonal matrix. 
 * Golub & van Loan Algorithm 8.3.3
 */
static void
tdeig(Real d[], Real f[], int n)
{
  int i, p, q;
  Real eps = __DBL_EPSILON__;

  q = 0;
  while (q < n) {

    for (i = 0; i < n-1; i++)
      if (fabs(f[i]) <= eps * (fabs(d[i]) + fabs(d[i+1])))
        f[i] = 0.0;

    // q is the size of the largest diagonal matrix, T33, at lower right of T.
    // diagonals may contain zeros.
    for (i = n-2; i >= 0 && (f[i] == 0.0); i--)
      ;
    if (i < 0) {
      q = n; p = 0;
    } else {
      q = n - i - 2;
      for (; i >= 0 && (f[i] != 0.0); i--)
        ;
      // (n-p-q) is the size of a unreduced tridiagonal adjacent to T33.
      if (i < 0)
        p = 0;
      else
        p = i + 1;
    }
    if (q < n)
      isQRstep(d + p,  f + p,  n - p - q, p); // isQRstep(d + p,  f + p,  n - p - q);
  }
  return;
}



/*
 * sort eigen(singular)values
 * and rearrange columns of the associated orthogonal matrices.
 * dims: d(n), U(m,m), V(n,n), w(2n)
 * w is an work area of 2n.
 */
static void
rearrange(Real d[], Real *U, Real *V, int m, int n, int *w)
{
  Real s;
  int i, j, k;
  int *ind, *indto;

  ind = w;
  indto = ind + n;

  for (i = 0; i < n; i++)
    ind[i] = i;
  // sort
  for (i = 0; i < n-1; i++) {
    k = i;
    s = d[ind[k]];		// s = fabs(d[ind[k]]);
    for (j = i+1; j < n; j++) {
      if (d[ind[j]] > s) {		// if (fabs(d[ind[j]]) > s) {
        s = d[ind[j]];			// s = fabs(d[ind[j]]);
        k = j;
      }
    }
    if (k > i) {
      j = ind[i];
      ind[i] = ind[k];
      ind[k] = j;
    }
  }
  // replace
  for (i = 0; i < n; i++)
    indto[ind[i]] = i;
  for (k = 0; k < n-1; k++) {
    if (ind[k] > k) {
      s = d[k];
      d[k] = d[ind[k]];
      d[ind[k]] = s;
      for (i = 0; i < m; i++) {
	s = U[i*m+k];
	U[i*m+k] = U[i*m+ind[k]];
	U[i*m+ind[k]] = s;
      }
      for (i = 0; i < n; i++) {
        s = V[i*n+k];
        V[i*n+k] = V[i*n+ind[k]];
        V[i*n+ind[k]] = s;
      }
      //
      indto[ind[k]] = indto[k];
      ind[indto[k]] = ind[k];
    }
  }
  return;
}


/*
 * SVD
 * singular values are returned in w[0..n-1]
 * w is also used as a work area of 2m+2n.
 */
int
svdgl(Real *A, Real *U, Real *V, int m, int n, Real *w)
{
  Real *B, *bd, *bf, *betR, *betC, *v;
  int i;

  bd = w;
  betR = bd + n;
  betC = betR + m;
  v = betC + n;			// v[0..m-1]
  bf = betC + n;

  B = A;
  bidiag(B, m, n, betR, v, 0);
  for (i = 0; i < n - 1; i++) {
    bd[i] = B[i*n+i];
    bf[i] = B[i*n+i+1];
  }
  bd[n-1] = B[n*n-1];
  grs_init();
  bdsvd(bd, bf, n);
  if (U) identityvectmatrix(U, m);
  if (V) identityvectmatrix(V, n);
  givensMatricesUV(U, V, m, n);
  if (U) houseMatrixU(A, betR, m, n, n, 0, U, v);
  if (V) houseMatrixV(A, betC, n, n-2, 0, V, v);
  rearrange(bd, U, V, (U != NULL ? m : 0), (V != NULL ? n : 0), (int *)(w + n));
  
  return 0;
}


/*
 * EVD of a symmetric matrix A = V Diag(w) V^T.
 * eigenvalues are returned in w[0..n-1]
 * w is also used as a work area of 4n
 */
int
seiggl(Real *A, Real *V, int n, Real *w)
{
  Real *T, *td, *tf, *beta, *v;
  int i, *ind;

  td = w;
  tf = td + n;
  beta = tf + n;
  v = beta + n;
  ind = (int *)tf;

  T = A;
  tridiag(T, n, beta, w, 0);		// w[0..3n-1]
  for (i = 0; i < n - 1; i++) {
    td[i] = T[i*n+i];
    tf[i] = T[i*n+i+1];
  }
  td[n-1] = T[n*n-1];
  grs_init();
  tdeig(td, tf, n);
  if (V) identityvectmatrix(V, n);
  givensMatricesUV(0, V, 0, n);
  if (V) houseMatrixV(A, beta, n, n-2, 0, V, v);

  rearrange(td, NULL, V, 0, (V != NULL ? n : 0), ind);
  return 0;
}
