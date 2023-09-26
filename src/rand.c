#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "atypes.h"
#include "memory.h"
#include "rand.h"
#include "matrix.h"
#include "lapack.h"

extern int UseMKL;

extern void abnormalexit(char *, int);

/* use random() as RNG */

/* other possiblly better RNGs: Mersenne-Twister, xorshift+, ... */

void rand_seed(int seed)
{
  srandom((unsigned int) seed);
}

/* a random number in [0,1) */
Real
rand_uniformc()
{
  return (double) random() / ((double) RAND_MAX + 1.0);
}

/* a random number in (0,1) */
Real
rand_uniformo()
{
  return ((double) random() + 1.0) / ((double) RAND_MAX + 2.0);
}

/* a random number drawn from N(0,1). Box-Muller. */
Real
rand_stdnorm()
{
  static Real z1;
  double u1, u2, r;

  u1 = rand_uniformo();
  u2 = rand_uniformo();
  r = sqrt(-2.0 * log(u1));
  z1 = r * cos(2.0 * M_PI * u2);
  // z2 = r * sin(2.0 * M_PI * u2); 	// alway discard z2
  return z1;
}

/*
 * a random vector xp drawn from the multivariate normal dist N(mp, cov)
 * cov : n x n matrix stored in the row major order of n columns.
 */
void
rand_mvnorm_evd(Real *mp, Real *cov, int n, Real *xp)
{
  Real *ap, *up, *sp, s, *zp, *work;
  int lwork, info, i, j, k;
  
  ap = get_workarea("A", n * n, sizeof(Real));
  up = get_workarea("U", n * n, sizeof(Real));
  if (ap == 0 || up == 0)
    abnormalexit("rand_mvnorm_evd(): cannot alloc A/U\n", 0);
  memcpy(ap, cov, n * n * sizeof(Real));	// Covariance

  if (UseMKL) {
    lwork = -1;
    sp = 0;
    dsyev_("V", "U", &n, ap, &n, sp, &s, &lwork, &info);
    lwork = (int) s;
    sp = get_workarea("S", n + lwork, sizeof(Real));
    work = sp + n;
    if (sp == 0)
      abnormalexit("rand_mvnorm_evd(): cannot alloc workarea\n", 0);
    /* A = A^T = U S U^T in cmo. */
    dsyev_("V", "U", &n, ap, &n, sp, work, &lwork, &info);
    if (info != 0)
      abnormalexit("rand_mvnorm_evd(): failure in eigendecomp\n", 0);
    /* transpose & turn into descending order (of eig values) column wise */
    for (i = 0; i < n; i++)
      for (j = 0; j <= i; j++) {
	up[i*n+n-1-j] = ap[i+j*n];
	up[j*n+n-1-i] = ap[j+i*n];
      }
    for (i = 0; i < n / 2; i++) {
      s = sp[n-1-i];
      sp[n-1-i] = sp[i];
      sp[i] = s;
    }
    /* strictly, sp[] is in non-acending order. If degenerated,
       up[] is not unique in the additional freedoms of rotations. */
  } else {
      sp = get_workarea("S", 4 * n, sizeof(Real));
    if (sp == 0)
      abnormalexit("rand_mvnorm_evd(): cannot alloc workarea\n", 0);
    seiggl(ap, up, n, sp);
  }

  for (k = 0; k < n; k++)
    if (sp[k] < 0.0)
      abnormalexit("rand_mvnorm_evd(): A is not PSD\n", 0);
    else
      sp[k] = sqrt(sp[k]);
  /* A <- U*S^1/2*U^T */ 
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      s = 0.0;
      for (k = 0; k < n; k++)
	s += up[i*n + k] * sp[k] * up[j*n + k];
      ap[i*n + j] = s;
    }

  zp = sp;
  for (i = 0; i < n; i++) 
    zp[i] = rand_stdnorm();
  for (i = 0; i < n; i++) {
    s = 0.0;
    for (k = 0; k < n; k++)
      s += ap[i*n + k] * zp[k];
    xp[i] = s + mp[i];
  }
  remove_workarea("S");
  remove_workarea("U");
  remove_workarea("A");
  return;
}

int
rand_mvnorm(Real *mp, Real *cov, int n, Real *xp, int fact)
{
  Real *ap, s, *zp;
  int i, k;
  
  ap = get_workarea("A", n * n, sizeof(Real));
  zp = get_workarea("Z", n, sizeof(Real));
  if (ap == 0 || zp == 0)
    abnormalexit("rand_mvnorm(): cannot alloc A/z\n", 0);

  memcpy(ap, cov, n * n * sizeof(Real));	// Covariance or its decomp

  if (!fact) {
    i = cholesky_decomp(ap, n);
    if (i > 0) return i;
  }

  for (i = 0; i < n; i++) 
    zp[i] = rand_stdnorm();
  for (i = 0; i < n; i++) {
    s = 0.0;
    for (k = 0; k <= i; k++)
      s += ap[i*n + k] * zp[k];
    xp[i] = s + mp[i];
  }

  remove_workarea("Z");
  remove_workarea("A");
  return 0;
}

/*
 * a random n x 3 matrix X drawn from the matrix normal dist MN(mp, cov, I3) 
 * mp : n x 3 matrix stored in the row major order of 3 columns.
 * cov : n x n matrix stored in the row major order of n columns.
 * xp : (on exit) an nx3 random matrix ~ MN(M,C,I3), in the row major order
 * fact: ==1 cov is actually a factorized matrix obtained using chol_decomp()
 */
int
rand_m3norm(Real *mp, Real *cov, int n, Real *xp, int fact)
{
  Real *rx, *ry, *rz, *mx, *my, *mz;
  int i, j;

  rx = get_workarea("rxyz", n * 3, sizeof(Real));
  ry = rx + n;
  rz = ry + n;
  mx = get_workarea("mxyz", n * 3, sizeof(Real));
  my = mx + n;
  mz = my + n;
  for (i = 0; i < n; i++)
    for (j = 0; j < 3; j++)
      mx[i+j*n] = mp[i*3+j];
  if ((i = rand_mvnorm(mx, cov, n, rx, fact)) > 0)
    return i;
  if ((i = rand_mvnorm(my, cov, n, ry, fact)) > 0)
    return i;
  if ((i = rand_mvnorm(mz, cov, n, rz, fact)) > 0)
    return i;
  for (i = 0; i < n; i++)
    for (j = 0; j < 3; j++)
      xp[i*3+j] = rx[i+j*n];
  remove_workarea("mxyz");
  remove_workarea("rxyz");
  return 0;
}


#define RAND_ROTMAT_EPS (1.0e-8)

/* 
 * generates a 3x3 rotation matrix R via a random quaternion
 * https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
void
orig_rand_rotmat(Real *rp)
{
  Real a, b, c, d, s, phi;
  Matrix3 R;

  /* pick a point in N^3(0,1) */
  while (1) {
    b = rand_stdnorm();
    c = rand_stdnorm();
    d = rand_stdnorm();
    s = b*b + c*c + d*d;
    if (s >= RAND_ROTMAT_EPS)
      break;
  }
  s = sqrt(s);
  b /= s;
  c /= s;
  d /= s;
  /* choose a angle phi in [0,\pi] to form a versor (a,b,c,d). */
  phi = rand_uniformc() * M_PI;
  a = cos(phi);
  s = sin(phi);
  b *= s;
  c *= s;
  d *= s;
  
  R = (Matrix3) rp;
  /* covert the versor to the orthogonal rotation matrix */
  R[0][0] = a*a + b*b - c*c - d*d;
  R[0][1] = 2*b*c - 2*a*d;
  R[0][2] = 2*b*d + 2*a*c; 
  R[1][0] = 2*b*c + 2*a*d;
  R[1][1] = a*a - b*b + c*c - d*d;
  R[1][2] = 2*c*d - 2*a*b;
  R[2][0] = 2*b*d - 2*a*c;
  R[2][1] = 2*c*d + 2*a*b;
  R[2][2] = a*a - b*b - c*c + d*d;

  return;
}

void
rand_rotmat(Real *rp)
{
  long double a, b, c, d, s, phi;
  Matrix3 R;

  /* pick a point in N^3(0,1) */
  while (1) {
    b = rand_stdnorm();
    c = rand_stdnorm();
    d = rand_stdnorm();
    s = b*b + c*c + d*d;
    if (s >= RAND_ROTMAT_EPS)
      break;
  }
  s = sqrtl(s);
  b /= s;
  c /= s;
  d /= s;
  /* choose a angle phi in [0,\pi] to form a versor (a,b,c,d). */
  phi = rand_uniformc() * M_PI;
  a = cosl(phi);
  s = sinl(phi);
  b *= s;
  c *= s;
  d *= s;
  
  R = (Matrix3) rp;
  /* covert the versor to the orthogonal rotation matrix */
  R[0][0] = a*a + b*b - c*c - d*d;
  R[0][1] = 2*b*c - 2*a*d;
  R[0][2] = 2*b*d + 2*a*c; 
  R[1][0] = 2*b*c + 2*a*d;
  R[1][1] = a*a - b*b + c*c - d*d;
  R[1][2] = 2*c*d - 2*a*b;
  R[2][0] = 2*b*d - 2*a*c;
  R[2][1] = 2*c*d + 2*a*b;
  R[2][2] = a*a - b*b - c*c + d*d;

  return;
}


/*

main()
{
  Real C[] = {1.0,-0.0,0.0,-0.4,-0.4,0.1,-0.1,0.0,-0.1,0.1,-0.1,
	      -0.0,1.0,-0.1,0.0,0.1,-0.1,0.0,-0.0,-0.0,-0.0,-0.0,
	      0.0,-0.1,1.0,-0.2,-0.0,-0.1,0.0,0.0,0.0,0.1,0.1,
	      -0.4,0.0,-0.2,1.0,0.2,-0.4,-0.0,-0.2,0.0,-0.4,-0.2,
	      -0.4,0.1,-0.0,0.2,1.0,-0.4,0.1,0.0,0.0,-0.1,0.0,
	      0.1,-0.1,-0.1,-0.4,-0.4,4.0,0.0,-0.1,-0.1,0.4,0.4,
	      -0.1,0.0,0.0,-0.0,0.1,0.0,1.0,-0.0,0.0,0.1,0.0,
	      0.0,-0.0,0.0,-0.2,0.0,-0.1,-0.0,4.0,-1.0,-0.2,0.9,
	      -0.1,-0.0,0.0,0.0,0.0,-0.1,0.0,-1.0,1.0,-0.0,0.0,
	      0.1,-0.0,0.1,-0.4,-0.1,0.4,0.1,-0.2,-0.0,1.0,0.0,
	      -0.1,-0.0,0.1,-0.2,0.0,0.4,0.0,0.9,0.0,0.0,1.0
};
  Real M[] = {0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0,
0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0,
0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0
  };
  Real x[11*3];
  int i, j, k, n = 11;

  init_workarea(1000000);
  for (i = 0; i < 100; i++)
    (void)rand_stdnorm();

  for (i = 0; i < n*3; i++)
    M[i] = 0;


  rand_rotmat(x);
  checkVT(x);
  for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
    M[i+j*3]=x[i*3+j];
  checkVT(M);
  exit(0);


  for (i = 0; i < 30000; i++) {
    rand_m3norm(M, C, n, x);
    // rand_mvnorm(M, C, n, n, x);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < n-1; k++)
	printf("%g,", x[j+k*3]);
      //printf("%g, %g, %g\n", x[k*3],x[k*3+1],x[k*3+2]);
      printf("%g\n", x[j+k*3]);
    }
  }
}
    
*/
