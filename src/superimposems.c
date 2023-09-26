/*
 * superimpose;
 *   complete-case (no missing atoms nor disorder) and
 *   diagonal isotropic (working) covariance Sigma.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "atypes.h"
#include "matrix.h"
#include "memory.h"
#include "msfit.h"
#include "lapack.h"

extern int UseMKL;

/*
 * translate M so that M^T S^{-1} 1 = 0 if diagS != 0, otherwise to M^T 1 = 0,
 * where M is the sample average of {Yj}, i.e., (sum Yj) / N.
 * There must be no missing atom on M.
 */
int
centerM(Real *mp, Real *diagS, int natoms, Real *t)
{
  Real d, s, shift[3];
  int iatom;

  shift[0] = shift[1] = shift[2] = d = 0.0;
  if (diagS != 0)
    for (iatom = 0; iatom < natoms; iatom++) {
      s = 1.0 / diagS[iatom];
      shift[0] += mp[iatom*3+0] * s;
      shift[1] += mp[iatom*3+1] * s;
      shift[2] += mp[iatom*3+2] * s;
      d += s;		/* 1^T Sigma^{-1} 1 */
    }
  else
    for (iatom = 0; iatom < natoms; iatom++) {
      shift[0] += mp[iatom*3+0];
      shift[1] += mp[iatom*3+1];
      shift[2] += mp[iatom*3+2];
      d += 1.0;
    }
  shift[0] /= d;
  shift[1] /= d;
  shift[2] /= d;
  /* translate */
  for (iatom = 0; iatom < natoms; iatom++) {
    mp[iatom*3+0] -= shift[0];
    mp[iatom*3+1] -= shift[1];
    mp[iatom*3+2] -= shift[2];
  }
  if (t) (void) memcpy(t, shift, sizeof(Real) * 3);
  return 0;
}


/*
 * complete-case version
 * tj is calculated as the center of Xj, assuming M has been centerized.
 * mp : the pointer to Mbuf
 */
int
estimate_t(Real *xp, int *op, Real *diagS, int natoms, Real *t)
{
  Real d;
  Real s, xs1[3];
  int i;

  d = 0.0;

  /* xs1 = 1^T S^{-1} Xj */
  xs1[0] = xs1[1] = xs1[2] = 0.0;
  for (i = 0; i < natoms; i++) {
    if (op == NULL || op[i] == 0)
      continue;
    s = 1.0 / diagS[i];
    xs1[0] += s * xp[i*3+0];
    xs1[1] += s * xp[i*3+1];
    xs1[2] += s * xp[i*3+2];
    d += s;		/* = 1^T S^{-1} 1 */
  }

  for (i = 0; i < 3; i++)
    t[i] = xs1[i] / d;

  return 0;
}


#ifdef DEBUG
#define WKSIZ_DGEEV3x3 200
/* dgeev_("V", "V", ...) reported opt worksiz = 102, on a 3x3 matrix */

/*
 * prints eigenvalues of a 3x3 matrix vp
 * also prints the det which is calculated assuming as a rotation matrix
 */
/*
void checkVT(Real *vp)
{
  int n, i, j, k, lwork, info;
  Real A[9], s, wr[3], wi[3], VL[9], VR[9], work[WKSIZ_DGEEV3x3];

  n = 3;
  for (i = 0; i < n * n; i++)
    A[i] = vp[i];
  lwork=WKSIZ_DGEEV3x3;
  dgeev_("V", "V", &n, A, &n, wr, wi, VL, &n, VR, &n, work, &lwork, &info);
  printf("vp eig(%d) = %e %e %e\n", info, wr[0], wr[1], wr[2]);
  printf("            %e %e %e\n", wi[0], wi[1], wi[2]);
  for (k = 0; k < n; k++)
    if (fabs(fabs(wr[k]) - 1.) < 1.e-10) break;
  s = wr[k];
  i = (k+1) % 3;
  j = (i+1) % 3;
  s *= wr[i] * wr[j] - wi[i] * wi[j];
  printf("  det info  %e %d\n", s, info);
}
*/
#undef WKSIZ_DGEEV3x3
#endif


/*
 * aligns M along its pricipal axes.
 */
int
alignM(Real *mp, int natoms, Real *vout)
{
  Real *xp, vp[9], s, *work, lambda[3], cn[3];
  int i, j, k, ndim, lwork, info;
  Real *matM;
  Real up[1];
  Matrix3 V;

  ndim = 3;

  matM = get_workarea("matM", natoms*3, sizeof(Real));
  if (matM == 0)
    abnormalexit("alignM() cannot creat an arrray matM", 0);
  if (UseMKL) {
    /* store M on matM in the column-major order */
    for (i = 0; i < natoms; i++)
      for (j = 0; j < ndim; j++)
	matM[i+j*natoms] = mp[i*ndim+j]; // / sqrt(natoms);
    lwork = -1;
    dgesvd_("N", "A", &natoms, &ndim, matM, &natoms, lambda, up, &natoms,
	    vp, &ndim, &s, &lwork, &info);
    lwork = (int) s;
    work = get_workarea("dgesvd_work", lwork, sizeof(Real));
    if (work == 0)
      abnormalexit("alignM() cannot creat an arrray dgesvd_work", 0);
    /* Compute SVD: M = U lam V^T  in the column-major order. */
    dgesvd_("N", "A", &natoms, &ndim, matM, &natoms, lambda, up, &natoms,
	    vp, &ndim, work, &lwork, &info);
    remove_workarea("dgesvd_work");
  } else {
    work = get_workarea("svdgl_work", natoms*2+ndim*2, sizeof(Real));
    if (work == 0)
      abnormalexit("alignM() cannot creat an arrray svdgl_work", 0);
    for (i = 0; i < natoms * ndim; i++)
      matM[i] = mp[i];
    svdgl(matM, NULL, vp, natoms, ndim, work);
    lambda[0] = work[0]; lambda[1] = work[1]; lambda[2] = work[2];
    remove_workarea("svdgl_work");
  }
  V = (Matrix3) vp;	/* V^T in the col-major order is V in the row-major */
  if (det3x3(V) < 0) {
    i = 2;
    if (lambda[1] < lambda[i]) i = 1;
    if (lambda[0] < lambda[i]) i = 0;
    V[0][i] *= -1.0;   V[1][i] *= -1.0;   V[2][i] *= -1.0;
  }
  cn[0] = cn[1] = cn[2] = 0.0;
  for (i = natoms-1; fabs(cn[0] * cn[1]) < 1e-8 && i > 0; i--) {
    cn[0] = mp[0*3+0] - mp[i*3+0];
    cn[1] = mp[0*3+1] - mp[i*3+1];
    cn[2] = mp[0*3+2] - mp[i*3+2];
  }
  s = cn[0]*V[0][0] + cn[1]*V[1][0] + cn[2]*V[2][0];
  if (s < 0.0)
    for (i = 0; i < 3; i++)
      for (j = 0; j < 2; j++)
	V[i][j] *= -1;
  s = cn[0]*V[0][2] + cn[1]*V[1][2] + cn[2]*V[2][2];
  if (s < 0.0)
    for (i = 0; i < 3; i++)
      for (j = 1; j < 3; j++)
  	V[i][j] *= -1;
  /*
  printf("V in alignM()\n");
  for (i = 0; i < 3; i++) 
    printf("\t%e %e %e\n", V[i][0], V[i][1], V[i][2]);
  */
  xp = matM;
  for (i = 0; i < natoms; i++)
    for (j = 0; j < 3; j++) {
      s = 0.0;
      for (k = 0; k < 3; k++)
	s += mp[i*3+k] * vp[k*3+j];
      xp[i*3+j] = s;
    }
  for (i = 0; i < natoms * 3; i++)
    mp[i] = xp[i];
  if (vout)
    for (i = 0; i < 9; i++)
      vout[i] = vp[i];

  /*
#ifdef DEBUG
  checkVT(vp);
#endif
  */
  remove_workarea("matM");
  return 0;
}

#define WKSIZ_DGESVD3x3 30		/* 15 would suffice */
					/* w of svdgl() >= 3 x 4 */
/*
 * R = U V^T s.t. A = U S V^T
 * vpout (if not 0) : V^T in the row-major order
 */
void
rotation_matrix_svd(Real *ap, Real *rp, Real *vpout)
{
  Real up[9], vp[9], lambda[3], wk[WKSIZ_DGESVD3x3], s;
  int i, j, k, n, lwork, info;
  Matrix3 R, U, VT;

  R = (Matrix3) rp;
  U = (Matrix3) up;
  VT = (Matrix3) vp;

  n = 3;
  if (UseMKL) {
    lwork = -1;
    dgesvd_("A","A", &n, &n, ap, &n, lambda, vp, &n, up, &n, &s, &lwork, &info);
    lwork = (int) s;
    if (lwork > WKSIZ_DGESVD3x3)
      abnormalexit("rotation_matrix_svd(): increase wksiz", 0);
    /* Compute SVD: A^T = V Sig U^T in the column-major order */  
    dgesvd_("A","A", &n, &n, ap, &n, lambda, vp, &n, up, &n, wk, &lwork, &info);
    if (info != 0)
      fprintf(stderr, "rotation_matrix_svd(): dgesvd returned %d\n", info);
  } else {
    svdgl(ap, up, vp, n, n, wk);
    s = vp[3]; vp[3] = vp[1]; vp[1] = s;
    s = vp[6]; vp[6] = vp[2]; vp[2] = s;
    s = vp[7]; vp[7] = vp[5]; vp[5] = s;
    lambda[0] = wk[0]; lambda[1] = wk[1]; lambda[2] = wk[2];
  }
  /* up and vp are U and V^T in the row-major order, respectively. */
  if (det3x3(U) * det3x3(VT) < 0) {
    i = 2;
    if (lambda[1] < lambda[i]) i = 1;
    if (lambda[0] < lambda[i]) i = 0;
    U[0][i] = -U[0][i];
    U[1][i] = -U[1][i];
    U[2][i] = -U[2][i];
  }
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      s = 0.0;
      for (k = 0; k < 3; k++)
        s += U[i][k] * VT[k][j];
      R[i][j] = s;
    }

  if (vpout)
    for (i = 0; i < 3*3; i++)
      vpout[i] = vp[i];

  return;
}

/*
 * complete-case version
 * estimates R for the conformer Xj(xp).
 * R is the covar-weighted cross-covariance between the Xj and M.
 * mp : the pointer to Mbuf
 * work: natoms * 7
 */
int
estimate_R(Real *xp, Real *t, int *op, Real *mp, Real *diagS, int natoms,
	   Real *rp)
{
  Real ap[9], s;
  Matrix3 A, X;
  int iatom, i, j;

  A = (Matrix3) ap;
  zerovectmatrix(ap, 3, 3);

  X = (Matrix3) xp;

  /* (Xj - 1 \otimes tj^T)^T S^{-1} M */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      s = 0.0;
      for (iatom = 0; iatom < natoms; iatom++) {
	if (op == NULL || op[iatom] == 0)
	  continue;
	s += (X[iatom][i] - t[i]) * mp[iatom*3+j] / diagS[iatom];
      }
      A[i][j] = s;
  }

  rotation_matrix_svd(ap, rp, 0);

  return 0;
}

/*
 * complete-case version
 * calculate the sample average of {Yj} excluding Yiexcl and stores on mp.
 * To calculate avarage over all {Yj}, call with iexcl == -1.
 * work: sizeof natoms
 */
int
average_Y(Conformer *cp, int nconfms, int natoms, int iexcl,
	  int r1weight, Real *mp, Real *work)
{
  Real s, *wp, *yp;
  Matrix3 M, Y;
  int iconfm, iatom, *op;

  M = (Matrix3) mp;
  wp = work;

  /* M <- 0, w <- 0 */
  zerovectmatrix(mp, natoms, 3);
  zerovectmatrix(wp, natoms, 1);

  /* M <- sum Yj / rmsdj, w <- sum 1/rmsdj */
  for (iconfm = 0; iconfm < nconfms; iconfm++, cp++) {
    if (iconfm == iexcl)
      continue;
    op = cp->o;
    yp = cp->y;
    Y = (Matrix3) yp;
    for (iatom = 0; iatom < natoms; iatom++)  {
      if (op == NULL || op[iatom] == 0)
	continue;
      s = 1.0;
      if (r1weight) s /= cp->rmsd;
      wp[iatom] += s;
      M[iatom][0] += s * Y[iatom][0];
      M[iatom][1] += s * Y[iatom][1];
      M[iatom][2] += s * Y[iatom][2];
    }
  }
    
  /* M <- W^{-1} M. */
  for (iatom = 0; iatom < natoms; iatom++) {
    M[iatom][0] /= wp[iatom];
    M[iatom][1] /= wp[iatom];
    M[iatom][2] /= wp[iatom];
  }

  return 0;
}



void
calculate_Y(Real *vectX, Real *vectR, Real *t, int natoms, Real *vectY)
{
  int i;
  Real sx[3];
  Matrix3 X, R, Y;

  X = (Matrix3) vectX;
  R = (Matrix3) vectR;
  Y = (Matrix3) vectY;
  for (i = 0; i < natoms; i++) {
    sx[0] = X[i][0] - t[0];
    sx[1] = X[i][1] - t[1];
    sx[2] = X[i][2] - t[2];
    Y[i][0] = sx[0]*R[0][0] + sx[1]*R[1][0] + sx[2]*R[2][0];
    Y[i][1] = sx[0]*R[0][1] + sx[1]*R[1][1] + sx[2]*R[2][1];
    Y[i][2] = sx[0]*R[0][2] + sx[1]*R[1][2] + sx[2]*R[2][2];
  }
}
