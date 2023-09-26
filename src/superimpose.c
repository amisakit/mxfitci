#include <math.h>
#include <string.h>
#include <stdio.h>
#include "atypes.h"
#include "lapack.h"
#include "matrix.h"
#include "memory.h"
#include "molec.h"

extern int UseMKL;
extern void abnormalexit(char *message, int code);

/*
 * translate M so that M^T S^{-1} 1 = 0 if diagS != 0, otherwise to M^T 1 = 0.
 * Notice: M must be complete, i.e., o[i] == 1 for all 0 <= i < natoms
 * work: natoms * 7
 */
int
center_M(Real *mp, Real diagS[], int natoms, Real *t)
{
  Real d, shift[3], s;
  int i;

  // inv_covar_weighted_matrix(mp, covp, natoms, dodiag, smp, work);
  shift[0] = shift[1] = shift[2] = d = 0.0;
  for (i = 0; i < natoms; i++) {
    s = 1.0 / diagS[i];
    shift[0] += mp[i*3+0] * s;
    shift[1] += mp[i*3+1] * s;
    shift[2] += mp[i*3+2] * s;
    d += s;           /* 1^T Sigma^{-1} 1 */
  }
  shift[0] /= d;
  shift[1] /= d;
  shift[2] /= d;
  /* translate */
  for (i = 0; i < natoms; i++) {
    mp[i*3+0] -= shift[0];
    mp[i*3+1] -= shift[1];
    mp[i*3+2] -= shift[2];
  }
  if (t) (void) memcpy(t, shift, sizeof(Real) * 3);
  return 0;
}


/*
 * translate M so that sum(M-miZi)^T S^{-1} 1 = 0. 
 * Notice: M must be complete, i.e., o[i] == 1 for all 0 <= i < natoms
 * work: natoms * 7
 */
int
center_MZ(Real *mp, Real Zbuf[], Real diagS[], int nconfms[], int ngroups, int nconfms_tot, int natoms, Real *t, Real *work)
{
  Real d, s, *mz, *zp;
  int i, iens;
  
  // MZ <- M + (1/m) sum_i mi Zi
  mz = work;
  for (i = 0; i < natoms*3; i++)
    mz[i] = mp[i];
  zp = Zbuf;
  for (iens = 0; iens < ngroups; iens++, zp += natoms*3) {
    s = nconfms[iens] / nconfms_tot;
    for (i = 0; i < natoms; i++) {
      mz[i*3+0] += zp[i*3+0] * s;
      mz[i*3+1] += zp[i*3+1] * s;
      mz[i*3+2] += zp[i*3+2] * s;
    }
  }
  
  // t <- (M + (1/m) sum_i mi Zi)^T Sigma^{-1} 1
  t[0] = t[1] = t[2] = 0.0;
  for (i = 0; i < natoms; i++) {
    s = 1.0 / diagS[i];
    d += s;           /* 1^T Sigma^{-1} 1 */
    t[0] += mz[i*3+0] * s;
    t[1] += mz[i*3+1] * s;
    t[2] += mz[i*3+2] * s;
  }
  // t <- 1/(1 Sigma^{-1} 1) * t
  t[0] /= d;
  t[1] /= d;
  t[2] /= d;
  /* translate */
  for (i = 0; i < natoms; i++) {
    mp[i*3+0] -= t[0];
    mp[i*3+1] -= t[1];
    mp[i*3+2] -= t[2];
  }
  return 0;
}


void
estimate_M_EM(Real *Ybuf, Real *Zbuf, int nconfms[], int ngroups, int natoms,
	      Real *mp)
{
  Real s, *yp, *zp;
  Matrix3 M, Y, Z;
  int iconfm, i, nconfms_tot, iens;

  M = (Matrix3) mp;

  /* M <- 0 */
  zerovectmatrix(mp, natoms, 3);

  /* M <- sum_ij Yij - sum_i mi Zi */
  nconfms_tot = 0;
  yp = Ybuf;
  zp = Zbuf;
  for (iens = 0; iens < ngroups; iens++) {
    for (iconfm = 0; iconfm < nconfms[iens]; iconfm++) {
      Y = (Matrix3) (yp);
      for (i = 0; i < natoms; i++) {
        M[i][0] += Y[i][0];
        M[i][1] += Y[i][1];
        M[i][2] += Y[i][2];
      }
      yp += natoms * 3;
    }
    Z = (Matrix3) zp;
    for (i = 0; i < natoms; i++) {
      s = nconfms[iens];
      M[i][0] -= s * Z[i][0];
      M[i][1] -= s * Z[i][1];
      M[i][2] -= s * Z[i][2];
    }
    zp += natoms * 3;
    nconfms_tot += nconfms[iens];
  }
  /* M <- 1/m (sum_ij Yij - sum_i mi Zi) */
  s = 1.0 / (Real)(nconfms_tot);
  for (i = 0; i < natoms * 3; i++)
    mp[i] *= s;
  return;
}

/*
 * diagV: contents are not used. NULL/non-NULL is relavant.
 * wokr: sizeof natoms * 3 * 2
 */
int
estimate_M_ML(Real *Ybuf, Ensemble *Bbuf, int ngroups, int natoms,
	      int njvars, Real *diagV, Real *mp, Real *work)
{
  Real *yp, *dp, *svi, *hp, *wp, s;
  int iconfm, i, j, k, iens, m_tot;
  const int ndim = 3;

  dp = work;
  wp = work + natoms * ndim;

  if (njvars == 1) {			// balanced design
    /* M <-  sum_i sum_j Yij */
    zerovectmatrix(mp, natoms, ndim);
    m_tot = 0;
    yp = Ybuf;
    for (iens = 0; iens < ngroups; iens++) {
      for (iconfm = 0; iconfm < Bbuf[iens].nconfms; iconfm++) {
	for (i = 0; i < natoms; i++)
	  for (j = 0; j < ndim; j++)
	    mp[i*ndim+j] += yp[i*ndim+j];
	yp += natoms * 3;
      }
      m_tot += Bbuf[iens].nconfms;
    }
    /* M <- M / m */
    for (i = 0; i < natoms; i++)
      for (j = 0; j < ndim; j++)
	mp[i*ndim+j] /= m_tot;

  } else if (diagV) {

    /* [M]ij = [h]ii sum_k [svi]ii [Dk]ij, Dk = sum_j Yij */
    zerovectmatrix(mp, natoms, ndim);
    yp = Ybuf;
    for (iens = 0; iens < ngroups; iens++) {
      svi = Bbuf[iens].jvp->svi;
      zerovectmatrix(dp, natoms, ndim);
      for (iconfm = 0; iconfm < Bbuf[iens].nconfms; iconfm++) {
	for (i = 0; i < natoms; i++)
	  for (j = 0; j < ndim; j++)
	    dp[i*ndim+j] += yp[i*ndim+j];
	yp += natoms * 3;
      }
      for (i = 0; i < natoms; i++)
	for (j = 0; j < ndim; j++)
	  mp[i*ndim+j] += svi[i*natoms+i] * dp[i*ndim+j];
    }
    hp = Bbuf[0].jvp->h;
    for (i = 0; i < natoms; i++)
      for (j = 0; j < ndim; j++)
	mp[i*ndim+j] *= hp[i*natoms+i];

  } else {		// unbalanced & full V

    zerovectmatrix(wp, natoms, ndim);
    yp = Ybuf;
    for (iens = 0; iens < ngroups; iens++) {
      /* dp <- sum_j Yij */
      zerovectmatrix(dp, natoms, ndim);
      for (iconfm = 0; iconfm < Bbuf[iens].nconfms; iconfm++) {
	for (i = 0; i < natoms; i++)
	  for (j = 0; j < ndim; j++)
	    dp[i*ndim+j] += yp[i*ndim+j];
	yp += natoms * 3;
      }
      /* wp += (S + mi V)^{-1} dp */
      svi = Bbuf[iens].jvp->svi;
      for (i = 0; i < natoms; i++)
	for (j = 0; j < ndim; j++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += svi[i*natoms+k] * dp[k*ndim+j];
	  wp[i*ndim+j] += s;
	}
    }
    /* mp <- (sum_i mi (S + mi V)^{-1})^{-1} wp */
    hp = Bbuf[0].jvp->h;
    for (i = 0; i < natoms; i++)
      for (j = 0; j < ndim; j++) {
	s = 0.0;
	for (k = 0; k < natoms; k++)
	  s += hp[i*natoms+k] * wp[k*ndim+j];
	mp[i*ndim+j] = s;
      }
  }

  return 0;
}


/*
 * aligns M along its pricipal axes.
 */
int
align_M(Real *mp, int natoms, Real *vout)
{
  Real *xp, vp[9], s, *work, lambda[3], cn[3];
  int i, j, k, ndim, lwork, info;
  Real *matM;
  Real up[1];
  Matrix3 V;

  ndim = 3;

  matM = get_workarea("matM", natoms*3, sizeof(Real));
  if (matM == 0)
    abnormalexit("align_M() cannot creat an arrray matM", 0);
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
      abnormalexit("align_M() cannot creat an arrray dgesvd_work", 0);
    /* Compute SVD: M = U lam V^T  in the column-major order. */
    dgesvd_("N", "A", &natoms, &ndim, matM, &natoms, lambda, up, &natoms,
            vp, &ndim, work, &lwork, &info);
    remove_workarea("dgesvd_work");
  } else {
    work = get_workarea("svdgl_work", natoms*2+ndim*2, sizeof(Real));
    if (work == 0)
      abnormalexit("align_M() cannot creat an arrray svdgl_work", 0);
    for (i = 0; i < natoms * ndim; i++)
      matM[i] = mp[i];
    svdgl(matM, NULL, vp, natoms, ndim, work);
    lambda[0] = work[0]; lambda[1] = work[1]; lambda[2] = work[2];
    remove_workarea("svdgl_work");
  }
  V = (Matrix3) vp;     /* V^T in the col-major order is V in the row-major */
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

  remove_workarea("matM");
  return 0;
}


#define WKSIZ_DGESVD3x3 30              /* 15 would suffice */
                                        /* w of svdgl() >= 3 x 4 */
/*
 * R = U V^T s.t. A = U S V^T
 * vpout (if not 0) : V^T in the row-major order
 */
static void
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
 * work: natoms x 3
 */
void
estimate_Rt_EM(Ensemble Bbuf[], Real *mp, Real diagS[], int ngroups, int natoms, int m_tot, int Mcentered, int YinM, Real *work)
{
  Matrix3 A, M, R, W, X, Y, Z;
  Real ap[9], s, sumS1, g[3], *t, *rp;
  int i, j, k, iens, iconf, mi;
  Ensemble *ensp;
  Conformer *confp;

  A = (Matrix3) ap;
  M = (Matrix3) mp;
  W = (Matrix3) work;

  /* sumS1 <- 1^T S^{-1} 1 */
  s = 0.0;
  for (i = 0; i < natoms; i++)
    s += 1.0 / diagS[i];
  sumS1 = s;
  
  ensp = Bbuf;
  for (iens = 0; iens < ngroups; iens++, ensp++) {
    Z = (Matrix3) ensp->z;
    mi = ensp->nconfms;
    confp = ensp->confmbp;
    for (iconf = 0; iconf < mi; iconf++, confp++) {
      t = confp->t;
      R = (Matrix3) (rp = confp->r);
      X = (Matrix3) confp->x;
      Y = (Matrix3) confp->y;

      /* g <- (1/sumS1) * X^T S^{-1} 1 */
      for (i = 0; i < 3; i++) {
	s = 0.0;
	for (k = 0; k < natoms; k++)
	  s += X[k][i] / diagS[k];
	g[i] = s / sumS1;
      }

      /* W <- S^{-1} (X - 1 \otimes g) */
      for (i = 0; i < natoms; i++) {
	W[i][0] = (X[i][0] - g[0]) / diagS[i];
	W[i][1] = (X[i][1] - g[1]) / diagS[i];
	W[i][2] = (X[i][2] - g[2]) / diagS[i];
      }
      if (YinM) {
      /* A <- W^T (M + Z - Y/m) */
	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++) {
	    s = 0.0;
	    for (k = 0; k < natoms; k++)
	      s += W[k][i] * (M[k][j] + Z[k][j] - Y[k][j]/m_tot);
	    A[i][j] = s;
	  }
      } else {		// if (!YinM)
	/* A <- W^T (M + Z) */
	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++) {
	    s = 0.0;
	    for (k = 0; k < natoms; k++)
	      s += W[k][i] * (M[k][j] + Z[k][j]);
	    A[i][j] = s;
	  }
      }

      rotation_matrix_svd(ap, rp, 0);

      if (Mcentered && YinM) {
	/* t <- R (Z - Y/m)^T S^{-1} 1 * (sumS1 * (1-1/m)) */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += ( R[i][0] * (Z[k][0] - Y[k][0]/m_tot)
		   + R[i][1] * (Z[k][1] - Y[k][1]/m_tot)
		   + R[i][2] * (Z[k][2] - Y[k][2]/m_tot) ) / diagS[k];
	  t[i] = s / (sumS1 * (1.0 - 1.0/m_tot));
	}
      } else if (!Mcentered && YinM) {
	/* t <- - R (M + Z - Y/m)^T S^{-1} 1 * (sumS1 * (1-1/m)) */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += ( R[i][0] * (M[k][0] + Z[k][0] - Y[k][0]/m_tot)
		   + R[i][1] * (M[k][0] + Z[k][1] - Y[k][1]/m_tot)
		   + R[i][2] * (M[k][0] + Z[k][2] - Y[k][2]/m_tot) ) / diagS[k];
	  t[i] = s / (sumS1 * (1.0 - 1.0/m_tot));
	}
      } else if (Mcentered && !YinM) {
	/* t <- R Z^T S^{-1} 1 * (sumS1 * (1-1/m)) */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += ( R[i][0] * (Z[k][0])
		   + R[i][1] * (Z[k][1])
		   + R[i][2] * (Z[k][2]) ) / diagS[k];
	  t[i] = s / sumS1;
	}
      } else {		// if (!Mcentered && !YinM) {
	/* t <- R (M + Z)^T S^{-1} 1 * sumS1 */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++) {
	    s += ( R[i][0] * (M[k][0] + Z[k][0])
		   + R[i][1] * (M[k][1] + Z[k][1])
		   + R[i][2] * (M[k][2] + Z[k][2]) ) / diagS[k];
	  }
	  t[i] = s / sumS1;
	}
      }

      /* t <- g - t */
      t[0] = g[0] - t[0];
      t[1] = g[1] - t[1];
      t[2] = g[2] - t[2];
    }
  }

  return;
}



void
estimate_Rt_ML(Ensemble Bbuf[], Real *mp, Real *diagV, Real diagS[], int ngroups, int natoms, int Mcentered, Real *work)
{
  int i, j, k, iens, iconf, mi, fullV;
  Real *ap, *cp, *dp, *sc, *t, *wp, *yp, *zp, s, sumSC, g[3];
  Real *MTS, *ZTC;
  Matrix3 A, D, M, W, R, X, Y, Z;
  Ensemble *ensp;
  Conformer *confp;

  fullV = (diagV == NULL);

  ap = work;
  dp = ap + 9;		// natoms x 3。これは main で管理して、loglikと共有するべき。
  wp = dp + natoms * 3;
  zp = wp + natoms * 3;
  MTS = zp + natoms * 3;
  ZTC = MTS + natoms * 3;

  A = (Matrix3) ap;
  D = (Matrix3) dp;
  M = (Matrix3) mp;
  W = (Matrix3) wp;
  Z = (Matrix3) zp;

  /* MTS == M^T S^{-1} */
  for (i = 0; i < natoms; i++) {
    s = 1.0 / diagS[i];
    MTS[0*natoms+i] = M[i][0] * s;
    MTS[1*natoms+i] = M[i][1] * s;
    MTS[2*natoms+i] = M[i][2] * s;
  }

  ensp = Bbuf;
  for (iens = 0; iens < ngroups; iens++, ensp++) {

    cp = ensp->jvp->c;			// Ci
    sc = ensp->jvp->sc;			// S^{-1} - Ci
    mi = ensp->nconfms;

    /* sumSC = 1^T SC 1 == 1^T (S^{-1} - C) 1 */
    if (fullV) {
      s = 0.0;
      for (i = 0; i < natoms*natoms; i++)
	s += sc[i];
      sumSC = s;
    } else {
      s = 0.0;
      for (i = 0; i < natoms; i++)
	s += sc[i*natoms+i];
      sumSC = s;
    }

    // D == sum_j (Yij - M)
    for (i = 0; i < natoms*3; i++)
      dp[i] = -mi * mp[i];
    yp = ensp->confmbp->y;
    for (iconf = 0; iconf < mi; iconf++, yp += natoms*3)
      for (i = 0; i < natoms*3; i++)
	dp[i] += yp[i];

    confp = ensp->confmbp;
    for (iconf = 0; iconf < mi; iconf++, confp++) {
      t = confp->t;
      R = (Matrix3) confp->r;
      X = (Matrix3) confp->x;
      Y = (Matrix3) confp->y;

      /* g <- (1/sumSC) * X^T (S^{-1}-C) 1 */
      if (fullV) {
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (j = 0; j < natoms; j++)
	    for (k = 0; k < natoms; k++)
	      s += X[k][i] * sc[k*natoms+j];
	  g[i] = s / sumSC;
	}
      } else {
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (j = 0; j < natoms; j++)
	    s += X[j][i] * sc[j*natoms+j];
	  g[i] = s / sumSC;
	}
      }

      // Z == D - Y
      for (i = 0; i < natoms; i++) {
	Z[i][0] = D[i][0] - Y[i][0];
	Z[i][1] = D[i][1] - Y[i][1];
	Z[i][2] = D[i][2] - Y[i][2];
      }
      // ZTC == Z^T C,  (3 x n)
      if (fullV) {
	for (i = 0; i < 3; i++)
	  for (j = 0; j < natoms; j++) {
	    s = 0.0;
	    for (k = 0; k < natoms; k++)
	      s += Z[k][i] * cp[k*natoms+j];
	    ZTC[i*natoms+j] = s;
	  }
      } else {
	for (i = 0; i < 3; i++)
	  for (j = 0; j < natoms; j++)
	    ZTC[i*natoms+j] = Z[j][i] * cp[j*natoms+j];
      }

      /* W <- (X - 1 \otimes g) */
      for (i = 0; i < natoms; i++) {
        W[i][0] = X[i][0] - g[0];
        W[i][1] = X[i][1] - g[1];
        W[i][2] = X[i][2] - g[2];
      }

      /* A <- (MTS + ZTC) W */
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
          s = 0.0;
          for (k = 0; k < natoms; k++)
            s += (MTS[i*natoms+k] + ZTC[i*natoms+k]) * W[k][j];
          A[j][i] = s;
        }

      rotation_matrix_svd(ap, confp->r, 0);

      /* shift == R M^T S^{-1} 1.  これは centerM() しておればゼロのはず。 */

      if (Mcentered) {
	/* t <- g - R Z^T C 1 */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += R[i][0] * (ZTC[0*natoms+k])
	      + R[i][1] * (ZTC[1*natoms+k])
	      + R[i][2] * (ZTC[2*natoms+k]);
	  t[i] = g[i] - s / sumSC;
	}
      } else {		// if (!Mcentered)
	/* t <- g - R (M^T S + Z^T C) 1 */
	for (i = 0; i < 3; i++) {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += R[i][0] * (MTS[0*natoms+k] + ZTC[0*natoms+k])
	      + R[i][1] * (MTS[1*natoms+k] + ZTC[1*natoms+k])
	      + R[i][2] * (MTS[2*natoms+k] + ZTC[2*natoms+k]);
	  t[i] = g[i] - s / sumSC;
	}
      }
    }
  }
  return;
}

/*
 * Y = (X - 1 x t^T) R
 * missing data are imputed with expected values if op != NULL
 */
void
calculate_Y(Real *xp, int *op, Real *rp, Real *t, Real *mp, Real *zp,
            int natoms, Real *yp)
{
  int i;
  Real sx[3];
  Matrix3 M, X, R, Y, Z;

  X = (Matrix3) xp;
  R = (Matrix3) rp;
  Y = (Matrix3) yp;
  Z = (Matrix3) zp;
  M = (Matrix3) mp;
  for (i = 0; i < natoms; i++) {
    if (op && !op[i]) {
      Y[i][0] = M[i][0] + Z[i][0];
      Y[i][1] = M[i][1] + Z[i][1];
      Y[i][2] = M[i][2] + Z[i][2];
    } else {
      sx[0] = X[i][0] - t[0];
      sx[1] = X[i][1] - t[1];
      sx[2] = X[i][2] - t[2];
      Y[i][0] = sx[0]*R[0][0] + sx[1]*R[1][0] + sx[2]*R[2][0];
      Y[i][1] = sx[0]*R[0][1] + sx[1]*R[1][1] + sx[2]*R[2][1];
      Y[i][2] = sx[0]*R[0][2] + sx[1]*R[1][2] + sx[2]*R[2][2];
    }
  }
}



