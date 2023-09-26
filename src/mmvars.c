#include <math.h>
#include <string.h>
#include "atypes.h"
#include "matrix.h"
#include "molec.h"
#include "mmfit.h"
#include "defs.h"

extern void abnormalexit(char *message, int code);
extern void gather_aniso(Real *dest, Real *src, int natoms, int diag);
extern void scatter_aniso(Real *dest, Real *src, int natoms, int diag);


/*
 * obtain the MAP estimate Ui and its auxiliary matrices for each i,
 * including the joint covariance matrix Sigma + mi V.
 * work area: natoms^2
 */
void
estimate_Uaux(struct jvar Jbuf[], Real V[], Real *diagV, Real diagS[], int natoms, Real *work)
{
  Real *up, *svp, *bp, *cp, *scp, s, w, *svip, *hp, *sumsvi;
  int i, j, k, mi, fullV, njvars, ngroups, nconfms_tot;
  struct jvar *jvp;

  sumsvi = work;

  fullV = (diagV == NULL);

  /*
  if (diagV)
    for (i = 0; i < natoms; i++)
      diagV[i] = V[i*natoms+i];
  */

  njvars = ngroups = nconfms_tot = 0;
  k = 0;
  jvp = Jbuf;
  while (jvp) {
    mi = jvp->mi;
    up = jvp->u;
    svp = jvp->sv;
    bp = jvp->b;
    if (fullV) {
      /* svp <- mi V + diagS */
      for (i = 0; i < natoms; i++) {
	for (j = 0; j < i; j++)
	  svp[i*natoms+j] = svp[j*natoms+i] = mi * V[i*natoms+j];
	svp[i*natoms+i] = mi * V[i*natoms+i] + diagS[i];
      }
      /* svp <- L s.t. L L^T == mi V + diagS */
      i = cholesky_decomp(svp, natoms);
      if (i != 0)
	abnormalexit("cholesky_decomp(svp)", k);
      /* log determinant */
      s = 0.0;
      for (i = 0; i < natoms; i++)
	s += log(svp[i*natoms+i]);
      jvp->logdet = 2.0 * s;
      /* solve (S + mi V) X = V for X. */
      memcpy(bp, V, sizeof(Real) * natoms * natoms);
      /* call with bp=V^T, returns with bp=X^T. note V == V^T.   */
      i = cholesky_solve(svp, bp, natoms, natoms);
      if (i != 0)
	abnormalexit("cholesky_solve(svp,bp)", k);
      /* up <- (X^T)_ij Sigma_j == (Sigma_i X_ij)^T == Ui^T. note U^T == U */
      for (i = 0; i < natoms; i++) {
	for (j = 0; j < i; j++)
	  up[i*natoms+j] = up[j*natoms+i]
	    = bp[i*natoms+j] * diagS[j];
	up[i*natoms+i] = bp[i*natoms+i] * diagS[i];
      }
    } else { // diagV
      /* svp <- sqrt(mi V + sig^2 diagS) */
      zerovectmatrix(svp, natoms, natoms);
      s = 0.0;
      for (i = 0; i < natoms; i++) {
	w = mi * diagV[i] + diagS[i];
	svp[i*natoms+i] = sqrt(w);
	s += log(w);
      }
      jvp->logdet = s;
      /* solve (S + mi V) X = V for X.  */
      zerovectmatrix(bp, natoms, natoms);
      for (i = 0; i < natoms; i++)
	bp[i*natoms+i] = diagV[i] / (diagS[i] + mi * diagV[i]);
      /* up <- (X^T)_ij Sigma_j == (Sigma_i X_ij)^T == Ui^T. note U^T == U */
      zerovectmatrix(up, natoms, natoms);
      for (i = 0; i < natoms; i++)
	up[i*natoms+i] = bp[i*natoms+i] * diagS[i];
    }

    /* Ci and S^{-1} - Ci */
    cp = jvp->c;
    scp = jvp->sc;
    if (cp) {
      zerovectmatrix(cp, natoms, natoms);
      zerovectmatrix(scp, natoms, natoms);
      if (fullV) {
        /* cp <- Sigma^{-1} X^T */
        for (i = 0; i < natoms; i++) {
          for (j = 0; j < i; j++)
            cp[i*natoms+j] = cp[j*natoms+i]
	      = bp[i*natoms+j] / diagS[i];
	  cp[i*natoms+i] = bp[i*natoms+i] / diagS[i];
	}
        /* scp <- Sigma^{-1} - C */
        for (i = 0; i < natoms*natoms; i++)
	  scp[i] = -cp[i];
        for (i = 0; i < natoms; i++)
          scp[i*natoms+i] += 1.0 / diagS[i];
      } else {
        /* cp <- Sigma^{-1} X^T */
        for (i = 0; i < natoms; i++)
          cp[i*natoms+i] = bp[i*natoms+i] / diagS[i];
        /* scp <- Sigma^{-1} - C */
        for (i = 0; i < natoms; i++)
          scp[i*natoms+i] = 1.0 / diagS[i] - cp[i*natoms+i];
      }
    }

    /* (S+miV)^{-1} */
    svip = jvp->svi;
    if (svip) {			// METHOD_FIX_ML or METHOD_VAR_REML
      zerovectmatrix(svip, natoms, natoms);
      if (fullV) {
	/* svip <- (S + mi V)^{-1} */
	i = cholesky_invert(svp, svip, natoms);
	if (i != 0)
	  abnormalexit("cholesky_invert(svp,svip)", k);
      } else {
	/* svip <- (S + mi V)^{-1} */
	for (i = 0; i < natoms; i++)
	  svip[i*natoms+i] = 1.0 / (diagS[i] + mi * diagV[i]);
      }
    }

    k++;
    ngroups += jvp->nens;
    nconfms_tot += jvp->mi * jvp->nens;
    jvp = jvp->next;
  }

  /* (\sum mi(S+miV)^{-1})^{-1}, i.e., weighted harmonic mean of (S+miV) m */
  njvars = k;		/* == nmi */
  hp = Jbuf[0].h;
  if (hp) {		// (METHOD_FIX_ML and njvars > 1) or METHOD_VAR_REML
    if (njvars == 1) {			/* balanced design */
      if (fullV) {
	// hp <- V / r
	for (i = 0; i < natoms*natoms; i++)
	  hp[i] = V[i] / ngroups;
	// hp <- V / r + S / m
	for (i = 0; i < natoms; i++)
	  hp[i*natoms+i] += diagS[i] / nconfms_tot;
      } else {
	zerovectmatrix(hp, natoms, natoms);
	for (i = 0; i < natoms; i++)
	  hp[i*natoms+i] = diagS[i] / nconfms_tot + diagV[i] / ngroups;
      }
    } else {				/* unbalanced design */
      if (fullV) {
	// sumsvi += mi_k * (S+mi_k V) * sizeof(k) for all mi-group k
	zerovectmatrix(sumsvi, natoms, natoms);
	for (jvp = Jbuf; jvp != NULL; jvp = jvp->next)
	  for (i = 0; i < natoms*natoms; i++)
	    sumsvi[i] += jvp->mi * jvp->svi[i] * jvp->nens;
	// sumsvi <- L s.t. L L^T = sum_i mi(S+miV)^{-1}
	i = cholesky_decomp(sumsvi, natoms);
	if (i != 0)
	  abnormalexit("cholesky_decomp(sumsvi)", 0);
	/* hp <- (sum_i mi(S+miV)^{-1})^{-1} */
	i = cholesky_invert(sumsvi, hp, natoms);
	if (i != 0)
          abnormalexit("cholesky_invert(sumsvi,hp)", 0);
      } else {
	zerovectmatrix(hp, natoms, natoms);
	for (i = 0; i < natoms; i++) {
	  s = 0.0;
	  for (jvp = Jbuf; jvp != NULL; jvp = jvp->next)
	    s += jvp->nens * jvp->mi / (diagS[i] + jvp->mi * diagV[i]);
	  hp[i*natoms+i] = 1.0 / s;
	}
      }
      //      for (jvp = Jbuf[0].next; jvp != NULL; jvp = jvp->next)
      for (jvp = Jbuf; jvp != NULL; jvp = jvp->next)
	memcpy(jvp->h, hp, sizeof(Real)*natoms*natoms);
    }
  }

  return;
}


/*
 * workarea: (natoms*3)^2 + natoms*3
 */
void
estimate_S(Real *Ybuf, Real *mp, Real *Zbuf, struct jvar Jbuf[],
	   int nconfms_tot, int nconfms[], int ngroups, int natoms, int method,
	   Real *diagV, Real diagS[], Real *sp, Real *ani, Real *work)
{
  Real *var, *up, *wp, *svi, *hp, *yij, *zi, s, *ep, *reml;
  int i, j, k, iens, iconfm, fullV, mi, nens, n3, ii, jj;
  int *ixm;
  Matrix3 M, Z, Y;
  struct jvar *jvp;

  ixm = nconfms + 3*ngroups + 1;	/* ixm[iens] is the index for Jbuf */

  n3 = natoms * 3;

  var = sp;
  wp = work;
  ep = wp + natoms*natoms;
  reml = ep + natoms * 3;

  fullV = (diagV == NULL);

  M = (Matrix3) mp;
  yij = Ybuf;
  zi = Zbuf;

  if (ani)
    zerovectmatrix(ani, n3, n3);
  else
    zerovectmatrix(var, natoms, natoms);

  for (iens = 0; iens < ngroups; iens++) {
    up = Jbuf[ixm[iens]].u;
    if (ani) {
      //scatter_aniso(wp, up, natoms, 0);
      for (iconfm = 0; iconfm < nconfms[iens]; iconfm++) {
	Y = (Matrix3) yij;
	for (i = 0; i < n3; i++)
	  ep[i] = yij[i] - mp[i] - zi[i];
	for (i = 0; i < n3; i++)
	  for (j = 0; j <= i; j++) {
	    ii = i / 3;
	    jj = j / 3;
	    ani[i*n3+j] += up[ii*natoms+jj] + ep[i] * ep[j];
	  }
	yij += natoms * 3;
      }
    } else {
      Z = (Matrix3) zi;
      for (iconfm = 0; iconfm < nconfms[iens]; iconfm++) {
	Y = (Matrix3) yij;
	for (i = 0; i < natoms; i++) {
	  for (j = 0; j <= i; j++) {
	    var[i*natoms+j] +=
	      (Y[i][0] - M[i][0] - Z[i][0]) * (Y[j][0] - M[j][0] - Z[j][0])
	      + (Y[i][1] - M[i][1] - Z[i][1]) * (Y[j][1] - M[j][1] - Z[j][1])
	      + (Y[i][2] - M[i][2] - Z[i][2]) * (Y[j][2] - M[j][2] - Z[j][2]);
	    var[i*natoms+j] += 3.0 * up[i*natoms+j]; 
	  }
	}
	yij += natoms * 3;
      }
    }
    zi += natoms * 3;
  }

  if (ani)
    gather_aniso(var, ani, natoms, 0);

  /* * mi を加えた。要検証。[2022-08-03 Wed 08:52] */
  if ((method & METHOD_VAR) == METHOD_VAR_REML) {
    zerovectmatrix(reml, natoms, natoms);
    if (fullV) {
      jvp = Jbuf;
      while (jvp) {
	svi = jvp->svi;
        hp = jvp->h;
	mi = jvp->mi;
	nens = jvp->nens;
        /* W <- S (S+miV)^{-1} H_i */
        for (i = 0; i < natoms; i++)
          for (j = 0; j < natoms; j++) {
            s = 0.0;
            for (k = 0; k < natoms; k++)
              s += diagS[i] * svi[i*natoms+k] * hp[k*natoms+j];
            wp[i*natoms+j] = s;
          }
        /* var += 3 mi W (S+miV)^{-1} S */
        for (i = 0; i < natoms; i++)
          for (j = 0; j <= i; j++) {
            s = 0.0;
            for (k = 0; k < natoms; k++)
              s += wp[i*natoms+k] * svi[k*natoms+j] * diagS[j];
            s = nens * mi * s;
            var[i*natoms+j] += 3.0 * s;
            reml[i*natoms+j] += s;
          }
	jvp = jvp->next;
      }
    } else {			// diagV
      jvp = Jbuf;
      while (jvp) {
	svi = jvp->svi;
        hp = jvp->h;
	mi = jvp->mi;
	nens = jvp->nens;
        /* W <- S (S+miV)^{-1} H_i */
        for (i = 0; i < natoms; i++)
          wp[i*natoms+i] = diagS[i] * svi[i*natoms+i] * hp[i*natoms+i];
        /* var += 3 mi W (S+miV)^{-1} S */
        for (i = 0; i < natoms; i++) {
          s = nens * mi * wp[i*natoms+i] * svi[i*natoms+i] * diagS[i];
          var[i*natoms+i] += 3.0 * s;
          reml[i*natoms+i] += s;
	}
	jvp = jvp->next;
      }
    }
    if (ani) {
      if (fullV) {
        for (i = 0; i < n3; i++) 
          for (j = 0; j <= i; j++) {
            ii = i / 3;
            jj = j / 3;
            ani[i*n3+j] += reml[ii*natoms+jj];
	  }
      } else {
	for (i = 0; i < n3; i++) {
          ii = i / 3;
          ani[i*n3+i] += reml[ii*natoms+ii];
        }
      }
    }
  }

  for (i = 0; i < natoms; i++) {
    for (j = 0; j < i; j++) {
      var[i*natoms+j] /= 3.0 * nconfms_tot;
      var[j*natoms+i] = var[i*natoms+j];
    }
    var[i*natoms+i] /= 3.0 * nconfms_tot;
  }

  for (i = 0; i < natoms; i++)
    diagS[i] = var[i*natoms+i];

  if (ani)
    for (i = 0; i < n3; i++)
      for (j = 0; j <= i; j++) {
	ani[i*n3+j] /= nconfms_tot;
	ani[j*n3+i] = ani[i*n3+j];
      }

  return;
}

/*
 * workarea: natoms^2 x 2
 * Output: V on vp
 */
int
estimate_V(Real *Zbuf, struct jvar Jbuf[], int nconfms[], int ngroups, int natoms, int method, Real *diagV, Real *vp, Real *ani, Real *work)
{
  Real *up, *wp, *bp, *hp, s, *zi, *reml;
  Matrix3 Z;
  int iens, i, j, k, kk, fullV, mi, *ixm, nens, n3, ii, jj;
  struct jvar *jvp;

  ixm = nconfms + 3*ngroups + 1;	/* ixm[iens] is the index for Jbuf */

  n3 = natoms * 3;

  wp = work;
  reml = wp + natoms * natoms;

  fullV = (diagV == NULL);

  if (ani)
    zerovectmatrix(ani, n3, n3);
  else
    zerovectmatrix(vp, natoms, natoms);

  zi = Zbuf;
  for (iens = 0; iens < ngroups; iens++) {
    // V += Zi Zi^T + 3 Ui
    up = Jbuf[ixm[iens]].u;
    if (ani) {
      for (i = 0; i < n3; i++) 
	for (j = 0; j <= i; j++) {
	  ii = i / 3;
	  jj = j / 3;
	  ani[i*n3+j] += up[ii*natoms+jj] + zi[i] * zi[j];
	}
    } else {
      Z = (Matrix3) zi;
      for (i = 0; i < natoms; i++) {
	for (j = 0; j <= i; j++)
	  vp[i*natoms+j] += 3.0 * up[i*natoms+j]
	    + Z[i][0]*Z[j][0] + Z[i][1]*Z[j][1] + Z[i][2]*Z[j][2];
      }
    }
    zi += natoms * 3;
  }

  if (ani)
    gather_aniso(vp, ani, natoms, 0);

  if ((method & METHOD_VAR) == METHOD_VAR_REML) {
    zerovectmatrix(reml, natoms, natoms);
    if (fullV) {
      jvp = Jbuf;
      while (jvp) {
        bp = jvp->b;
        hp = jvp->h;
        mi = jvp->mi;
        nens = jvp->nens;
	/* W <- B_i H_i */
	for (i = 0; i < natoms; i++)
	  for (j = 0; j < natoms; j++) {
	    s = 0.0;
	    for (k = 0; k < natoms; k++)
	      s += bp[i*natoms+k] * hp[k*natoms+j];
	    wp[i*natoms+j] = s;
	  }
	/* V += W B_i^T */
	for (i = 0; i < natoms; i++)
	  for (j = 0; j <= i; j++) {
	    s = 0.0;
	    for (k = 0; k < natoms; k++)
	      s += wp[i*natoms+k] * bp[j*natoms+k];
	    s *= nens * mi * mi;
	    vp[i*natoms+j] += 3.0 * s;
	    reml[i*natoms+j] += s;
	  }
	jvp = jvp->next;
      }
    } else {		// diagV
      jvp = Jbuf;
      while (jvp) {
        bp = jvp->b;
        hp = jvp->h;
        mi = jvp->mi;
        nens = jvp->nens;
	/* W <- B_i H_i */
	for (i = 0; i < natoms; i++)
	  wp[i*natoms+i] = bp[i*natoms+i] * hp[i*natoms+i];
	/* V += W B_i^T */
	for (i = 0; i < natoms; i++) {
	  s = nens * mi * mi * wp[i*natoms+i] * bp[i*natoms+i];
	  vp[i*natoms+i] += 3.0 * s;
	  reml[i*natoms+i] += s;
	}
	jvp = jvp->next;
      }
    }
    if (ani) {
      if (fullV) {
	for (i = 0; i < n3; i++) 
	  for (j = 0; j <= i; j++) {
	    ii = i / 3;
	    jj = j / 3;
	    ani[i*n3+j] += reml[ii*natoms+jj];
	  }
      } else {
	for (i = 0; i < n3; i++) {
	  ii = i / 3;
	  ani[i*n3+i] += reml[ii*natoms+ii];
	}
      }
    }
  }

  for (i = 0; i < natoms; i++)
    for (j = 0; j <= i; j++) {
      vp[i*natoms+j] /= (3.0 * ngroups);
      vp[j*natoms+i] = vp[i*natoms+j];
    }
  if (diagV)
    for (i = 0; i < natoms; i++)
      diagV[i] = vp[i*natoms+i];
  if (ani)
    for (i = 0; i < n3; i++)
      for (j = 0; j <= i; j++) {
	ani[i*n3+j] /= ngroups;
	ani[j*n3+i] = ani[i*n3+j];
      }

  return 0;
}



/*
 * 
 */
Real
estimate_SV_BML(Real Ybuf[], Real *mp, int natoms, int nconfms_tot,
	       int nconfms[], int ngroups, Real *diagS, Real *vp, Real *work)
{
  Real *sp, *dp, dk[3], s, sig0sq, *dij, *yp;
  int iens, iconf, mi, k, l, ndim=3;
  
  sp = work;
  dp = sp + natoms*natoms;
  dij = dp + natoms * 3;
  
  zerovectmatrix(sp, natoms, natoms);
  zerovectmatrix(vp, natoms, natoms);

  yp = Ybuf;
  for (iens = 0; iens < ngroups; iens++) {

    mi = nconfms[iens];

    // D <- 0
    zerovectmatrix(dp, natoms, 3);
    // S <- S + sum_j (Y - M) (Y - M)^T
    // D <- sum_j (Y - M)
    for (iconf = 0; iconf < mi; iconf++) {

      // Dij <- Yij - M
      // D <- D + Dij
      for (k = 0; k < natoms*3; k++) {
        dij[k] = yp[k] - mp[k];
	dp[k] += dij[k];
      }
      // S <- S + (Yij - M) (Yij - M)^T
      for (k = 0; k < natoms; k++) {
	dk[0] = dij[k*ndim+0];
	dk[1] = dij[k*ndim+1];
	dk[2] = dij[k*ndim+2];
	for (l = 0; l < k; l++) {
	  s = dk[0]*dij[l*ndim+0] + dk[1]*dij[l*ndim+1] + dk[2]*dij[l*ndim+2];
	  sp[k*natoms + l] += s; 
	}
	sp[k*natoms + k] += dk[0]*dk[0] + dk[1]*dk[1] + dk[2]*dk[2];
      }
      
      yp += natoms * 3;
    }
    
    // S <- S - (1/m1) (D D^T),  it must be mi == m1 
    // V <- V + (1/m1) (D D^T),  it must be mi == m1 
    for (k = 0; k < natoms; k++) {
      for (l = 0; l < k; l++) {
	s = dp[k*ndim+0] * dp[l*ndim+0]
	  + dp[k*ndim+1] * dp[l*ndim+1]
	  + dp[k*ndim+2] * dp[l*ndim+2];
	sp[k*natoms+l] -= s / mi;
	vp[k*natoms+l] += s / mi;
      }
	s = dp[k*ndim+0] * dp[k*ndim+0]
	  + dp[k*ndim+1] * dp[k*ndim+1]
	  + dp[k*ndim+2] * dp[k*ndim+2];
	sp[k*natoms+k] -= s / mi;
	vp[k*natoms+k] += s / mi;
    }
  }
    
  // S <- S / (3(m-r))
  s = 1.0 / (3.0 * (nconfms_tot - ngroups));
  for (k = 0; k < natoms*natoms; k++)
    sp[k] *= s;

  sig0sq = 0.0;
  for (k = 0; k < natoms; k++) {
    for (l = 0; l < k; l++)
      sp[l*natoms+k] = sp[k*natoms+l];
    diagS[k] = sp[k*natoms+k];
    sig0sq += diagS[k];
  }
  sig0sq /= natoms;

  // V <- V / 3m - S * r / m,   it must be mi == m1 
  for (k = 0; k < natoms; k++) {
    for (l = 0; l < k; l++)
      vp[k*natoms+l] = vp[l*natoms+k]
	= vp[k*natoms+l] / (3.0 * nconfms_tot);
    vp[k*natoms+k] = vp[k*natoms+k] / (3.0 * nconfms_tot)
      - diagS[k] * ngroups / nconfms_tot;
  }

  return sig0sq;
}
