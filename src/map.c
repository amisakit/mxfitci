#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atypes.h"
#include "molec.h"
#include "matrix.h"

/*
 * estimates Z_i (using U_i)
 * work: natoms x 3
 */
void
estimate_Z(Ensemble Bbuf[], Real *diagV, Real *mp, int natoms, int ngroups,
	   Real *work)
{
  int i, j, k, iens, iconfm, nconfms;
  const int ndim=3;
  Real s, *yij, *dp, *bp, *zp;

  dp = work;
  for (iens = 0; iens < ngroups; iens++) {
    /* dp(Di) <- \sum_j (Y_{ij} - M) */
    zerovectmatrix(dp, natoms, ndim);
    yij = Bbuf[iens].y;
    nconfms = Bbuf[iens].nconfms;
    for (iconfm = 0; iconfm < nconfms; iconfm++) {
      for (i = 0; i < natoms * ndim; i++)
	dp[i] += yij[i] - mp[i];
      yij += natoms * ndim;
    }

    /* Z <- Bi Di */
    bp = Bbuf[iens].jvp->b;
    zp = Bbuf[iens].z;
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < ndim; j++) {
	if (diagV) {
	  zp[i*ndim+j] = bp[i*natoms+i] * dp[i*ndim+j];
	} else {
	  s = 0.0;
	  for (k = 0; k < natoms; k++)
	    s += bp[i*natoms+k] * dp[k*ndim+j];
	  zp[i*ndim+j] = s;
	}
      }
    }
  }
  return;
}
