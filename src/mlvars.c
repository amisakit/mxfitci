#include "atypes.h"
#include "matrix.h"
#include "molec.h"

/*
 * 
 */
Real
estimate_SV_MLold(Ensemble Bbuf[], Real *mp, int natoms,
	       int nconfms_tot, int nconfms[], int ngroups,
	       Real *diagS, Real *vp, Real *work)
{
  Real *sp, *dp, di[3], dj[3], s, sig0sq;
  Matrix3 M, Y;
  Ensemble *ensp;
  Conformer *confp;
  int iens, iconf, mi, i, j, ndim=3;
  
  sp = work;
  dp = sp + natoms*natoms;

  M = (Matrix3) mp;

  zerovectmatrix(sp, natoms, natoms);
  zerovectmatrix(vp, natoms, natoms);

  ensp = Bbuf;
  for (iens = 0; iens < ngroups; iens++, ensp++) {

    mi = ensp->nconfms;

    // D <- 0
    zerovectmatrix(dp, natoms, 3);

    // S <- S + sum_j (Y - M) (Y - M)^T
    // D <- sum_j (Y - M)
    confp = ensp->confmbp;
    for (iconf = 0; iconf < mi; iconf++, confp++) {
      Y = (Matrix3) confp->y;
      for (i = 0; i < natoms; i++) {
	di[0] = Y[i][0] - M[i][0];
	di[1] = Y[i][1] - M[i][1];
	di[2] = Y[i][2] - M[i][2];
	for (j = 0; j < natoms; j++) {
	  dj[0] = Y[j][0] - M[j][0];
	  dj[1] = Y[j][1] - M[j][1];
	  dj[2] = Y[j][2] - M[j][2];
	  s = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
	  sp[i*natoms + j] += s; 
	}
	dp[i*ndim+0] += di[0];
	dp[i*ndim+1] += di[1];
	dp[i*ndim+2] += di[2];
      }
    }
    
    // S <- S - (1/m1) (D D^T),  it must be mi == m1 
    // V <- V + (1/m1) (D D^T),  it must be mi == m1 
    for (i = 0; i < natoms; i++)
      for (j = 0; j < natoms; j++) {
	s = dp[i*ndim+0] * dp[j*ndim+0]
	  + dp[i*ndim+1] * dp[j*ndim+1]
	  + dp[i*ndim+2] * dp[j*ndim+2];
	sp[i*natoms+j] -= s / mi;
	vp[i*natoms+j] += s / mi;
      }
  }
    
  // S <- S / (3(m-r))
  for (i = 0; i < natoms*natoms; i++)
    sp[i] /= 3.0 * (nconfms_tot - ngroups);

  sig0sq = 0.0;
  for (i = 0; i < natoms; i++) {
    diagS[i] = sp[i*natoms+i];
    sig0sq += diagS[i];
  }
  sig0sq /= natoms;

  // V <- V / 3m - S * r / m,   it must be mi == m1 
  for (i = 0; i < natoms; i++)
    for (j = 0; j < natoms; j++)
      if (i == j)
	vp[i*natoms+j] = vp[i*natoms+j] / (3.0 * nconfms_tot)
	  - diagS[i] * ngroups / nconfms_tot;
      else
	vp[i*natoms+j] = vp[i*natoms+j] / (3.0 * nconfms_tot);

  return sig0sq;
}
