#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include "atypes.h"
#include "matrix.h"
#include "msfit.h"
#include "lapack.h"
#include "memory.h"

  
/*
 * complete-case. (but returns pairwise covar on incomplete-case)
 * calculates sample covariances, and writes
 * Saniso[0..n3^2-1] and Siso[0..natoms^2-1] into sp and siso, respectively.
 * Saniso == sum vec(Yj^T) vec(Yj^T)^T / N
 * Siso   == sum Yj Yj^T / N
 * cp->yp points Ybuf which contains vec(Yj^T), j=1,...,N. (xyzxyzxyz...)
 * work ... workarea of size natoms^2+natoms*3
 */
int
sample_covar(Conformer *cp, Real *mp, int natoms, int nconfms,
	     Real *sp, Real *siso, Real *work)
{
  int n3, i, j, k, l, iconfm, *op;
  Real *dp, *wp, *yp;

  n3 = natoms * 3;
  wp = work;
  dp = wp + natoms * natoms;

  zerovectmatrix(sp, n3, n3);
  zerovectmatrix(siso, natoms, natoms);
  zerovectmatrix(wp, natoms, natoms);
  
  for (iconfm = 0; iconfm < nconfms; iconfm++, cp++) {
    yp = cp->y;
    op = cp->o;

    for (i = 0; i < n3; i++)
      dp[i] = yp[i] - mp[i];
    for (i = 0; i < natoms; i++) {
      if (op == NULL || op[i] == 0)
	continue;
      for (j = 0; j < natoms; j++) {
	if (op == NULL || op[j] == 0)
	  continue;
	wp[i*natoms+j] += 1;
	for (k = 0; k < 3; k++)
	  for (l = 0; l < 3; l++)
	    sp[(i*3+k)*n3+j*3+l] += dp[i*3+k] * dp[j*3+l];
	siso[i*natoms+j] += dp[i*3+0] * dp[j*3+0]
	  + dp[i*3+1] * dp[j*3+1] + dp[i*3+2] * dp[j*3+2];
      }
    }
  }
  /*
  if (impute0)
    for (i = 0; i < n3 * n3; i++)
      sp[i] /= sumgam;
  else
  */
  for (i = 0; i < natoms* natoms; i++)
    siso[i] /= 3.0 * wp[i];
  for (i = 0; i < natoms; i++)
    for (j = 0; j < natoms; j++)
      for (k = 0; k < 3; k++)
	for (l = 0; l < 3; l++)
	  sp[(i*3+k)*n3+j*3+l] /= wp[i*natoms+j];

  return 0;
}



/*
 * Output:
 *   Saniso[0..(3*natoms)^-1]   ( hat(Omega) == Saniso )
 *   Siso[0..(natoms)^-1]  ( hat(Omega) == Siso \otimes I_3 )
 * work ... workarea of size natoms^2+natoms*3
 */
int
estimate_S(Conformer *cp, Real *mp, int natoms, int nconfms,
	   Real *diagS, Real *saniso, Real *siso, Real *work)
{
  int i;

  sample_covar(cp, mp, natoms, nconfms, saniso, siso, work);
  for (i = 0; i < natoms; i++)
    diagS[i] = siso[i*natoms+i];

  return 0;
}
