#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
       #include <sys/types.h>
       #include <sys/stat.h>
       #include <fcntl.h>
#include <unistd.h>
#include "atypes.h"
#include "matrix.h"
#include "defs.h"



void
usage()
{
  fprintf(stderr, "Usage: anisovar [-b] [-z] [-e eps] [-n niters] [-d] [-r] iin din dout\n"
	  "iin din: binary files from mmfit\n"
	  "dout: estimates for V & S (binary file)\n"
	  "-b: in balanced-designs, use dedicated estimators\n"
	  "-d: ignore off-diagonal terms of V\n"
	  "-r: REML (default: ML)\n"
	  "-z: update Zi (default: do not update Zi)\n"
	  
	  );
  exit(-1);
}

void abnormalexit(char *message, int code)
{
  if (message && code)
    fprintf(stderr, "(%d) %s\n", code, message);
  else if (message)
    fprintf(stderr, "%s\n", message);
  exit(code);
}

void print_binary(unsigned int n)
{
    if (n >> 1) {
        print_binary(n >> 1);
    }
    putc((n & 1) ? '1' : '0', stdout);
}

int
parse_args(int argc, char *argv[],
	   char **iinfname, char **dinfname, char **outfname,
	   int *niters, Real *eps, int *updateZi, unsigned int *method)
{
  int opt;

  *niters = 1000;
  *eps = 0.5e-6;
  *updateZi = 0;
  while ((opt = getopt(argc, argv, "bde:n:rz")) != -1)
    switch (opt) {
    case 'b':
      *method &= ~METHOD_D;
      *method |= METHOD_D_BALANCED;
      break;
    case 'd':
      *method &= ~METHOD_V;
      *method |= METHOD_V_DIAG;
      break;
    case 'e':
      *eps = atof(optarg);
      break;
    case 'n':
      *niters = atoi(optarg);
      break;
    case 'r':
      *method &= ~METHOD_VAR;
      *method |= METHOD_VAR_REML;
      break;
    case 'z':
      *updateZi = 1;
      break;
    default:
      usage();
    }
  if ((argc - optind) != 3)
    usage();
  *iinfname = strdup(argv[optind++]);
  *dinfname = strdup(argv[optind++]);
  *outfname = strdup(argv[optind++]);
  return 0;
}


int main(int argc, char *argv[])
{
  char *iinfname, *dinfname, *outfname;
  int fd, num_atoms, ngroups, *nconfms, nconfms_tot, iens, nmi, mi;
  int *nc, *ix, *ixm, *ne;
  int i, j, k, n, iter, niters;
  Real *wp, *vp, *sp, *vpp, *spp, *yp, *mp, *zp, *up, *dp, *ep, *svp;
  Real *Dbuf, *Ebuf, *Mbuf, *Ubuf, *Ybuf, *Zbuf, *covbuf, *crdbuf;
  Real *isovp, *isosp;
  Real d, s, s1, eps;
  Real *svi, *bvp, *bsp, *hp, *bp;
  int fullV, kmi, updateZi;
  unsigned int method;

  parse_args(argc, argv, &iinfname, &dinfname, &outfname,
	     &niters, &eps, &updateZi, &method);
  fullV = (method & METHOD_V) != METHOD_V_DIAG;


  /*
   * int params
   */
  fd = open(iinfname, O_RDONLY);
  if (fd < 0) {
    fprintf(stderr, "cannot open %s.\n", iinfname);
    exit(-1);
  }
  read(fd, &num_atoms, sizeof(int));

  n = num_atoms * 3;

  read(fd, &ngroups, sizeof(int));
  nconfms = (int *) malloc(sizeof(int) * (ngroups * 5 + 1));
  read(fd, nconfms, sizeof(int) * ngroups);
  close(fd);
  nconfms_tot = 0;
  for (iens = 0; iens < ngroups; iens++)
    nconfms_tot += nconfms[iens];

  /*
   * Given V and S, Ui and related matrices are only functions of mi.
   * The matrices are prepared for each distinct mi.
   */
                                // nc[-1] == num of different mi's, nmi
  nc = nconfms + ngroups + 1;   // nc[] <- unique(sort(nconfms[]))
  ix = nc + ngroups;            // ix[] s.t. sort(nconfms)[k] == nconfms[ix[k]]
  ixm = ix + ngroups;           // ixm[] s.t. nconfoms[k] == nc[ixm[k]]
  ne = ixm + ngroups;		// ne[k] == # of ensembles with size of nc[k] 

  for (iens = 0; iens < ngroups; iens++) {
    nc[iens] = nconfms[iens];
    ix[iens] = iens;
  }
  for (i = 0; i < ngroups; i++)
    for (j = i+1; j < ngroups; j++)
      if (nc[j] < nc[i]) {
        mi = nc[j];
        nc[j] = nc[i];
        nc[i] = mi;
        mi = ix[j];
        ix[j] = ix[i];
        ix[i] = mi;
      }
  j = 0;
  for (i = 0; i < ngroups; i++) {
    if (nc[i] > nc[j])
      nc[++j] = nc[i];
    ixm[ix[i]] = j;
  }
  nmi = j + 1;          // num of different mi's (== # of S+mi*V)
  nconfms[ngroups] = nmi;
  for (i = 0; i < nmi; i++)
    ne[i] = 0;
  for (i = 0; i < ngroups; i++)
    ne[ixm[i]]++;

  if (nmi > 1 && (method & METHOD_D) == METHOD_D_BALANCED) {
    method &= ~METHOD_D;
    method |= METHOD_D_UNBALANCED;
    printf("option for balanced-estimator is turned off since nmi = %d",nmi);
  }
  printf("method = ");
  print_binary(method);
  printf("\n");

  /*
  printf("      mi   ix  ixm\n");
  for (i = 0; i < ngroups; i++)
    printf("%3d %4d %4d %4d\n", i, nconfms[i], ix[i], ixm[i]);
  printf("\n");
  */
  printf("      mi  # ensembles\n");
  for (i = 0; i < nmi; i++)
    printf("%3d %4d %4d\n", i, nc[i], ne[i]);


  /*
   * the coordinates and covariances
   */
  fd = open(dinfname, O_RDONLY);
  if (fd < 0) {
    fprintf(stderr, "cannot open %s.\n", dinfname);
    exit(-1);
  }

  /* coordinates */
  crdbuf = (Real *) malloc(sizeof(Real)
			   * ( 1		// M
			       + nconfms_tot	// Yij
			       + ngroups * 2	// Zi, Di,
			       + 1		// E
			       ) * n);
  Mbuf = crdbuf;
  Ybuf = Mbuf + n;
  Zbuf = Ybuf + n * nconfms_tot;
  Dbuf = Zbuf + n * ngroups;
  Ebuf = Dbuf + n * ngroups;
  i = sizeof(Real) * n;
  if (read(fd, Mbuf, i) < i)
    fprintf(stderr, "Mbuf\n"), exit(-1);
  if (read(fd, Ybuf, i*nconfms_tot) < i*nconfms_tot)
    fprintf(stderr, "Ybuf\n"), exit(-1);
  if (read(fd, Zbuf, i*ngroups) < i*ngroups)
    fprintf(stderr, "Zbuf\n"), exit(-1);
 
  /* covariances */
  covbuf = (Real *)
    malloc(sizeof(Real) *
	   ((2) * num_atoms*num_atoms	// iso V,S
	    + ( 4			// V,S, prev V,S
		+ 1			// wp
		+ 1		// (sum mi (S+miV)^{-1})^{-1}
		+ nmi * 6		// Ui,S+miV,V(S+miV)^{-1}
		// (S+miV)^{-1},mi(S+miV)^{-1}V,(S+miV)^{-1}S
		) * n * n));
  vp = covbuf + 2*num_atoms*num_atoms;
  sp = vp + n*n;
  vpp = sp + n*n;
  spp = vpp + n*n;
  wp = spp + n*n;
  hp = wp + n*n;
  Ubuf = hp + n*n;
  i = sizeof(Real) * num_atoms * num_atoms;
  if (read(fd, &covbuf[0], i) < i)
    fprintf(stderr, "iso V\n"), exit(-1);
  if (read(fd, &covbuf[num_atoms*num_atoms], i) < i)
    fprintf(stderr, "iso S\n"), exit(-1);

  i = sizeof(Real) * n * n;
  j = read(fd, &vp[0], i);
  k = read(fd, &sp[0], i);

  close(fd);

  if (j != i || k != i) {
    printf("Starting from isotropic matrices...\n");
    /*
     * expand isotropic variances over diagonal
     */
    isovp = covbuf;
    isosp = covbuf + 1 * num_atoms*num_atoms;
    zerovectmatrix(vp, n, n);
    zerovectmatrix(sp, n, n);
    for (i = 0; i < num_atoms; i++) {
      vp[(i*3+0)*n+(i*3+0)] = vp[(i*3+1)*n+(i*3+1)] = vp[(i*3+2)*n+(i*3+2)]
	= isovp[i*num_atoms+i];
      sp[(i*3+0)*n+(i*3+0)] = sp[(i*3+1)*n+(i*3+1)] = sp[(i*3+2)*n+(i*3+2)]
	= isosp[i*num_atoms+i];
    }
  } else
    printf("Starting from anisotropic matrices...\n");


  /*
   * Di <- sum_j (Yij-M)
   */
  dp = Dbuf;
  mp = Mbuf;
  yp = Ybuf;
  for (iens = 0; iens < ngroups; iens++) {
    mi = nconfms[iens];
    zerovectmatrix(dp, n, 1);
    for (j = 0; j < mi; j++) {
      for (k = 0; k < n; k++)
	dp[k] += yp[k] - mp[k];
      yp += n;
    }
    dp += n;
  }


  /* main loop */
  for (iter = 0; iter < niters; iter++) {

    /*
     * update Ui and related matrices
     * Off-diagonal components of V are ignored if `-d` given
     */
    for (kmi = 0; kmi < nmi; kmi++) {
      mi = nc[kmi];
      up = Ubuf + n*n * kmi;
      svp = up + n*n * nmi;
      bp =  up + n*n * nmi * 2;

      if (fullV) {
	/* svp <- mi V + diagS */
	for (i = 0; i < n; i++) {
	  for (j = 0; j < i; j++)
	    svp[i*n+j] = svp[j*n+i] = mi * vp[i*n+j];
	  svp[i*n+i] = mi * vp[i*n+i] + sp[i*n+i];
	}

	/* svp <- L s.t. L L^T == mi V + diagS */
	i = cholesky_decomp(svp, n);
	if (i != 0)
	  fprintf(stderr, "cholesky_decomp(svp) %d\n", iens), exit(-1);
	/* solve (S + mi V) X = V^T for X^T in wp. note V = V^T.  */
	memcpy(bp, vp, sizeof(Real) * n * n);
	i = cholesky_solve(svp, bp, n, n);
	if (i != 0)
	  fprintf(stderr, "cholesky_solve(svp,bp) %d\n", iens), exit(-1);
	/* now, bp == [ (S + mi V)^{-1} V^T ]^T = X^T */
	/* up <- sum_k (X^T)ik Skj = (X^T)ij Sjj since diagS */
	/* Ui symmetric @@ */
	for (i = 0; i < n; i++) {
	  for (j = 0; j < i; j++) {
	    up[i*n+j] = up[j*n+i] = bp[i*n+j] * sp[j*n+j];
	  }
	  up[i*n+i] = bp[i*n+i] * sp[i*n+i];
	// up[i*n+i] = vp[i*n+i] * sp[i*n+i] / (sp[i*n+i] + mi*vp[i*n+i]);
	}
      } else {
	for (i = 0; i < n; i++) {
	  for (j = 0; j < i; j++) {
	    up[i*n+j] = up[j*n+i] = 0.0;
	    bp[i*n+j] = bp[j*n+i] = 0.0;
	  }
	  up[i*n+i] = vp[i*n+i] * sp[i*n+i] / (sp[i*n+i] + mi*vp[i*n+i]);
	  bp[i*n+i] = vp[i*n+i] / (sp[i*n+i] + mi*vp[i*n+i]);
	}
      }
    }

    /*
     * update zi <- V (S + mi V)^{-1} di == bi di
     */
    if (updateZi) {
      dp = Dbuf;
      yp = Ybuf;
      zp = Zbuf;
      for (iens = 0; iens < ngroups; iens++) {
	mi = nconfms[iens];
	up = Ubuf + n*n * ixm[iens];
	bp =  up + n*n * nmi * 2;
	for (i = 0; i < n; i++) {
	  s = 0.0;
	  for (k = 0; k < n; k++)
	    s += bp[i*n+k] * dp[k];
	  zp[i] = s;
	}
	zp += n;
	dp += n;
      }
    }

    /*
     * update REML-specific matrices.
     * Off-diagonal components of V are ignored if '-d' given.
     */
    if ((method & METHOD_VAR) == METHOD_VAR_REML) {
      zerovectmatrix(hp, n, n);
      for (kmi = 0; kmi < nmi; kmi++) {
        mi = nc[kmi];
        up = Ubuf + n*n * kmi;
        svp = up + n*n*nmi * 1;
        bp  = up + n*n*nmi * 2;
        svi = up + n*n*nmi * 3;
        bvp = up + n*n*nmi * 4;
        bsp = up + n*n*nmi * 5;
	
	/* svi <- (S + mi V)^{-1} */
	if (fullV) {
	  i = cholesky_invert(svp, svi, n);
	  if (i != 0)
	    abnormalexit("cholesky_invert(svp,svi)", iens);
	} else {
	  zerovectmatrix(svi, n, n);
	  for (i = 0; i < n; i++)
	    svi[i*n+i] = 1.0 / (sp[i*n+i] + mi * vp[i*n+i]);
	}

	/* bvp == mi (S + mi V)^{-1} V == mi Bi^T  */
	if (fullV) {
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
	      bvp[i*n+j] = mi * bp[j*n+i];
	} else {
	  zerovectmatrix(bvp, n, n);
	  for (i = 0; i < n; i++)
	    bvp[i*n+i] = mi * vp[i*n+i] / (sp[i*n+i] + mi * vp[i*n+i]);
	}

	/* bsp <- (S + mi V)^{-1} S */
	if (fullV) {
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
	      bsp[i*n+j] = svi[i*n+j] * sp[j*n+j];
	} else {
	  zerovectmatrix(bsp, n, n);
	  for (i = 0; i < n; i++)
	    bsp[i*n+i] = sp[i*n+i] / (sp[i*n+i] + mi * vp[i*n+i]);
	}

	/* hp <- hp + mi * svi * nens with that mi */
	if (fullV) {
	  for (i = 0; i < n*n; i++)
	    hp[i] += mi * svi[i] * ne[kmi];
	} else {
	  for (i = 0; i < n; i++)
	    hp[i*n+i] += mi * svi[i*n+i] * ne[kmi];
	}
      }
      
      if ((method & METHOD_D) == METHOD_D_BALANCED && nmi == 1) {
	if (fullV)
	  for (i = 0; i < n; i++) {
	    for (j = 0; j < i; j++)
	      hp[i*n+j] = hp[j*n+i] = vp[i*n+j] / ngroups;
	    hp[i*n+i] = sp[i*n+i] / nconfms_tot + vp[i*n+i] / ngroups;
	  }
	else
	  for (i = 0; i < n; i++)
	    hp[i*n+i] = sp[i*n+i] / nconfms_tot + vp[i*n+i] / ngroups;
      } else {
	// hp == sum_i^r mi (S+miV)^{-1}, currently.
	/* hp <- (sum mi (S+miV)^{-1})^{-1} */
	if (fullV) {
	  memcpy(wp, hp, sizeof(Real) * n * n);
	  i = cholesky_decomp(wp, n);
	  if (i != 0)
	    abnormalexit("cholesky_decomp(wp)", 0);
	  /* hp <- (sum_i mi(S+miV)^{-1})^{-1} */
	  i = cholesky_invert(wp, hp, n);
	  if (i != 0)
	    abnormalexit("cholesky_invert(wp,hp)", 0);
	} else
	  for (i = 0; i < n; i++)
	    hp[i*n+i] = 1.0 / hp[i*n+i];
      }
    }


    memcpy(vpp, vp, sizeof(Real) * n * n);
    memcpy(spp, sp, sizeof(Real) * n * n);

    memset(vp, 0, sizeof(Real) * n * n);
    memset(sp, 0, sizeof(Real) * n * n);

    dp = 0;
    ep = Ebuf;
    mp = Mbuf;
    yp = Ybuf;
    zp = Zbuf;
    for (iens = 0; iens < ngroups; iens++) {
      mi = nconfms[iens];
      up = Ubuf + n*n * ixm[iens];
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  vp[i*n+j] += up[i*n+j] + zp[i] * zp[j];
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  sp[i*n+j] += mi * up[i*n+j];
      for (k = 0; k < mi; k++) {
	for (i = 0; i < n; i++)
	  ep[i] = yp[i] - mp[i] - zp[i];
	for (i = 0; i < n; i++)
	  for (j = 0; j < n; j++)
	    sp[i*n+j] += ep[i] * ep[j];
	yp += n;
      }
      zp += n;
    }

    if ((method & METHOD_VAR) == METHOD_VAR_REML) {
      for (kmi = 0; kmi < nmi; kmi++) {
        mi = nc[kmi];
        up = Ubuf + n*n * kmi;
        svp = up + n*n*nmi * 1;
	bp  = up + n*n*nmi * 2;
        svi = up + n*n*nmi * 3;
        bvp = up + n*n*nmi * 4;
        bsp = up + n*n*nmi * 5;

        // wp <- H Bv == (sum mi(s+miV)^{-1})^{-1} mi(s+miV)^{-1} V
	if (fullV)
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++) {
	      s = 0.0;
	      for (k = 0; k < n; k++)
		s += hp[i*n+k] * bvp[k*n+j];
	      wp[i*n+j] = s;
	    }
	else
	  for (i = 0; i < n; i++)
	    wp[i*n+i] = hp[i*n+i] * bvp[i*n+i];
        // vp <- (vp + Bv^T H Bv) * nens with that mi
	if (fullV) {
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++) {
	      s = 0.0;
	      for (k = 0; k < n; k++)
		s += bvp[k*n+i] * wp[k*n+j];
	      vp[i*n+j] += s * ne[kmi];
	    }
	} else
	  for (i = 0; i < n; i++)
	    vp[i*n+i] += bvp[i*n+i] * wp[i*n+i] * ne[kmi];
        // wp <- H Bs == (sum mi(s+miV)^{-1})^{-1} (s+miV)^{-1} S
	if (fullV) {
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++) {
	      s = 0.0;
	      for (k = 0; k < n; k++)
		s += hp[i*n+k] * bsp[k*n+j];
	      wp[i*n+j] = s;
	    }
	} else
	  for (i = 0; i < n; i++)
	    wp[i*n+i] = hp[i*n+i] * bsp[i*n+i];
        // sp <- (sp + Bs^T H Bs) * mi * nens with that mi
	if (fullV) {
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++) {
	      s = 0.0;
	      for (k = 0; k < n; k++)
		s += bsp[k*n+i] * wp[k*n+j];
	      sp[i*n+j] += s * mi * ne[kmi];
	    }
	} else
	  for (i = 0; i < n; i++)
	    sp[i*n+i] += bsp[i*n+i] * wp[i*n+i] * mi * ne[kmi];
      }
    }
 
    for (i = 0; i < n * n; i++) {
      vp[i] /= ngroups;
      sp[i] /= nconfms_tot;
    }
    for (i = 0; i < n; i++)
      for (j = 0; j < i; j++) {
	vp[i*n+j] = vp[j*n+i] = (vp[i*n+j] + vp[j*n+i]) / 2.0;
	sp[i*n+j] = sp[j*n+i] = (sp[i*n+j] + sp[j*n+i]) / 2.0;
      }

    printf("%5d\t", iter);

    s = 0.0;
    for (i = 0; i < n * n; i++) {
      d = fabs(vp[i] - vpp[i]) / (1.0 + fabs(vp[i]));
      if (d > s) s = d;
    }
    printf("%e\t", s);
    s1 = s;
    s = 0.0;
    for (i = 0; i < n * n; i++) {
      d = fabs(sp[i] - spp[i]) / (1.0 + fabs(sp[i]));
      if (d > s) s = d;
    }
    printf("%e\n", s);
    if (s1 < eps && s < eps)
      break;

    // use diagonal S also in the next iteration
    /*
    for (i = 0; i < n; i++)
      for (j = 0; j < i; j++)
	sp[i*n+j] = sp[j*n+i] = 0.0;
    */
  }
  
  fd = open(outfname, O_WRONLY | O_CREAT, 0644);
  if (fd < 0) {
    fprintf(stderr, "cannot open %s.\n", outfname);
    exit(-1);
  }
  write(fd, vp, sizeof(Real) * n * n);
  write(fd, sp, sizeof(Real) * n * n);
  close(fd);

  printf("Normal termination.\n");
  return 0;
}



