extern void rand_seed(int);
extern Real rand_uniformc(void);
extern Real rand_uniformo(void);
extern Real rand_stdnorm(void);
extern int rand_mvnorm(Real *mp, Real *cov, int n, Real *xp, int fact);
extern int rand_m3norm(Real *mp, Real *cov, int n, Real *xp, int fact);
extern void rand_rotmat(Real *rp);
