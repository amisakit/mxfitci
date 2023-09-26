#ifndef _MMFIT_H
#define _MMFIT_H

extern void abnormalexit(char *message, int code);

/* superimpose.c */
extern int center_M(Real *mp, Real diagS[], int natoms, Real *t);
extern int
center_MZ(Real *mp, Real Zbuf[], Real diagS[], int nconfms[], int ngroups, int nconfms_tot, int natoms, Real *t, Real *work);
extern void
estimate_M_EM(Real *Ybuf, Real *Zbuf, int nconfms[], int ngroups,
	      int natoms, Real *mp);
extern int
estimate_M_ML(Real *Ybuf, Ensemble *Bbuf, int ngroups, int natoms,
	      int njvars, Real *diagV, Real *mp, Real *work);
extern int
align_M(Real *mp, int natoms, Real *vout);
extern void
estimate_Rt_EM(Ensemble Bbuf[], Real *mp, Real diagS[], int ngroups, int natoms, int m_tot, int Mcentered, int YinM, Real *work);
extern void
estimate_Rt_ML(Ensemble Bbuf[], Real *mp, Real *diagV, Real diagS[], int ngroups, int natoms, int Mcentered, Real *work);
extern void
calculate_Y(Real *xp, int *op, Real *rp, Real *t, Real *mp, Real *zp,
            int natoms, Real *yp);

/* map.c */
extern void
estimate_Z(Ensemble Bbuf[], Real *diagV, Real *mp, int natoms, int ngroups,
	   Real *work);

/* mmvars.c */
extern void
estimate_Uaux(struct jvar Jbuf[], Real V[], Real *diagV, Real diagS[], int natoms, Real *work);

extern int
estimate_V(Real *Zbuf, struct jvar Jbuf[], int nconfms[], int ngroups, int natoms, int method, Real *diagV, Real *vp, Real *ani, Real *work);

extern void
estimate_S(Real *Ybuf, Real *mp, Real *Zbuf, struct jvar Jbuf[],
	   int nconfms_tot, int nconfms[], int ngroups, int natoms, int method,
	   Real *diagV, Real diagS[], Real *sp, Real *ani, Real *work);

extern Real
estimate_SV_BML(Real Ybuf[], Real *mp, int natoms, int nconfms_tot,
	       int nconfms[], int ngroups, Real *diagS, Real *vp, Real *work);
Real
estimate_SV_MLold(Ensemble Bbuf[], Real *mp, int natoms,
               int nconfms_tot, int nconfms[], int ngroups,
               Real *diagS, Real *vp, Real *work);

#endif	/* _MMFIT_H */
