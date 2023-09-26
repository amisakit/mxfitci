#include "atypes.h"
#include "lapack.h"

#ifndef MKLROOT
#define MKL_INT int
void dgesvd_( char* jobu, char* jobvt, MKL_INT* m,
              MKL_INT* n, Real* a, MKL_INT* lda, Real* s,
              Real* u, MKL_INT* ldu, Real* vt, MKL_INT* ldvt,
              Real* work, MKL_INT* lwork, MKL_INT* info ) {}

void dpotrf_( char* uplo, MKL_INT* n, Real* a,
              MKL_INT* lda, MKL_INT* info ) {}

void dpotri_( char* uplo, MKL_INT* n, Real* a,
              MKL_INT* lda, MKL_INT* info ) {}
void dpotrs_( char* uplo, MKL_INT* n, MKL_INT* nrhs,
              Real* a, MKL_INT* lda, Real* b,
              MKL_INT* ldb, MKL_INT* info ) {}
void dgeev_( char* jobvl, char* jobvr, MKL_INT* n, Real* a,
             MKL_INT* lda, Real* wr, Real* wi, Real* vl,
             MKL_INT* ldvl, Real* vr, MKL_INT* ldvr,
             Real* work, MKL_INT* lwork, MKL_INT* info ) {}
void dsyev_( char* jobz, char* uplo, MKL_INT* n, Real* a,
             MKL_INT* lda, Real* w, Real* work, MKL_INT* lwork,
             MKL_INT* info ) {}

void dsyevr_( char* jobz, char* range, char* uplo,
              MKL_INT* n, Real* a, MKL_INT* lda,
              Real* vl, Real* vu, MKL_INT* il,
              MKL_INT* iu, Real* abstol, MKL_INT* m, Real* w,
              Real* z, MKL_INT* ldz, MKL_INT* isuppz, Real* work,
              MKL_INT* lwork, MKL_INT* iwork, MKL_INT* liwork,
              MKL_INT* info ) {}
#endif
