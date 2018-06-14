#ifndef FORMIC_LAPACK_DECLARE_OTHER_HEADER
#define FORMIC_LAPACK_DECLARE_OTHER_HEADER

#include<complex>

#include<formic/utils/mangle.h>

// Clang issues a warning if the C return type is std::complex<double>
// Use the C return type instead
// #define complex_ret double _Complex

extern "C" {

  // BLAS
  void FC_GLOBAL(dscal,DSCAL)(int *n, double *a, double *x, int *incx);
  void FC_GLOBAL(zscal,ZSCAL)(int *n, std::complex<double> *a, std::complex<double> *x, int *incx);
  void FC_GLOBAL(dcopy,DCOPY)(int *n, const double *x, int *incx, double *y, int *incy);
  void FC_GLOBAL(zcopy,ZCOPY)(int *n, const std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
  void FC_GLOBAL(daxpy,DAXPY)(int *n, double *a, const double *x, int *incx, double *y, int *incy);
  void FC_GLOBAL(zaxpy,ZAXPY)(int *n, std::complex<double> *a, const std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
  double FC_GLOBAL(ddot,DDOT)(int *n, const double *x, int *incx, const double *y, int *incy);
  // complex_ret FC_GLOBAL(zdotu,ZDOTU)(int *n, const std::complex<double> *x, int *incx, const std::complex<double> *y, int *incy);
  // complex_ret FC_GLOBAL(zdotc,ZDOTC)(int *n, const std::complex<double> *x, int *incx, const std::complex<double> *y, int *incy);
  void FC_GLOBAL(dgemv,DGEMV)(const char *trans, const int *m, const int *n,
                              const double *alpha, const double *a, const int *lda,
                              const double *x, const int *incx,
                              const double *beta, double *y, const int *incy);
  void FC_GLOBAL(zgemv,ZGEMV)(const char *trans, const int *m, const int *n,
                              const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                              const std::complex<double> *x, const int *incx,
                              const std::complex<double> *beta, std::complex<double> *y, const int *incy);
  void FC_GLOBAL(dgemm,DGEMM)(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                              const double *alpha, const double *a, const int *lda,
                              const double *b, const int *ldb,
                              const double *beta, double *c, const int *ldc);
  void FC_GLOBAL(zgemm,ZGEMM)(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                              const std::complex<double> *alpha, const std::complex<double> *a, const int *lda,
                              const std::complex<double> *b, const int *ldb,
                              const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

  // LAPACK
  void FC_GLOBAL(dgesv,DGESV)(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv,
                              double *b, const int *ldb, int *info);
  void FC_GLOBAL(zgesv,ZGESV)(const int *n, const int *nrhs, std::complex<double> *a, const int *lda, int *ipiv,
                              std::complex<double> *b, const int *ldb, int *info);
  void FC_GLOBAL(dsyev,DSYEV)(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
                              double *work, int *lwork, int *info);
  void FC_GLOBAL(zheev,ZHEEV)(char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda, double *w,
                              std::complex<double> *work, int *lwork, double *rwork, int *info);
  void FC_GLOBAL(dsygv,DSYGV)(int * itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb,
                              double *w, double *work, int *lwork, int *info);
  void FC_GLOBAL(zhegv,ZHEGV)(int * itype, char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda,
                              std::complex<double> *b, int *ldb, double *w,
                              std::complex<double> *work, int *lwork, double *rwork, int *info);
  void FC_GLOBAL(dgeev,DGEEV)(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi,
                              double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
  void FC_GLOBAL(zgeev,ZGEEV)(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda, std::complex<double> *w,
                              std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr,
                              std::complex<double> *work, int *lwork, double *rwork, int *info);
  void FC_GLOBAL(dgesvd,DGESVD)(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda,
                                double *s, double *u, int *ldu, double *vt, int *ldvt,
                                double *work, int *lwork, int *info);
  void FC_GLOBAL(zgesvd,ZGESVD)(char *jobu, char *jobvt, int *m, int *n, std::complex<double> *a, int *lda,
                                double *s, std::complex<double> *u, int *ldu, std::complex<double> *vt, int *ldvt,
                                std::complex<double> *work, int *lwork, double *rwork, int *info);
  void FC_GLOBAL(dgetrf,DGETRF)(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
  void FC_GLOBAL(zgetrf,ZGETRF)(int *m, int *n, std::complex<double> *a, int *lda, int *ipiv, int *info);
  void FC_GLOBAL(dggev,DGGEV)(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *b, int *ldb,
                              double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr,
                              double *work, int *lwork, int *info);

}

typedef std::complex<double> xcomplex;

#endif
