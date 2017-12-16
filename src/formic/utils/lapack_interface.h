#ifndef FORMIC_UTILS_LAPACK_INTERFACE_HEADER
#define FORMIC_UTILS_LAPACK_INTERFACE_HEADER

#include<complex>

namespace formic {

  void xconj(int n, double *x, int incx);
  void xconj(int n, std::complex<double> *x, int incx);

  void dscal(int n, double a, double *x, int incx);
  void xscal(int n, double a, double *x, int incx);
  void zscal(int n, std::complex<double> a, std::complex<double> *x, int incx);
  void xscal(int n, std::complex<double> a, std::complex<double> *x, int incx);

  void dcopy(int n, const double *x, int incx, double *y, int incy);
  void xcopy(int n, const double *x, int incx, double *y, int incy);
  void zcopy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy);
  void xcopy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy);

  void daxpy(int n, double a, const double *x, int incx, double *y, int incy);
  void xaxpy(int n, double a, const double *x, int incx, double *y, int incy);
  void zaxpy(int n, std::complex<double> a, const std::complex<double> *x, int incx, std::complex<double> *y, int incy);
  void xaxpy(int n, std::complex<double> a, const std::complex<double> *x, int incx, std::complex<double> *y, int incy);

  double ddot(int n, const double *x, int incx, const double *y, int incy);
  double xdot(int n, const double *x, int incx, const double *y, int incy);
  // std::complex<double> zdot(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy);
  // std::complex<double> xdot(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy);

  double ddotc(int n, const double *x, int incx, const double *y, int incy);
  double xdotc(int n, const double *x, int incx, const double *y, int incy);
  // std::complex<double> zdotc(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy);
  std::complex<double> xdotc(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy);

  void dgemv(char trans, int m, int n,
             double alpha, const double *a, int lda, const double *x, int incx,
             double beta, double *y, int incy);
  void xgemv(char trans, int m, int n,
             double alpha, const double *a, int lda, const double *x, int incx,
             double beta, double *y, int incy);
  void zgemv(char trans, int m, int n,
             std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *x, int incx,
             std::complex<double> beta, std::complex<double> *y, int incy);
  void xgemv(char trans, int m, int n,
             std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *x, int incx,
             std::complex<double> beta, std::complex<double> *y, int incy);

  void dgemm(char transa, char transb, int m, int n, int k,
             double alpha, const double *a, int lda, const double *b, int ldb,
             double beta, double *c, int ldc);
  void xgemm(char transa, char transb, int m, int n, int k,
             double alpha, const double *a, int lda, const double *b, int ldb,
             double beta, double *c, int ldc);
  void zgemm(char transa, char transb, int m, int n, int k,
             std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *b, int ldb,
             std::complex<double> beta, std::complex<double> *c, int ldc);
  void xgemm(char transa, char transb, int m, int n, int k,
             std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *b, int ldb,
             std::complex<double> beta, std::complex<double> *c, int ldc);

  void dgesv(int n, int nrhs, double *a, int lda,
                    int *ipiv, double *b, int ldb, int & info);
  void xgesv(int n, int nrhs, double *a, int lda,
                    int *ipiv, double *b, int ldb, int & info);
  void zgesv(int n, int nrhs, std::complex<double> *a, int lda,
                    int *ipiv, std::complex<double> *b, int ldb, int & info);
  void xgesv(int n, int nrhs, std::complex<double> *a, int lda,
                    int *ipiv, std::complex<double> *b, int ldb, int & info);

  void dsyev(char jobz, char uplo, int n, double *a, int lda,
             double * w, double * work, int lwork, int & info);
  void xsyev(char jobz, char uplo, int n, double *a, int lda,
             double * w, double * work, int lwork, double * rwork, int & info);
  void zheev(char jobz, char uplo, int n, std::complex<double> *a, int lda,
             double * w, std::complex<double> * work, int lwork, double * rwork, int & info);
  void xsyev(char jobz, char uplo, int n, std::complex<double> *a, int lda,
             double * w, std::complex<double> * work, int lwork, double * rwork, int & info);

  void dsygv(int itype, char jobz, char uplo, int n, double *a, int lda,
             double *b, int ldb, double * w, double * work, int lwork, int & info);
  void xsygv(int itype, char jobz, char uplo, int n, double *a, int lda,
             double *b, int ldb, double * w, double * work, int lwork,
             double * rwork, int & info);
  void zhegv(int itype, char jobz, char uplo, int n, std::complex<double> *a, int lda,
             std::complex<double> *b, int ldb, double * w,
             std::complex<double> * work, int lwork, double * rwork, int & info);
  void xsygv(int itype, char jobz, char uplo, int n, std::complex<double> *a, int lda,
             std::complex<double> *b, int ldb, double * w,
             std::complex<double> * work, int lwork, double * rwork, int & info);

  void dgeev(char jobvl, char jobvr, int n, double *a, int lda,
             double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr,
             double *work, int lwork, int & info);
  void xgeev(char jobvl, char jobvr, int n, double *a, int lda,
             double *w, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr,
             double *work, int lwork, double *rwork, int & info);
  void zgeev(char jobvl, char jobvr, int n, std::complex<double> *a, int lda,
             std::complex<double> *w, std::complex<double> *vl, int ldvl, std::complex<double> *vr, int ldvr,
             std::complex<double> *work, int lwork, double *rwork, int & info);
  void xgeev(char jobvl, char jobvr, int n, std::complex<double> *a, int lda,
             std::complex<double> *w, double *wr, double *wi, std::complex<double> *vl, int ldvl,
             std::complex<double> *vr, int ldvr, std::complex<double> *work, int lwork,
             double *rwork, int & info);

  void dgesvd(char jobu, char jobvt, int m, int n, double *a, int lda,
              double *s, double *u, int ldu, double *vt, int ldvt,
              double *work, int lwork, int & info);
  void xgesvd(char jobu, char jobvt, int m, int n, double *a, int lda,
              double *s, double *u, int ldu, double *vt, int ldvt,
              double *work, int lwork, double *rwork, int & info);
  void zgesvd(char jobu, char jobvt, int m, int n, std::complex<double> *a, int lda,
              double *s, std::complex<double> *u, int ldu, std::complex<double> *vt, int ldvt,
              std::complex<double> *work, int lwork, double *rwork, int & info);
  void xgesvd(char jobu, char jobvt, int m, int n, std::complex<double> *a, int lda,
              double *s, std::complex<double> *u, int ldu, std::complex<double> *vt, int ldvt,
              std::complex<double> *work, int lwork, double *rwork, int & info);

  void dgetrf(int m, int n, double *a, int lda, int * ipiv, int & info);
  void xgetrf(int m, int n, double *a, int lda, int * ipiv, int & info);
  void zgetrf(int m, int n, std::complex<double> *a, int lda, int * ipiv, int & info);
  void xgetrf(int m, int n, std::complex<double> *a, int lda, int * ipiv, int & info);

  void dggev(char jobvl, char jobvr, int n, double *a, int lda, double *b, int ldb,
             double *alphar, double *alphai, double *beta, double *vl, int ldvl, double *vr, int ldvr,
             double *work, int lwork, int & info);

}

#endif
