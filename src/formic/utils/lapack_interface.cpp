#include<formic/utils/lapack_interface.h>
#include<formic/utils/declare.h>
#include<formic/utils/mangle.h>

// complex conjugate of a vector
void formic::xconj(const int n, double * const y, const int incy) {
}
void formic::xconj(const int n, std::complex<double> * const y, const int incy) {
  const int n5 = n - n % 5;
  int i = 0;
  if (incy == 1) {
    for ( ; i < n5; i += 5) {
      y[i+0] = std::conj(y[i+0]);
      y[i+1] = std::conj(y[i+1]);
      y[i+2] = std::conj(y[i+2]);
      y[i+3] = std::conj(y[i+3]);
      y[i+4] = std::conj(y[i+4]);
    }
    for ( ; i < n; i += 1)
      y[i+0] = std::conj(y[i+0]);
  } else {
    for ( ; i < n5; i += 5) {
      y[(i+0)*incy] = std::conj(y[(i+0)*incy]);
      y[(i+1)*incy] = std::conj(y[(i+1)*incy]);
      y[(i+2)*incy] = std::conj(y[(i+2)*incy]);
      y[(i+3)*incy] = std::conj(y[(i+3)*incy]);
      y[(i+4)*incy] = std::conj(y[(i+4)*incy]);
    }
    for ( ; i < n; i += 1)
      y[(i+0)*incy] = std::conj(y[(i+0)*incy]);
  }
}

void formic::dscal(int n, double a, double *x, int incx) {
  ::FC_GLOBAL(dscal,DSCAL)(&n, &a, x, &incx);
}
void formic::xscal(int n, double a, double *x, int incx) {
  ::FC_GLOBAL(dscal,DSCAL)(&n, &a, x, &incx);
}
void formic::zscal(int n, std::complex<double> a, std::complex<double> *x, int incx) {
  ::FC_GLOBAL(zscal,ZSCAL)(&n, (xcomplex*)&a, (xcomplex*)x, &incx);
}
void formic::xscal(int n, std::complex<double> a, std::complex<double> *x, int incx) {
  ::FC_GLOBAL(zscal,ZSCAL)(&n, (xcomplex*)&a, (xcomplex*)x, &incx);
}

void formic::dcopy(int n, const double *x, int incx, double *y, int incy) {
  ::FC_GLOBAL(dcopy,DCOPY)(&n, x, &incx, y, &incy);
}
void formic::xcopy(int n, const double *x, int incx, double *y, int incy) {
  ::FC_GLOBAL(dcopy,DCOPY)(&n, x, &incx, y, &incy);
}
void formic::zcopy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zcopy,ZCOPY)(&n, (const xcomplex*)x, &incx, (xcomplex*)y, &incy);
}
void formic::xcopy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zcopy,ZCOPY)(&n, (const xcomplex*)x, &incx, (xcomplex*)y, &incy);
}

void formic::daxpy(int n, double a, const double *x, int incx, double *y, int incy) {
  ::FC_GLOBAL(daxpy,DAXPY)(&n, &a, x, &incx, y, &incy);
}
void formic::xaxpy(int n, double a, const double *x, int incx, double *y, int incy) {
  ::FC_GLOBAL(daxpy,DAXPY)(&n, &a, x, &incx, y, &incy);
}
void formic::zaxpy(int n, std::complex<double> a, const std::complex<double> *x, int incx, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zaxpy,ZAXPY)(&n, (xcomplex*)&a, (const xcomplex*)x, &incx, (xcomplex*)y, &incy);
}
void formic::xaxpy(int n, std::complex<double> a, const std::complex<double> *x, int incx, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zaxpy,ZAXPY)(&n, (xcomplex*)&a, (const xcomplex*)x, &incx, (xcomplex*)y, &incy);
}

double formic::ddot(int n, const double *x, int incx, const double *y, int incy) {
  return ::FC_GLOBAL(ddot,DDOT)(&n, x, &incx, y, &incy);
}
double formic::xdot(int n, const double *x, int incx, const double *y, int incy) {
  return ::FC_GLOBAL(ddot,DDOT)(&n, x, &incx, y, &incy);
}

/* commented to avoid compiler error.
std::complex<double> formic::zdot(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy) {
  #ifdef FORMIC_HAVE_MKL
  std::complex<double> retval;
  ::FC_GLOBAL(zdotu,ZDOTU)((xcomplex*)&retval, &n, (const xcomplex*)x, &incx, (const xcomplex*)y, &incy);
  //ZDOTU((MKL_Complex16 *)&retval, &n, (const MKL_Complex16 *)x, &incx, (const MKL_Complex16 *)y, &incy);
  return retval;
  #else
  return ::FC_GLOBAL(zdotu,ZDOTU)(&n, (const xcomplex*)x, &incx, (const xcomplex*)y, &incy);
  #endif
}
std::complex<double> formic::xdot(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy) {
  return formic::zdot(n, x, incx, y, incy);
}
*/

double formic::ddotc(int n, const double *x, int incx, const double *y, int incy) {
  return ::FC_GLOBAL(ddot,DDOT)(&n, x, &incx, y, &incy);
}
double formic::xdotc(int n, const double *x, int incx, const double *y, int incy) {
  return ::FC_GLOBAL(ddot,DDOT)(&n, x, &incx, y, &incy);
}

/* commented to avoid compiler error.
std::complex<double> formic::zdotc(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy) {
  #ifdef FORMIC_HAVE_MKL
  std::complex<double> retval;
  ::FC_GLOBAL(zdotc,ZDOTC)((xcomplex*)&retval, &n, (const xcomplex*)x, &incx, (const xcomplex*)y, &incy);
  return retval;
  #else
  return ::FC_GLOBAL(zdotc,ZDOTC)(&n, (const xcomplex*)x, &incx, (const xcomplex*)y, &incy);
  #endif
}
std::complex<double> formic::xdotc(int n, const std::complex<double> *x, int incx, const std::complex<double> *y, int incy) {
  return formic::zdotc(n, x, incx, y, incy);
}
*/

void formic::dgemv(char trans, int m, int n,
                   double alpha, const double *a, int lda, const double *x, int incx,
                   double beta, double *y, int incy) {
  ::FC_GLOBAL(dgemv,DGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
void formic::xgemv(char trans, int m, int n,
                   double alpha, const double *a, int lda, const double *x, int incx,
                   double beta, double *y, int incy) {
  ::FC_GLOBAL(dgemv,DGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
void formic::zgemv(char trans, int m, int n,
                   std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *x, int incx,
                   std::complex<double> beta, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zgemv,ZGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
void formic::xgemv(char trans, int m, int n,
                   std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *x, int incx,
                   std::complex<double> beta, std::complex<double> *y, int incy) {
  ::FC_GLOBAL(zgemv,ZGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void formic::dgemm(char transa, char transb, int m, int n, int k,
                   double alpha, const double *a, int lda, const double *b, int ldb,
                   double beta, double *c, int ldc) {
  ::FC_GLOBAL(dgemm,DGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
void formic::xgemm(char transa, char transb, int m, int n, int k,
                   double alpha, const double *a, int lda, const double *b, int ldb,
                   double beta, double *c, int ldc) {
  ::FC_GLOBAL(dgemm,DGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
void formic::zgemm(char transa, char transb, int m, int n, int k,
                   std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *b, int ldb,
                   std::complex<double> beta, std::complex<double> *c, int ldc) {
  ::FC_GLOBAL(zgemm,ZGEMM)(&transa, &transb, &m, &n, &k,
                           (xcomplex*)&alpha, (const xcomplex*)a, &lda, (const xcomplex*)b, &ldb,
                           (xcomplex*)&beta, (xcomplex*)c, &ldc);
}
void formic::xgemm(char transa, char transb, int m, int n, int k,
                   std::complex<double> alpha, const std::complex<double> *a, int lda, const std::complex<double> *b, int ldb,
                   std::complex<double> beta, std::complex<double> *c, int ldc) {
  ::FC_GLOBAL(zgemm,ZGEMM)(&transa, &transb, &m, &n, &k,
                           (xcomplex*)&alpha, (const xcomplex*)a, &lda, (const xcomplex*)b, &ldb,
                           (xcomplex*)&beta, (xcomplex*)c, &ldc);
}

void formic::dgesv(int n, int nrhs, double *a, int lda,
                   int *ipiv, double *b, int ldb, int & info) {
  ::FC_GLOBAL(dgesv,DGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}
void formic::xgesv(int n, int nrhs, double *a, int lda,
                   int *ipiv, double *b, int ldb, int & info) {
  formic::dgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
}
void formic::zgesv(int n, int nrhs, std::complex<double> *a, int lda,
                   int *ipiv, std::complex<double> *b, int ldb, int & info) {
  ::FC_GLOBAL(zgesv,ZGESV)(&n, &nrhs, (xcomplex*)a, &lda, ipiv, (xcomplex*)b, &ldb, &info);
}
void formic::xgesv(int n, int nrhs, std::complex<double> *a, int lda,
                   int *ipiv, std::complex<double> *b, int ldb, int & info) {
  formic::zgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
}

void formic::dsyev(char jobz, char uplo, int n, double *a, int lda,
                   double * w, double * work, int lwork, int & info) {
  ::FC_GLOBAL(dsyev,DSYEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
}
void formic::xsyev(char jobz, char uplo, int n, double *a, int lda,
                   double * w, double * work, int lwork, double * rwork, int & info) {
  formic::dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
}
void formic::zheev(char jobz, char uplo, int n, std::complex<double> *a, int lda,
                   double * w, std::complex<double> * work, int lwork, double * rwork, int & info) {
  ::FC_GLOBAL(zheev,ZHEEV)(&jobz, &uplo, &n, (xcomplex*)a, &lda, w, (xcomplex*)work, &lwork, rwork, &info);
}
void formic::xsyev(char jobz, char uplo, int n, std::complex<double> *a, int lda,
                   double * w, std::complex<double> * work, int lwork, double * rwork, int & info) {
  formic::zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

void formic::dsygv(int itype, char jobz, char uplo, int n, double *a, int lda,
                   double *b, int ldb, double * w, double * work, int lwork, int & info) {
  ::FC_GLOBAL(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
}
void formic::xsygv(int itype, char jobz, char uplo, int n, double *a, int lda,
                   double *b, int ldb, double * w, double * work, int lwork,
                   double * rwork, int & info) {
  formic::dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}
void formic::zhegv(int itype, char jobz, char uplo, int n, std::complex<double> *a, int lda,
                   std::complex<double> *b, int ldb, double * w,
                   std::complex<double> * work, int lwork, double * rwork, int & info) {
  ::FC_GLOBAL(zhegv,ZHEGV)(&itype, &jobz, &uplo, &n, (xcomplex*)a, &lda, (xcomplex*)b, &ldb,
                           w, (xcomplex*)work, &lwork, rwork, &info);
}
void formic::xsygv(int itype, char jobz, char uplo, int n, std::complex<double> *a, int lda,
                   std::complex<double> *b, int ldb, double * w,
                   std::complex<double> * work, int lwork, double * rwork, int & info) {
  formic::zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}

void formic::dgeev(char jobvl, char jobvr, int n, double *a, int lda,
                   double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr,
                   double *work, int lwork, int & info) {
  ::FC_GLOBAL(dgeev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
}
void formic::xgeev(char jobvl, char jobvr, int n, double *a, int lda,
                   double *w, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr,
                   double *work, int lwork, double *rwork, int & info) {
  formic::dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
void formic::zgeev(char jobvl, char jobvr, int n, std::complex<double> *a, int lda,
                   std::complex<double> *w, std::complex<double> *vl, int ldvl, std::complex<double> *vr, int ldvr,
                   std::complex<double> *work, int lwork, double *rwork, int & info) {
  ::FC_GLOBAL(zgeev,ZGEEV)(&jobvl, &jobvr, &n, (xcomplex*)a, &lda, (xcomplex*)w, (xcomplex*)vl, &ldvl, (xcomplex*)vr, &ldvr,
                           (xcomplex*)work, &lwork, rwork, &info);
}
void formic::xgeev(char jobvl, char jobvr, int n, std::complex<double> *a, int lda,
                   std::complex<double> *w, double *wr, double *wi, std::complex<double> *vl, int ldvl,
                   std::complex<double> *vr, int ldvr, std::complex<double> *work, int lwork,
                   double *rwork, int & info) {
    formic::zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
    for (int i = 0; i < n; i++) {
      wr[i] = w[i].real();
      wi[i] = w[i].imag();
    }
}

void formic::dgesvd(char jobu, char jobvt, int m, int n, double *a, int lda,
                    double *s, double *u, int ldu, double *vt, int ldvt,
                    double *work, int lwork, int & info) {
  ::FC_GLOBAL(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
}
void formic::xgesvd(char jobu, char jobvt, int m, int n, double *a, int lda,
                    double *s, double *u, int ldu, double *vt, int ldvt,
                    double *work, int lwork, double *rwork, int & info) {
  ::FC_GLOBAL(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
}
void formic::zgesvd(char jobu, char jobvt, int m, int n, std::complex<double> *a, int lda,
                    double *s, std::complex<double> *u, int ldu, std::complex<double> *vt, int ldvt,
                    std::complex<double> *work, int lwork, double *rwork, int & info) {
  ::FC_GLOBAL(zgesvd,ZGESVD)(&jobu, &jobvt, &m, &n, (xcomplex*)a, &lda, s, (xcomplex*)u, &ldu, (xcomplex*)vt, &ldvt,
                             (xcomplex*)work, &lwork, rwork, &info);
}
void formic::xgesvd(char jobu, char jobvt, int m, int n, std::complex<double> *a, int lda,
                    double *s, std::complex<double> *u, int ldu, std::complex<double> *vt, int ldvt,
                    std::complex<double> *work, int lwork, double *rwork, int & info) {
  ::FC_GLOBAL(zgesvd,ZGESVD)(&jobu, &jobvt, &m, &n, (xcomplex*)a, &lda, s, (xcomplex*)u, &ldu, (xcomplex*)vt, &ldvt,
                             (xcomplex*)work, &lwork, rwork, &info);
}

void formic::dgetrf(int m, int n, double *a, int lda, int * ipiv, int & info) {
  ::FC_GLOBAL(dgetrf,DGETRF)(&m, &n, a, &lda, ipiv, &info);
}
void formic::xgetrf(int m, int n, double *a, int lda, int * ipiv, int & info) {
  ::FC_GLOBAL(dgetrf,DGETRF)(&m, &n, a, &lda, ipiv, &info);
}
void formic::zgetrf(int m, int n, std::complex<double> *a, int lda, int * ipiv, int & info) {
  ::FC_GLOBAL(zgetrf,ZGETRF)(&m, &n, (xcomplex*)a, &lda, ipiv, &info);
}
void formic::xgetrf(int m, int n, std::complex<double> *a, int lda, int * ipiv, int & info) {
  ::FC_GLOBAL(zgetrf,ZGETRF)(&m, &n, (xcomplex*)a, &lda, ipiv, &info);
}

void formic::dggev(char jobvl, char jobvr, int n, double *a, int lda, double *b, int ldb,
                   double *alphar, double *alphai, double *beta, double *vl, int ldvl, double *vr, int ldvr,
                   double *work, int lwork, int & info) {
  ::FC_GLOBAL(dggev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                           work, &lwork, &info);
}

