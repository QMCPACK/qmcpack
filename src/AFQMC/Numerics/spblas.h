
#ifndef AFQMC_SPBLAS_H
#define AFQMC_SPBLAS_H

#if defined(HAVE_MKL)

extern "C" {

void mkl_scsrmv(const char &transa, const MKL_INT &m, const MKL_INT &k, const float &alpha, const char *matdescra, const float *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const float *x, const float &beta, float *y);
void mkl_ccsrmv(const char &transa, const MKL_INT &m, const MKL_INT &k, const std::complex<float> &alpha, const char *matdescra, const std::complex<float> *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const std::complex<float> *x, const std::complex<float> &beta, std::complex<float> *y);
void mkl_dcsrmv(const char &transa, const MKL_INT &m, const MKL_INT &k, const double &alpha, const char *matdescra, const double *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const double *x, const double &beta, double *y);
void mkl_zcsrmv(const char &transa, const MKL_INT &m, const MKL_INT &k, const std::complex<double> &alpha, const char *matdescra, const std::complex<double> *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const std::complex<double> *x, const std::complex<double> &beta, std::complex<double> *y);

void mkl_scsrmm(const char &transa, const MKL_INT &m, const MKL_INT &n, const MKL_INT &k, const float &alpha, const char *matdescra, const float *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const float *b, const MKL_INT &ldb, const float &beta, float *c, const MKL_INT &ldc);
void mkl_ccsrmm(const char &transa, const MKL_INT &m, const MKL_INT &n, const MKL_INT &k, const std::complex<float> &alpha, const char *matdescra, const std::complex<float> *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const std::complex<float> *b, const MKL_INT &ldb, const std::complex<float> &beta, std::complex<float> *c, const MKL_INT &ldc);
void mkl_dcsrmm(const char &transa, const MKL_INT &m, const MKL_INT &n, const MKL_INT &k, const double &alpha, const char *matdescra, const double *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const double *b, const MKL_INT &ldb, const double &beta, double *c, const MKL_INT &ldc);
void mkl_zcsrmm(const char &transa, const MKL_INT &m, const MKL_INT &n, const MKL_INT &k, const std::complex<double> &alpha, const char *matdescra, const std::complex<double> *val, const MKL_INT *indx, const MKL_INT *pntrb, const MKL_INT *pntre, const std::complex<double> *b, const MKL_INT &ldb, const std::complex<double> &beta, std::complex<double> *c, const MKL_INT &ldc);



void mkl_dcsrmultd (const char *trans , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , double *a , MKL_INT *ja , MKL_INT *ia , double *b , MKL_INT *jb , MKL_INT *ib , double *c , MKL_INT *ldc );

void mkl_scsrmultd (const char *trans , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , float *a , MKL_INT *ja , MKL_INT *ia , float *b , MKL_INT *jb , MKL_INT *ib , float *c , MKL_INT *ldc );

void mkl_ccsrmultd (const char *trans , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , std::complex<float> *a , MKL_INT *ja , MKL_INT *ia , std::complex<float> *b , MKL_INT *jb , MKL_INT *ib , std::complex<float> *c , MKL_INT *ldc );

void mkl_zcsrmultd (const char *trans , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , std::complex<double> *a , MKL_INT *ja , MKL_INT *ia , std::complex<double> *b , MKL_INT *jb , MKL_INT *ib , std::complex<double> *c , MKL_INT *ldc );


}


#endif



#endif
