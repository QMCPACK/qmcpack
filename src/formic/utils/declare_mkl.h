#ifndef FORMIC_LAPACK_DECLARE_MKL_HEADER
#define FORMIC_LAPACK_DECLARE_MKL_HEADER

#include <mkl_types.h>
#define MKL_Complex16 std::complex<double>
typedef std::complex<double> xcomplex;
#include<mkl.h>

#define FORMIC_HAVE_MKL


#endif
