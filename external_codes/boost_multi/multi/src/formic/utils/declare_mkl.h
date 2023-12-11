#ifndef FORMIC_LAPACK_DECLARE_MKL_HEADER
#define FORMIC_LAPACK_DECLARE_MKL_HEADER

#include <mkl_types.h>
#define MKL_Complex16 std::complex<double>
using xcomplex = std::complex<double>;
#include<mkl.h>

#define FORMIC_HAVE_MKL


#endif
