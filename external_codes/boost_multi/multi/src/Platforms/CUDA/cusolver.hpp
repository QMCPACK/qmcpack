//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUSOLVER_H
#define QMCPLUSPLUS_CUSOLVER_H

#include <cusolverDn.h>
#include <complex>
#include <iostream>
#include <string>
#include <stdexcept>

#define cusolverErrorCheck(ans, cause)                \
  {                                                   \
    cusolverAssert((ans), cause, __FILE__, __LINE__); \
  }
/// prints cusolver error messages. Always use cusolverErrorCheck macro.
inline void cusolverAssert(cusolverStatus_t code,
                           const std::string& cause,
                           const char* file,
                           int line,
                           bool abort = true)
{
  if (code != CUSOLVER_STATUS_SUCCESS)
  {
    std::string cusolver_error;
    switch (code)
    {
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      cusolver_error = "CUSOLVER_STATUS_NOT_INITIALIZED";
      break;
    case CUSOLVER_STATUS_ALLOC_FAILED:
      cusolver_error = "CUSOLVER_STATUS_ALLOC_FAILED";
      break;
    case CUSOLVER_STATUS_INVALID_VALUE:
      cusolver_error = "CUSOLVER_STATUS_INVALID_VALUE";
      break;
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      cusolver_error = "CUSOLVER_STATUS_ARCH_MISMATCH";
      break;
    case CUSOLVER_STATUS_EXECUTION_FAILED:
      cusolver_error = "CUSOLVER_STATUS_EXECUTION_FAILED";
      break;
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      cusolver_error = "CUSOLVER_STATUS_INTERNAL_ERROR";
      break;
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      cusolver_error = "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
      break;
    default:
      cusolver_error = "<unknown>";
    }

    std::ostringstream err;
    err << "cusolverAssert: " << cusolver_error << ", file " << file << " , line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    //if (abort) exit(code);
    throw std::runtime_error(cause);
  }
}

namespace qmcplusplus
{
/** interface to cusolver calls for different data types S/C/D/Z
 */
namespace cusolver
{
inline cusolverStatus_t getrf_bufferSize(cusolverDnHandle_t& handle, int m, int n, double* A, int lda, int* lwork)
{
  return cusolverDnDgetrf_bufferSize(handle, m, n, A, lda, lwork);
}

inline cusolverStatus_t getrf_bufferSize(cusolverDnHandle_t& handle,
                                         int m,
                                         int n,
                                         std::complex<double>* A,
                                         int lda,
                                         int* lwork)
{
  return cusolverDnZgetrf_bufferSize(handle, m, n, (cuDoubleComplex*)A, lda, lwork);
}

inline cusolverStatus_t getrf(cusolverDnHandle_t& handle,
                              int m,
                              int n,
                              double* A,
                              int lda,
                              double* work,
                              int* ipiv,
                              int* info)
{
  return cusolverDnDgetrf(handle, m, n, A, lda, work, ipiv, info);
}

inline cusolverStatus_t getrf(cusolverDnHandle_t& handle,
                              int m,
                              int n,
                              std::complex<double>* A,
                              int lda,
                              std::complex<double>* work,
                              int* ipiv,
                              int* info)
{
  return cusolverDnZgetrf(handle, m, n, (cuDoubleComplex*)A, lda, (cuDoubleComplex*)work, ipiv, info);
}

inline cusolverStatus_t getrs(cusolverDnHandle_t& handle,
                              const cublasOperation_t& transa,
                              int m,
                              int n,
                              const double* A,
                              int lda,
                              int* ipiv,
                              double* B,
                              int ldb,
                              int* info)
{
  return cusolverDnDgetrs(handle, transa, m, n, A, lda, ipiv, B, ldb, info);
}

inline cusolverStatus_t getrs(cusolverDnHandle_t& handle,
                              const cublasOperation_t& transa,
                              int m,
                              int n,
                              const std::complex<double>* A,
                              int lda,
                              int* ipiv,
                              std::complex<double>* B,
                              int ldb,
                              int* info)
{
  return cusolverDnZgetrs(handle, transa, m, n, (const cuDoubleComplex*)A, lda, ipiv, (cuDoubleComplex*)B, ldb, info);
}
} // namespace cusolver

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_CUSOLVER_H
