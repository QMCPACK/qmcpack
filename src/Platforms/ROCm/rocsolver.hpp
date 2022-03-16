//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ROCSOLVER_H
#define QMCPLUSPLUS_ROCSOLVER_H

// Interface to rocSOLVER linear algebra library.
// File copied and modified from CUDA/cusolver.hpp

#include <rocsolver.h>
#include <complex>
#include <iostream>
#include <string>
#include <stdexcept>

#define rocsolverErrorCheck(ans, cause)                \
  {                                                    \
    rocsolverAssert((ans), cause, __FILE__, __LINE__); \
  }
/// prints rocsolver error messages. Always use rocsolverErrorCheck macro.
inline void rocsolverAssert(rocblas_status code,
                            const std::string& cause,
                            const char* file,
                            int line,
                            bool abort = true)
{
  if (code != rocblas_status_success)
  {
    std::string rocsolver_error;
    switch (code)
    {
    case rocblas_status_invalid_handle:
      rocsolver_error = "rocblas_status_invalid_handle";
      break;
    case rocblas_status_not_implemented:
      rocsolver_error = "rocblas_status_not_implemented";
      break;
    case rocblas_status_invalid_pointer:
      rocsolver_error = "rocblas_status_invalid_pointer";
      break;
    case rocblas_status_invalid_size:
      rocsolver_error = "rocblas_status_invalid_size";
      break;
    case rocblas_status_memory_error:
      rocsolver_error = "rocblas_status_memory_error";
      break;
    case rocblas_status_internal_error:
      rocsolver_error = "rocblas_status_internal_error";
      break;
    case rocblas_status_perf_degraded:
      rocsolver_error = "rocblas_status_perf_degraded";
      break;
    case rocblas_status_size_query_mismatch:
      rocsolver_error = "rocblas_status_size_query_mismatch";
      break;
    case rocblas_status_size_increased:
      rocsolver_error = "rocblas_status_size_increased";
      break;
    case rocblas_status_size_unchanged:
      rocsolver_error = "rocblas_status_size_unchanged";
      break;
    case rocblas_status_invalid_value:
      rocsolver_error = "rocblas_status_invalid_value";
      break;
    case rocblas_status_continue:
      rocsolver_error = "rocblas_status_continue";
      break;
    case rocblas_status_check_numerics_fail:
      rocsolver_error = "rocblas_status_check_numerics_fail";
      break;
    default:
      rocsolver_error = "<unknown>";
    }

    std::ostringstream err;
    err << "rocsolverAssert: " << rocsolver_error << ", file " << file << " , line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    //if (abort) exit(code);
    throw std::runtime_error(cause);
  }
}

namespace qmcplusplus
{
/** interface to rocsolver calls for different data types S/C/D/Z
 */
namespace rocsolver
{


inline rocblas_status getrf(rocblas_handle& handle, int m, int n, double* A, int lda, int* ipiv, int* info)
{
  return rocsolver_dgetrf(handle, m, n, A, lda, ipiv, info);
}

inline rocblas_status getrf(rocblas_handle& handle,
                            int m,
                            int n,
                            std::complex<double>* A,
                            int lda,
                            int* ipiv,
                            int* info)
{
  return rocsolver_zgetrf(handle, m, n, (rocblas_double_complex*)A, lda, ipiv, info);
}

inline rocblas_status getrs(rocblas_handle& handle,
                            const rocblas_operation& transa,
                            int m,
                            int n,
                            double* A,
                            int lda,
                            int* ipiv,
                            double* B,
                            int ldb)
{
  return rocsolver_dgetrs(handle, transa, m, n, A, lda, ipiv, B, ldb);
}

inline rocblas_status getrs(rocblas_handle& handle,
                            const rocblas_operation& transa,
                            int m,
                            int n,
                            std::complex<double>* A,
                            int lda,
                            int* ipiv,
                            std::complex<double>* B,
                            int ldb)
{
  return rocsolver_zgetrs(handle, transa, m, n, (rocblas_double_complex*)A, lda, ipiv, (rocblas_double_complex*)B, ldb);
}

inline rocblas_status getri(rocblas_handle& handle, int n, double* A, int lda, int* ipiv, int* info)
{
  return rocsolver_dgetri(handle, n, A, lda, ipiv, info);
}

inline rocblas_status getri(rocblas_handle& handle, int n, std::complex<double>* A, int lda, int* ipiv, int* info)
{
  return rocsolver_zgetri(handle, n, (rocblas_double_complex*)A, lda, ipiv, info);
}
} // namespace rocsolver

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_ROCSOLVER_H
