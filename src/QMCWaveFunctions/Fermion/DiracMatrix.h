//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_MATRIX_H
#define QMCPLUSPLUS_DIRAC_MATRIX_H

#include "Numerics/Blasf.h"
#include <OhmmsPETE/OhmmsMatrix.h>
#include "Numerics/BlasThreadingEnv.h"
#include <type_traits/scalar_traits.h>
#include "simd/simd.hpp"

namespace qmcplusplus
{
inline void Xgetrf(int n, int m, float* restrict a, int lda, int* restrict piv)
{
  int status;
  sgetrf(n, m, a, lda, piv, status);
}

inline void Xgetri(int n, float* restrict a, int lda, int* restrict piv, float* restrict work, int& lwork)
{
  int status;
  sgetri(n, a, lda, piv, work, lwork, status);
}

inline void Xgetrf(int n, int m, std::complex<float>* restrict a, int lda, int* restrict piv)
{
  int status;
  cgetrf(n, m, a, lda, piv, status);
}

/** inversion of a float matrix after lu factorization*/
inline void Xgetri(int n,
                   std::complex<float>* restrict a,
                   int lda,
                   int* restrict piv,
                   std::complex<float>* restrict work,
                   int& lwork)
{
  int status;
  cgetri(n, a, lda, piv, work, lwork, status);
}

inline void Xgetrf(int n, int m, double* restrict a, int lda, int* restrict piv)
{
  int status;
  dgetrf(n, m, a, lda, piv, status);
}

inline void Xgetri(int n, double* restrict a, int lda, int* restrict piv, double* restrict work, int& lwork)
{
  int status;
  dgetri(n, a, lda, piv, work, lwork, status);
}

inline void Xgetrf(int n, int m, std::complex<double>* restrict a, int lda, int* restrict piv)
{
  int status;
  zgetrf(n, m, a, lda, piv, status);
}

/** inversion of a std::complex<double> matrix after lu factorization*/
inline void Xgetri(int n,
                   std::complex<double>* restrict a,
                   int lda,
                   int* restrict piv,
                   std::complex<double>* restrict work,
                   int& lwork)
{
  int status;
  zgetri(n, a, lda, piv, work, lwork, status);
}


template<typename TIN, typename TOUT>
inline void TansposeSquare(const TIN* restrict in, TOUT* restrict out, size_t n, size_t lda)
{
#pragma omp simd
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      out[i * lda + j] = in[i + j * lda];
}


template<typename T, typename T_FP>
inline void computeLogDet(const T* restrict diag, int n, const int* restrict pivot, std::complex<T_FP>& logdet)
{
  logdet = std::complex<T_FP>();
  for (size_t i = 0; i < n; i++)
    logdet += std::log(std::complex<T_FP>((pivot[i] == i + 1) ? diag[i] : -diag[i]));
}

/** helper class to compute matrix inversion and the log value of determinant
 * @tparam T_FP the datatype used in the actual computation of matrix inversion
 */
template<typename T_FP>
class DiracMatrix
{
  typedef typename scalar_traits<T_FP>::real_type real_type_fp;
  aligned_vector<T_FP> m_work;
  aligned_vector<int> m_pivot;
  int Lwork;
  /// scratch space used for mixed precision
  Matrix<T_FP> psiM_fp;
  /// LU diagonal elements
  aligned_vector<T_FP> LU_diag;

  /// reset internal work space
  inline void reset(T_FP* invMat_ptr, const int lda)
  {
    m_pivot.resize(lda);
    Lwork = -1;
    T_FP tmp;
    real_type_fp lw;
    Xgetri(lda, invMat_ptr, lda, m_pivot.data(), &tmp, Lwork);
    convert(tmp, lw);
    Lwork = static_cast<int>(lw);
    m_work.resize(Lwork);
    LU_diag.resize(lda);
  }

  /** compute the inverse of invMat (in place) and the log value of determinant
   * @tparam TREAL real type
   * @param n invMat is n x n matrix
   * @param lda the first dimension of invMat
   * @param LogDet log determinant value of invMat before inversion
   */
  template<typename TREAL>
  inline void computeInvertAndLog(T_FP* invMat, const int n, const int lda, std::complex<TREAL>& LogDet)
  {
    BlasThreadingEnv knob(getNextLevelNumThreads());
    if (Lwork < lda)
      reset(invMat, lda);
    Xgetrf(n, n, invMat, lda, m_pivot.data());
    for(int i=0; i<n; i++)
      LU_diag[i] = invMat[i*lda+i];
    computeLogDet(LU_diag.data(), n, m_pivot.data(), LogDet);
    Xgetri(n, invMat, lda, m_pivot.data(), m_work.data(), Lwork);
  }

public:

  DiracMatrix() : Lwork(0) {}

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are the same
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL>
  inline std::enable_if_t<std::is_same<T_FP, TMAT>::value>
  invert_transpose(const Matrix<TMAT>& amat,
                   Matrix<TMAT>& invMat,
                   std::complex<TREAL>& LogDet)
  {
    const int n   = invMat.rows();
    const int lda = invMat.cols();
    simd::transpose(amat.data(), n, amat.cols(), invMat.data(), n, lda);
    computeInvertAndLog(invMat.data(), n, lda, LogDet);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are not the same and need scratch space psiM_fp
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL>
  inline std::enable_if_t<!std::is_same<T_FP, TMAT>::value>
  invert_transpose(const Matrix<TMAT>& amat,
                   Matrix<TMAT>& invMat,
                   std::complex<TREAL>& LogDet)
  {
    const int n   = invMat.rows();
    const int lda = invMat.cols();
    psiM_fp.resize(n,lda);
    simd::transpose(amat.data(), n, amat.cols(), psiM_fp.data(), n, lda);
    computeInvertAndLog(psiM_fp.data(), n, lda, LogDet);
    invMat = psiM_fp;
  }

};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DIRAC_MATRIX_H
