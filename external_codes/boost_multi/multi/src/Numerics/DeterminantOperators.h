//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 * @brief Define determinant operators
 */
#ifndef OHMMS_NUMERIC_DETERMINANT_H
#define OHMMS_NUMERIC_DETERMINANT_H

#include <algorithm>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CPU/BLAS.hpp"
#include "CPU/math.hpp"
#include "CPU/OMPThreadCountProtectorLA.hpp"
#include "CPU/SIMD/inner_product.hpp"
#include "Numerics/determinant_operators.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
/** LU factorization of double */
inline void LUFactorization(int n, int m, double* restrict a, int n0, int* restrict piv)
{
  int status;
  dgetrf(n, m, a, n0, piv, status);
}

/** LU factorization of float */
inline void LUFactorization(int n, int m, float* restrict a, const int& n0, int* restrict piv)
{
  int status;
  sgetrf(n, m, a, n0, piv, status);
}

/** LU factorization of std::complex<double> */
inline void LUFactorization(int n, int m, std::complex<double>* restrict a, int n0, int* restrict piv)
{
  int status;
  zgetrf(n, m, a, n0, piv, status);
}

/** LU factorization of complex<float> */
inline void LUFactorization(int n, int m, std::complex<float>* restrict a, int n0, int* restrict piv)
{
  int status;
  cgetrf(n, m, a, n0, piv, status);
}

/** Inversion of a double matrix after LU factorization*/
inline void InvertLU(int n, double* restrict a, int n0, int* restrict piv, double* restrict work, int n1)
{
  int status;
  dgetri(n, a, n0, piv, work, n1, status);
}

/** Inversion of a float matrix after LU factorization*/
inline void InvertLU(const int& n,
                     float* restrict a,
                     const int& n0,
                     int* restrict piv,
                     float* restrict work,
                     const int& n1)
{
  int status;
  sgetri(n, a, n0, piv, work, n1, status);
}

/** Inversion of a std::complex<double> matrix after LU factorization*/
inline void InvertLU(int n,
                     std::complex<double>* restrict a,
                     int n0,
                     int* restrict piv,
                     std::complex<double>* restrict work,
                     int n1)
{
  int status;
  zgetri(n, a, n0, piv, work, n1, status);
}

/** Inversion of a complex<float> matrix after LU factorization*/
inline void InvertLU(int n,
                     std::complex<float>* restrict a,
                     int n0,
                     int* restrict piv,
                     std::complex<float>* restrict work,
                     int n1)
{
  int status;
  cgetri(n, a, n0, piv, work, n1, status);
}

/** @}*/

/** inverse a matrix
 * @param x starting address of an n-by-m matrix
 * @param n rows
 * @param m cols
 * @param work workspace array
 * @param pivot integer pivot array
 * @return determinant
 */
template<class T>
inline T Invert(T* restrict x, int n, int m, T* restrict work, int* restrict pivot)
{
  T detvalue(1.0);
  LUFactorization(n, m, x, n, pivot);
  for (int i = 0, ip = 1; i < m; i++, ip++)
  {
    if (pivot[i] == ip)
      detvalue *= x[i * m + i];
    else
      detvalue *= -x[i * m + i];
  }
  InvertLU(n, x, n, pivot, work, n);
  return detvalue;
}

/** determinant of a matrix
 * @param x starting address of an n-by-m matrix
 * @param n rows
 * @param m cols
 * @param pivot integer pivot array
 * @return determinant
 */
template<class T>
inline T Determinant(T* restrict x, int n, int m, int* restrict pivot)
{
  T detvalue(1.0);
  LUFactorization(n, m, x, n, pivot);
  for (int i = 0, ip = 1; i < m; i++, ip++)
  {
    if (pivot[i] == ip)
      detvalue *= x[i * m + i];
    else
      detvalue *= -x[i * m + i];
  }
  return detvalue;
}

/** inverse a matrix
 * @param x starting address of an n-by-m matrix
 * @param n rows
 * @param m cols
 * @return determinant
 *
 * Workspaces are handled internally.
 */
template<class T>
inline T Invert(T* restrict x, int n, int m)
{
  std::vector<int> pivot(n);
  std::vector<T> work(n);
  return Invert(x, n, m, work.data(), pivot.data());
}

template<class T, class T1>
inline void InvertWithLog(T* restrict x, int n, int m, T* restrict work, int* restrict pivot, std::complex<T1>& logdet)
{
  LUFactorization(n, m, x, n, pivot);
  logdet = std::complex<T1>();
  for (int i = 0; i < n; i++)
    logdet += std::log(std::complex<T1>((pivot[i] == i + 1) ? x[i * m + i] : -x[i * m + i]));
  InvertLU(n, x, n, pivot, work, n);
}

/** invert a matrix
 * \param M a matrix to be inverted
 * \param getdet bool, if true, calculate the determinant
 * \return the determinant
 */
template<class MatrixA>
inline typename MatrixA::value_type invert_matrix(MatrixA& M, bool getdet = true)
{
  OMPThreadCountProtectorLA protector;
  using value_type = typename MatrixA::value_type;
  const int n      = M.rows();
  std::vector<int> pivot(n);
  std::vector<value_type> work(n);
  LUFactorization(n, n, M.data(), n, pivot.data());
  value_type det0 = 1.0;
  if (getdet)
  // calculate determinant
  {
    int sign = 1;
    for (int i = 0; i < n; ++i)
    {
      if (pivot[i] != i + 1)
        sign *= -1;
      det0 *= M(i, i);
    }
    det0 *= static_cast<value_type>(sign);
  }
  InvertLU(n, M.data(), n, pivot.data(), work.data(), n);
  return det0;
}

/** determinant ratio with a row substitution
 * @param Minv inverse matrix
 * @param newv row vector
 * @param rowchanged row index to be replaced
 * @return \f$ M^{new}/M\f$
 */
template<typename MatA, typename VecB>
inline typename MatA::value_type DetRatioByRow(const MatA& Minv, const VecB& newv, int rowchanged)
{
  return simd::dot(Minv[rowchanged], newv.data(), Minv.cols());
  //return BLAS::dot(Minv.cols(),Minv[rowchanged],newv.data());
}

/** determinant ratio with a column substitution
 * @param Minv inverse matrix
 * @param newv column vector
 * @param colchanged column index to be replaced
 * @return \f$ M^{new}/M\f$
 */
template<typename MatA, typename VecB>
inline typename MatA::value_type DetRatioByColumn(const MatA& Minv, const VecB& newv, int colchanged)
{
  //use BLAS dot since the stride is not uniform
  return simd::dot(Minv.cols(), Minv.data() + colchanged, Minv.cols(), newv.data(), 1);
}

/** update a inverse matrix by a row substitution
 * @param Minv in/out inverse matrix
 * @param newrow row vector
 * @param rvec workspace
 * @param rvecinv workspace
 * @param rowchanged row index to be replaced
 * @param c_ratio determinant-ratio with the row replacement
 */
template<typename T, typename ALLOC>
inline void InverseUpdateByRow(Matrix<T, ALLOC>& Minv,
                               Vector<T, ALLOC>& newrow,
                               Vector<T, ALLOC>& rvec,
                               Vector<T, ALLOC>& rvecinv,
                               int rowchanged,
                               T c_ratio)
{
  //using gemv+ger
  det_row_update(Minv.data(), newrow.data(), Minv.cols(), rowchanged, c_ratio, rvec.data(), rvecinv.data());
  //int ncols=Minv.cols();
  //typename MatA::value_type ratio_inv=1.0/c_ratio;
  //for(int j=0; j<ncols; j++) {
  //  if(j == rowchanged) continue;
  //  typename MatA::value_type temp = 0.0;
  //  for(int k=0; k<ncols; k++) temp += newrow[k]*Minv(j,k);
  //  temp *= -ratio_inv;
  //  for(int k=0; k<ncols; k++) Minv(j,k) += temp*Minv(rowchanged,k);
  //}
  //for(int k=0; k<ncols; k++) Minv(rowchanged,k) *= ratio_inv;
}

template<typename T, typename ALLOC>
inline void InverseUpdateByColumn(Matrix<T, ALLOC>& Minv,
                                  Vector<T, ALLOC>& newcol,
                                  Vector<T, ALLOC>& rvec,
                                  Vector<T, ALLOC>& rvecinv,
                                  int colchanged,
                                  T c_ratio)
{
  det_col_update(Minv.data(), newcol.data(), Minv.rows(), colchanged, c_ratio, rvec.data(), rvecinv.data());
  //int nrows=Minv.rows();
  //typename MatA::value_type ratio_inv=1.0/c_ratio;
  //for(int i=0; i<nrows; i++) {
  //  if(i == colchanged) continue;
  //  typename MatA::value_type temp = 0.0;
  //  for(int k=0; k<nrows; k++) temp += newcol[k]*Minv(k,i);
  //  temp *= -ratio_inv;
  //  for(int k=0; k<nrows; k++) Minv(k,i) += temp*Minv(k,colchanged);
  //}
  //for(int k=0; k<nrows; k++) Minv(k,colchanged) *= ratio_inv;
}
} // namespace qmcplusplus
#endif
