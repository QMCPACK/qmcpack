//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DETERMINANT_OPERATORS_FAST_H
#define QMCPLUSPLUS_DETERMINANT_OPERATORS_FAST_H
#include <complex>
#include <algorithm>
#include <cstring>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
//extern "C"
//{
//  void dger_(const int& m, const int& n, const double& alpha
//      , const double* x, const int& incx, const double* y, const int& incy
//      , double* a, const int& lda);
//}
namespace qmcplusplus
{
template<typename T>
struct const_traits
{};

template<>
struct const_traits<double>
{
  using value_type = double;
  inline static double zero() { return 0.0; }
  inline static double one() { return 1.0; }
  inline static double minus_one() { return -1.0; }
};

template<>
struct const_traits<std::complex<double>>
{
  using value_type = std::complex<double>;
  inline static std::complex<double> zero() { return value_type(); }
  inline static std::complex<double> one() { return value_type(1.0, 0.0); }
  inline static std::complex<double> minus_one() { return value_type(-1.0, 0.0); }
};

template<>
struct const_traits<float>
{
  using value_type = float;
  inline static float zero() { return 0.0f; }
  inline static float one() { return 1.0f; }
  inline static float minus_one() { return -1.0f; }
};

template<>
struct const_traits<std::complex<float>>
{
  using value_type = std::complex<float>;
  inline static std::complex<float> zero() { return value_type(); }
  inline static std::complex<float> one() { return value_type(1.0f, 0.0f); }
  inline static std::complex<float> minus_one() { return value_type(-1.0f, 0.0f); }
};

template<typename T>
inline void det_row_update(T* restrict pinv,
                           const T* restrict tv,
                           int m,
                           int rowchanged,
                           T c_ratio,
                           T* restrict temp,
                           T* restrict rcopy) //pass buffer
{
  //const T ratio_inv(1.0/c_ratio);
  c_ratio = T(1) / c_ratio;
  BLAS::gemv('T', m, m, c_ratio, pinv, m, tv, 1, const_traits<T>::zero(), temp, 1);
  temp[rowchanged] = const_traits<T>::one() - c_ratio;
  memcpy(rcopy, pinv + m * rowchanged, m * sizeof(T));
  BLAS::ger(m, m, const_traits<T>::minus_one(), rcopy, 1, temp, 1, pinv, m);
}

template<typename T>
inline void det_col_update(T* restrict pinv,
                           const T* restrict tv,
                           int m,
                           int colchanged,
                           T c_ratio,
                           T* restrict temp,
                           T* restrict rcopy)
{
  const T cone(1);
  c_ratio = cone / c_ratio;
  BLAS::gemv('N', m, m, c_ratio, pinv, m, tv, 1, T(), temp, 1);
  temp[colchanged] = cone - c_ratio;
  BLAS::copy(m, pinv + colchanged, m, rcopy, 1);
  BLAS::ger(m, m, -1.0, temp, 1, rcopy, 1, pinv, m);
}

template<typename T>
inline void getRatiosByRowSubstitution(const T* restrict tm_new,
                                       const T* restrict r_replaced,
                                       T* restrict ratios,
                                       int m,
                                       int howmany)
{
  BLAS::gemv('T', m, howmany, const_one(T()), tm_new, m, r_replaced, 1, T(), ratios, 1);
}

template<typename T, typename INDARRAY>
inline void getRatiosByRowSubstitution(const T* restrict tm_new,
                                       const T* restrict r_replaced,
                                       T* restrict ratios,
                                       int m,
                                       const INDARRAY& ind)
{
  for (int i = 0; i < ind.size(); ++i)
    ratios[i] = BLAS::dot(r_replaced, tm_new + ind[i] * m, m);
}

template<typename T>
inline void getRatiosByRowSubstitution_dummy(const T* restrict tm_new,
                                             const T* restrict r_replaced,
                                             T* restrict ratios,
                                             int m,
                                             int howmany)
{
  for (int i = 0; i < howmany; ++i)
    ratios[i] = BLAS::dot(r_replaced, tm_new + i * m, m);
}

/** evaluate the determinant ratio with a column substitution
*/
template<typename T>
inline T getRatioByColSubstitution(const T* restrict pinv, const T* restrict tc, int m, int colchanged)
{
  return BLAS::dot(m, pinv + colchanged, m, tc, 1);
}

template<typename MAT, typename VV>
inline typename MAT::value_type getRatioByColSubstitution(const MAT& pinv, const VV& tc, int colchanged)
{
  return BLAS::dot(pinv.cols(), pinv.data() + colchanged, pinv.cols(), tc.data(), 1);
}


/** evaluate the ratio with a column substitution and multiple row substitutions
 *
 * @param refinv reference inverse(m,m)
 * @param tcm new column whose size > m
 * @param ratios ratios of the row substitutions with the column substitution
 * @param m dimension
 * @param colchanged ratio with a column substitution of refinv
 * @param r_replaced row index for the excitations in ind
 * @param ind indexes for the excited states (rows)
 * @return the ratio when a column is replaced for the refinv
 */
template<typename MAT, typename VV, typename IV>
inline typename MAT::value_type getRatioByColSubstitution(const MAT& refinv,
                                                          const VV& tcm,
                                                          VV& ratios,
                                                          int m,
                                                          int colchanged,
                                                          int r_replaced,
                                                          IV& ind)
{
  using value_type = typename MAT::value_type;
  //save the value of refinv(r,c)
  value_type old_v              = refinv(r_replaced, colchanged);
  value_type pinned             = tcm[r_replaced] * old_v;
  typename MAT::value_type res0 = BLAS::dot(m, refinv.data() + colchanged, m, tcm.data(), 1);
  for (int i = 0; i < ind.size(); ++i)
    ratios[i] = res0 - pinned + old_v * tcm[ind[i] + m];
  return res0;
}

} // namespace qmcplusplus
#endif
