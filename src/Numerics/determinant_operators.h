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
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
//extern "C"
//{
//  void dger_(const int& m, const int& n, const double& alpha
//      , const double* x, const int& incx, const double* y, const int& incy
//      , double* a, const int& lda);
//}
namespace qmcplusplus
{
template<typename T> struct const_traits {};

template<>
struct const_traits<double>
{
  typedef double value_type;
  inline static double zero()
  {
    return 0.0;
  }
  inline static double one()
  {
    return 1.0;
  }
  inline static double minus_one()
  {
    return -1.0;
  }
};

template<>
struct const_traits<std::complex<double> >
{
  typedef std::complex<double> value_type;
  inline static std::complex<double> zero()
  {
    return value_type();
  }
  inline static std::complex<double> one()
  {
    return value_type(1.0,0.0);
  }
  inline static std::complex<double> minus_one()
  {
    return value_type(-1.0,0.0);
  }
};

template<>
struct const_traits<float>
{
  typedef float value_type;
  inline static float zero()
  {
    return 0.0f;
  }
  inline static float one()
  {
    return 1.0f;
  }
  inline static float minus_one()
  {
    return -1.0f;
  }
};

template<>
struct const_traits<std::complex<float> >
{
  typedef std::complex<float> value_type;
  inline static std::complex<float> zero()
  {
    return value_type();
  }
  inline static std::complex<float> one()
  {
    return value_type(1.0f,0.0f);
  }
  inline static std::complex<float> minus_one()
  {
    return value_type(-1.0f,0.0f);
  }
};

//template<typename T>
//  inline void det_row_update(T* restrict pinv,  const T* restrict tv, int m, int rowchanged, T c_ratio)
//  {
//    T ratio_inv=1.0/c_ratio;
//    double temp[m], rcopy[m];
//    BLAS::gemv('T', m, m, ratio_inv, pinv, m, tv, 1, T(), temp, 1);
//    int roffset=rowchanged*m;
//    temp[rowchanged]=1.0-ratio_inv;
//    copy(pinv+roffset,pinv+roffset+m,rcopy);
//    for(int i=0,ij=0;i<m;++i)
//    {
//      T t=temp[i];
//      for(int j=0; j<m; ++j) pinv[ij++] -= t*rcopy[j];
//    }
//  }

template<typename T>
inline void det_row_update(T* restrict pinv, const T* restrict tv
                           , int m, int rowchanged, T c_ratio
                           , T* restrict temp, T* restrict rcopy)//pass buffer
{
  //const T ratio_inv(1.0/c_ratio);
  c_ratio=T(1)/c_ratio;
  BLAS::gemv('T', m, m, c_ratio, pinv, m, tv, 1, const_traits<T>::zero(), temp, 1);
  temp[rowchanged]=const_traits<T>::one()-c_ratio;
  memcpy(rcopy,pinv+m*rowchanged,m*sizeof(T));
  BLAS::ger(m,m,const_traits<T>::minus_one(),rcopy,1,temp,1,pinv,m);
}

/** experimental: identical to det_row_update above but using temporary arrays
 */
template<typename T>
inline void det_row_update(T* restrict pinv,  const T* restrict tv, int m, int rowchanged, T c_ratio)
{
  T temp[m], rcopy[m];
  c_ratio=1.0/c_ratio;
  BLAS::gemv('T', m, m, c_ratio, pinv, m, tv, 1, const_traits<T>::zero(), temp, 1);
  temp[rowchanged]=const_traits<T>::one()-c_ratio;
  memcpy(rcopy,pinv+m*rowchanged,m*sizeof(T));
  BLAS::ger(m,m,const_traits<T>::minus_one(),rcopy,1,temp,1,pinv,m);
}


template<typename T>
inline void det_row_block_update
(const T* restrict inv0,  const T* restrict tv, const T* restrict tv_replaced
 , T* restrict inv
 ,int n, int m, int rowchanged, T inv_utv)
{
  T temp[n];
  BLAS::gemv('T', m, n, inv_utv, inv0, m, tv, 1, T(), temp, 1);
  temp[rowchanged]=1.0-inv_utv;
  BLAS::ger(m,n,-1.0,tv_replaced,1,temp,1,inv,m);
}

//template<typename T>
//  inline void det_col_update(T* pinv,  const T* tv, int m, int colchanged, T c_ratio)
//  {
//    T ratio_inv=1.0/c_ratio;
//    for(int i=0; i<m; ++i)
//    {
//      if(i==colchanged) continue;
//      T temp=T();
//      for(int k=0; k<m; ++k) temp += tv[k]*pinv[k*m+i];
//      temp *= -ratio_inv;
//      for(int k=0; k<m; ++k) pinv[k*m+i]+=temp*pinv[k*m+colchanged];
//    }
//    for(int k=0; k<m; ++k) pinv[k*m+colchanged] *= ratio_inv;
//  }
template<typename T>
inline void det_col_update(T* restrict pinv,  const T* restrict tv, int m, int colchanged, T c_ratio
                           , T* restrict temp, T* restrict rcopy)
{
  const T cone(1);
  c_ratio=cone/c_ratio;
  BLAS::gemv('N', m, m, c_ratio, pinv, m, tv, 1, T(), temp, 1);
  temp[colchanged]=cone-c_ratio;
  BLAS::copy(m,pinv+colchanged,m,rcopy,1);
  BLAS::ger(m,m,-1.0,temp,1,rcopy,1,pinv,m);
}

template<typename T>
inline void det_col_update(T* restrict pinv,  const T* restrict tv, int m, int colchanged, T c_ratio)
{
  T temp[m], rcopy[m];
  det_col_update(pinv,tv,m,colchanged,c_ratio,temp,rcopy);
}

template<typename T>
inline void multidet_row_update(const T* restrict pinv,  const T* restrict tm, const T* ratios
                                , T* restrict new_invs, int m, int rowchanged, int howmany)
{
  double temp[m], rcopy[m];
  int roffset=rowchanged*m;
  int m2=m*m;
  copy(pinv+roffset,pinv+roffset+m,rcopy);
  for(int r=0; r<howmany; ++r)
  {
    T ratio_inv=1.0/ratios[r];
    BLAS::gemv('T', m, m, ratio_inv, pinv, m, tm+r*m, 1, T(), temp, 1);
    temp[rowchanged]=1.0-ratio_inv;
    for(int i=0,ij=0,rtot=r*m2; i<m; ++i)
    {
      T t=temp[i];
      for(int j=0; j<m; ++j)
        new_invs[rtot++] = pinv[ij++] - t*rcopy[j];
    }
  }
}

template<typename T, typename INDARRAY>
inline void multidet_row_update(const T* restrict pinv,  const T* restrict tm, const T* ratios
                                , T* restrict new_invs, int m, int rowchanged, const INDARRAY& ind)
{
  double temp[m], rcopy[m];
  int roffset=rowchanged*m;
  int m2=m*m;
  copy(pinv+roffset,pinv+roffset+m,rcopy);
  for(int k=0; k<ind.size(); ++k)
  {
    int r=ind[k];
    T ratio_inv=1.0/ratios[k];
    BLAS::gemv('T', m, m, ratio_inv, pinv, m, tm+r*m, 1, T(), temp, 1);
    temp[rowchanged]=1.0-ratio_inv;
    for(int i=0,ij=0,rtot=r*m2; i<m; ++i)
      for(int j=0; j<m; ++j)
        new_invs[rtot++] = pinv[ij++] - temp[i]*rcopy[j];
  }
}

template<typename MAT, typename VV, typename INDARRAY>
inline void multidet_row_update(const MAT& pinv,  const MAT& tm, const VV& ratios
                                , std::vector<MAT*> new_invs, int m, int rowchanged, const INDARRAY& ind)
{
  typedef typename MAT::value_type value_type;
  value_type temp[m], rcopy[m];
  int roffset=rowchanged*m;
  copy(pinv.data()+roffset,pinv.data()+roffset+m,rcopy);
  for(int k=0; k<ind.size(); ++k)
  {
    value_type ratio_inv=1.0/ratios[k];
    //dgemv(transa, m, m, ratio_inv, pinv.data(), m, tm.data()+(ind[k]+m)*m, 1, zero, temp, 1);
    BLAS::gemv('T', m, m, ratio_inv, pinv.data(), m, tm.data()+(ind[k]+m)*m, 1, value_type(), temp, 1);
    temp[rowchanged]=1.0-ratio_inv;
    value_type* restrict t=new_invs[k]->data();
    for(int i=0,ij=0; i<m; ++i)
      for(int j=0; j<m; ++j,++ij)
        t[ij] = pinv(ij) - temp[i]*rcopy[j];
  }
}


template<typename T>
inline void getRatiosByRowSubstitution(const T* restrict tm_new, const T* restrict r_replaced, T* restrict ratios
                                       , int m, int howmany)
{
  BLAS::gemv('T', m, howmany, const_one(T()), tm_new, m, r_replaced, 1, T(), ratios, 1);
}

template<typename T, typename INDARRAY>
inline void getRatiosByRowSubstitution(const T* restrict tm_new, const T* restrict r_replaced, T* restrict ratios
                                       , int m, const INDARRAY& ind)
{
  for(int i=0; i<ind.size(); ++i)
    ratios[i]=BLAS::dot(r_replaced,tm_new+ind[i]*m,m);
}

template<typename T>
inline void getRatiosByRowSubstitution_dummy(const T* restrict tm_new, const T* restrict r_replaced, T* restrict ratios
    , int m, int howmany)
{
  for(int i=0; i<howmany; ++i)
    ratios[i]=BLAS::dot(r_replaced,tm_new+i*m,m);
}

/** evaluate the determinant ratio with a column substitution
*/
template<typename T>
inline
T getRatioByColSubstitution(const T* restrict pinv, const T* restrict tc, int m, int colchanged)
{
  return BLAS::dot(m,pinv+colchanged,m,tc,1);
}

template<typename MAT, typename VV>
inline typename MAT::value_type
getRatioByColSubstitution(const MAT& pinv, const VV& tc, int colchanged)
{
  return BLAS::dot(pinv.cols(),pinv.data()+colchanged,pinv.cols(),tc.data(),1);
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
inline typename MAT::value_type
getRatioByColSubstitution(const MAT& refinv,const VV& tcm, VV& ratios
                          ,int m , int colchanged, int r_replaced, IV& ind)
{
  typedef typename MAT::value_type value_type;
  //save the value of refinv(r,c)
  value_type old_v=refinv(r_replaced,colchanged);
  value_type pinned=tcm[r_replaced]*old_v;
  typename MAT::value_type res0=BLAS::dot(m,refinv.data()+colchanged,m,tcm.data(),1);
  for(int i=0; i<ind.size(); ++i)
    ratios[i]=res0-pinned+old_v*tcm[ind[i]+m];
  return res0;
}

}
#endif

