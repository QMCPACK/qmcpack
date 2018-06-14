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
    
    



/** @file DeterminantOperator.h
 * @brief Define determinant operators
 */
#ifndef OHMMS_NUMERIC_DETERMINANT_H
#define OHMMS_NUMERIC_DETERMINANT_H

#include <algorithm>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <Numerics/OhmmsBlas.h>
#include <config/stdlib/math.h>
#include <simd/simd.hpp>
#include <Numerics/determinant_operators.h>

namespace qmcplusplus
{

/** LU factorization of double */
inline void
LUFactorization(int n, int m, double* restrict a, int n0, int* restrict piv)
{
  int status;
  dgetrf(n,m,a,n0,piv,status);
}

/** LU factorization of float */
inline void
LUFactorization(int n, int m, float* restrict a, const int& n0, int* restrict piv)
{
  int status;
  sgetrf(n,m,a,n0,piv,status);
}

/** LU factorization of std::complex<double> */
inline void
LUFactorization(int n, int m, std::complex<double>* restrict a, int n0, int* restrict piv)
{
  int status;
  zgetrf(n,m,a,n0,piv,status);
}

/** LU factorization of complex<float> */
inline void
LUFactorization(int n, int m, std::complex<float>* restrict a, int n0, int* restrict piv)
{
  int status;
  cgetrf(n,m,a,n0,piv,status);
}

/** Inversion of a double matrix after LU factorization*/
inline void InvertLU(int n, double* restrict a, int n0 , int* restrict piv, double* restrict work, int n1)
{
  int status;
  dgetri(n,a,n0,piv,work,n1,status);
}

/** Inversion of a float matrix after LU factorization*/
inline void InvertLU(const int& n, float* restrict a, const int& n0,
                     int* restrict piv, float* restrict work, const int& n1)
{
  int status;
  sgetri(n,a,n0,piv,work,n1,status);
}

/** Inversion of a std::complex<double> matrix after LU factorization*/
inline void InvertLU(int n, std::complex<double>* restrict a, int n0
                     , int* restrict piv, std::complex<double>* restrict work, int n1)
{
  int status;
  zgetri(n,a,n0,piv,work,n1,status);
}

/** Inversion of a complex<float> matrix after LU factorization*/
inline void InvertLU(int n, std::complex<float>* restrict a, int n0
                     , int* restrict piv, std::complex<float>* restrict work, int n1)
{
  int status;
  cgetri(n,a,n0,piv,work,n1,status);
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
inline T
Invert(T* restrict x, int n, int m, T* restrict work, int* restrict pivot)
{
  T detvalue(1.0);
  LUFactorization(n,m,x,n,pivot);
  for(int i=0,ip=1; i<m; i++, ip++)
  {
    if(pivot[i]==ip)
      detvalue *= x[i*m+i];
    else
      detvalue *= -x[i*m+i];
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
inline T
Determinant(T* restrict x, int n, int m,int* restrict pivot)
{
  T detvalue(1.0);
  LUFactorization(n,m,x,n,pivot);
  for(int i=0,ip=1; i<m; i++, ip++)
  {
    if(pivot[i]==ip)
      detvalue *= x[i*m+i];
    else
      detvalue *= -x[i*m+i];
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
  int pivot[n];
  T work[n];
  return Invert(x,n,m,work,pivot);
}

template<class T>
inline T
InvertWithLog(T* restrict x, int n, int m
              , T* restrict work, int* restrict pivot, T& phase)
{
  T logdet(0.0);
  LUFactorization(n,m,x,n,pivot);
  int sign_det=1;
  for(int i=0; i<n; i++)
  {
    sign_det *= (pivot[i] == i+1)?1:-1;
    sign_det *= (x[i*m+i]>0)?1:-1;
    logdet += std::log(std::abs(x[i*m+i]));
  }
  InvertLU(n, x, n, pivot, work, n);
  phase=(sign_det>0)?0.0:M_PI;
  return logdet;
}

template<class T>
inline T
InvertWithLog(std::complex<T>* restrict x, int n, int m
              , std::complex<T>* restrict work, int* restrict pivot
              , T& phase)
{
  T logdet(0.0);
  LUFactorization(n,m,x,n,pivot);
  phase=0.0;
  for(int i=0; i<n; i++)
  {
    int ii=i*m+i;
    phase += std::arg(x[ii]);
    if(pivot[i]!=i+1)
      phase += M_PI;
    logdet+=std::log(x[ii].real()*x[ii].real()+x[ii].imag()*x[ii].imag());
//slightly smaller error with the following
//        logdet+=2.0*std::log(std::abs(x[ii]);
  }
  InvertLU(n, x, n, pivot, work, n);
  const T one_over_2pi=1.0/TWOPI;
  phase -= std::floor(phase*one_over_2pi)*TWOPI;
  return 0.5*logdet;
}

/** matrix inversion using log method to avoid numerical over/under-flow
 * @param x starting address of an n-by-m
 * @param n rows
 * @param m cols
 * @param phase output
 * @return |det|
 *
 * Use type-specific InvertWithLog defined above.
 */
template<typename T, typename RT>
inline RT
InvertWithLog(T* restrict x, int n, int m, RT& phase)
{
  int pivot[m];
  T work[m];
  return InvertWithLog(x,n,m,work,pivot,phase);
}

/** invert a matrix
 * \param M a matrix to be inverted
 * \param getdet bool, if true, calculate the determinant
 * \return the determinant
 */
template<class MatrixA>
inline typename MatrixA::value_type
invert_matrix(MatrixA& M, bool getdet=true)
{
  typedef typename MatrixA::value_type value_type;
  const int n=M.rows();
  int pivot[n];
  value_type work[n];
  LUFactorization(n,n,M.data(),n,pivot);
  value_type det0 = 1.0;
  if(getdet)
    // calculate determinant
  {
    int sign = 1;
    for(int i=0; i<n; ++i)
    {
      if(pivot[i] != i+1)
        sign *= -1;
      det0 *= M(i,i);
    }
    det0 *= static_cast<value_type>(sign);
  }
  InvertLU(n, M.data(), n, pivot, work, n);
  return det0;
}

/** determinant a matrix
 * \param M a matrix to be inverted
 * \return the determinant
 */
template<class MatrixA>
inline typename MatrixA::value_type
determinant_matrix(MatrixA& M)
{
  typedef typename MatrixA::value_type value_type;
  MatrixA N(M);
  const int n=N.rows();
  int pivot[n];
  value_type work[n];
  LUFactorization(n,n,N.data(),n,pivot);
  value_type det0 = 1.0;
  int sign = 1;
  for(int i=0; i<n; ++i)
  {
    if(pivot[i] != i+1)
      sign *= -1;
    det0 *= M(i,i);
  }
  det0 *= static_cast<value_type>(sign);
  return det0;
}
//
///** invert a matrix
// * \param M a matrix to be inverted
// * \param getdet bool, if true, calculate the determinant
// * \return the determinant
// */
//template<class MatrixA>
//  inline typename MatrixA::value_type
//invert_matrix_log(MatrixA& M, int &sign_det, bool getdet)
//{
//  typedef typename MatrixA::value_type value_type;
//  std::vector<int> pivot(M.rows());
//  std::vector<value_type> work(M.rows());
//  int n(M.rows());
//  int m(M.cols());
//  MatrixA Mcopy(M);
//  LUFactorization(n,m,M.data(),n,&pivot[0]);
//  value_type logdet = 0.0;
//  sign_det=1;
//  if(getdet)
//  {
//    for(int i=0; i<n; i++)
//    {
//      ////if(pivot[i] != i+1) sign_det *= -1;
//      sign_det *= (pivot[i] == i+1)?1:-1;
//      sign_det *= (M(i*m+i)>0)?1:-1;
//      logdet += std::log(std::abs(M(i*m+i)));
//    }
//  }
//  InvertLU(n,M.data(),m, &pivot[0], &work[0], n);
//  //value_type det0 = Invert(Mcopy.data(),n,m,&work[0], &pivot[0]);
//  //double expdetp = sign_det*std::exp(logdet);
//  //std::cerr <<"DETS ARE NOW "<<det0<<" "<<expdetp<<" "<<logdet<< std::endl;
//  return logdet;
//}
//
/** determinant ratio with a row substitution
 * @param Minv inverse matrix
 * @param newv row vector
 * @param rowchanged row index to be replaced
 * @return \f$ M^{new}/M\f$
 */
template<typename MatA, typename VecB>
inline
typename MatA::value_type
DetRatioByRow(const MatA& Minv, const VecB& newv, int rowchanged)
{
  return simd::dot(Minv[rowchanged],newv.data(),Minv.cols());
  //return BLAS::dot(Minv.cols(),Minv[rowchanged],newv.data());
}

/** determinant ratio with a column substitution
 * @param Minv inverse matrix
 * @param newv column vector
 * @param colchanged column index to be replaced
 * @return \f$ M^{new}/M\f$
 */
template<typename MatA, typename VecB>
inline
typename MatA::value_type
DetRatioByColumn(const MatA& Minv, const VecB& newv, int colchanged)
{
  //use BLAS dot since the stride is not uniform
  return simd::dot(Minv.cols(),Minv.data()+colchanged,Minv.cols(),newv.data(),1);
}

/** update a inverse matrix by a row substitution
 * @param Minv in/out inverse matrix
 * @param newrow row vector
 * @param rvec workspace
 * @param rvecinv workspace
 * @param rowchanged row index to be replaced
 * @param c_ratio determinant-ratio with the row replacement
 */
template<class MatA, class VecT>
inline void InverseUpdateByRow(MatA& Minv, VecT& newrow
                               , VecT& rvec, VecT& rvecinv
                               , int rowchanged, typename MatA::value_type c_ratio
                              )
{
  //using gemv+ger
  det_row_update(Minv.data(),newrow.data(),Minv.cols(),rowchanged,c_ratio,rvec.data(),rvecinv.data());
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

template<typename MatA, typename VecT>
inline void InverseUpdateByColumn(MatA& Minv, VecT& newcol
                                  , VecT& rvec, VecT& rvecinv
                                  , int colchanged, typename MatA::value_type c_ratio)
{
  det_col_update(Minv.data(),newcol.data(),Minv.rows(),colchanged,c_ratio
                 ,rvec.data(), rvecinv.data());
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

// template<class T1, class T2>
// inline T2
// dot(const T1* restrict a, const T2* restrict b, int n) {
//   T2 res;
//   for(int i=0; i<n; i++) res += a[i]*b[i];
//   return res;
// }

}
#endif
