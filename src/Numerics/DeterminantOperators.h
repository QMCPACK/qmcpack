//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
/** @file DeterminantOperator.h
 * @brief Define determinant operators 
 */
#ifndef OHMMS_NUMERIC_DETERMINANT_H
#define OHMMS_NUMERIC_DETERMINANT_H

#include <algorithm>
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/Blasf.h"

namespace qmcplusplus {

  /** LU factorization of double */
  inline void 
    LUFactorization(const int& n, const int& m, double* restrict a, const int& n0, 
        int* restrict piv) {
      int status;
      dgetrf(n,m,a,n0,piv,status);
    }

  /** LU factorization of complex<double> */
  inline void 
    LUFactorization(const int& n, const int& m, std::complex<double>* restrict a, 
        const int& n0, int* restrict piv) {
      int status;
      zgetrf(n,m,a,n0,piv,status);
    }

  /** Inversion of a double matrix after LU factorization*/
  inline void InvertLU(const int& n, double* restrict a, const int& n0, 
      int* restrict piv, double* restrict work, const int& n1){
    int status;
    dgetri(n,a,n0,piv,work,n1,status);
  }

  /** Inversion of a complex<double> matrix after LU factorization*/
  inline void InvertLU(const int& n, std::complex<double>* restrict a, const int& n0, 
      int* restrict piv, std::complex<double>* restrict work, const int& n1){
    int status;
    zgetri(n,a,n0,piv,work,n1,status);
  }

  template<class T>
    inline T 
    Invert(T* restrict x, int n, int m, T* restrict work, int* restrict pivot) {
      T detvalue(1.0);
      LUFactorization(n,m,x,n,pivot);
      for(int i=0,ip=1; i<m; i++, ip++) {
        if(pivot[i]==ip) 
          detvalue *= x[i*m+i];
        else 
          detvalue *= -x[i*m+i];
      }
      InvertLU(n, x, n, pivot, work, n);
      return detvalue;
    }

  template<class T>
    inline T Invert(T* restrict x, int n, int m) {
      T detvalue(1.0);
      vector<T> work(n);
      vector<int> pivot(n);
      LUFactorization(n,m,x,n,&pivot[0]);
      for(int i=0,ip=1; i<m; i++, ip++) {
        if(pivot[i]==ip) 
          detvalue *= x[i*m+i];
        else 
          detvalue *= -x[i*m+i];
      }
      InvertLU(n,x, n, &pivot[0], &work[0], n);
      return detvalue;
    }

    inline double
    InvertWithLog(double* restrict x, int n, int m, double* restrict work, int* restrict pivot, 
        double& phase) 
    {
      double logdet(0.0);
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

/** invert a matrix
 * \param M a matrix to be inverted
 * \param getdet bool, if true, calculate the determinant
 * \return the determinant
 */
template<class MatrixA>
inline double
invert_matrix(MatrixA& M, bool getdet=true) {

  typedef typename MatrixA::value_type value_type;
  Vector<int> pivot(M.rows());
  Vector<value_type> work(M.rows());
  int status;

  dgetrf(M.rows(), M.cols(), M.data(), M.rows(), pivot.data(), status);

  value_type det0 = 1.0; 

  if(getdet) {// calculate determinant
    int sign = 1;
    for(int i=0; i<M.rows(); ++i){
      if(pivot[i] != i+1) sign *= -1;
      det0 *= M(i,i);
    }
    det0 *= static_cast<value_type>(sign);
  }

  dgetri(M.rows(), M.data(),  M.rows(), pivot.data(), work.data(), 
      M.rows(), status);
  return det0;
}  

/** invert a matrix                                                                                                                                          
 * \param M a matrix to be inverted                                                                                                                          
 * \param getdet bool, if true, calculate the determinant                                                                                                    
 * \return the determinant                                                                                                                                   
 */
template<class MatrixA>
  inline typename MatrixA::value_type
invert_matrix_log(MatrixA& M, int &sign_det, bool getdet) 
{
  typedef typename MatrixA::value_type value_type;
  vector<int> pivot(M.rows());
  vector<value_type> work(M.rows());
  int n(M.rows());
  int m(M.cols());
  MatrixA Mcopy(M);
  LUFactorization(n,m,M.data(),n,&pivot[0]);
  value_type logdet = 0.0;
  sign_det=1;
  if(getdet) 
  {
    for(int i=0; i<n; i++)
    {
      ////if(pivot[i] != i+1) sign_det *= -1;
      sign_det *= (pivot[i] == i+1)?1:-1;
      sign_det *= (M(i*m+i)>0)?1:-1;
      logdet += std::log(std::abs(M(i*m+i)));
    }
  }
  InvertLU(n,M.data(),m, &pivot[0], &work[0], n);
  //value_type det0 = Invert(Mcopy.data(),n,m,&work[0], &pivot[0]);
  //double expdetp = sign_det*std::exp(logdet);
  //cerr<<"DETS ARE NOW "<<det0<<" "<<expdetp<<" "<<logdet<<endl;
  return logdet;
}

template<typename MatA, typename Iter>
inline
typename MatA::value_type 
DetRatio(const MatA& Minv, Iter newrow, int rowchanged) {
  typename MatA::value_type res = 0.0;
  for(int j=0; j<Minv.cols(); j++, newrow++)
    res += Minv(rowchanged,j)*(*newrow);
  return res;
}

template<typename MatA, typename Iter>
inline
typename MatA::value_type 
DetRatioTranspose(const MatA& Minv, Iter newcol, int colchanged) {
  typename MatA::value_type res = 0.0;
  for(int i=0; i<Minv.rows(); i++, newcol++)
    res += Minv(i,colchanged)*(*newcol);
  return res;
}

template<class T, class Iter2>
inline T
detRatio(T* restrict Minv, Iter2 newrow, int cols) {
  T res(0.0);
  for(int j=0; j<cols; j++) res += (*Minv++)*(*newrow++);
  return res;
}

template<class MatA, class VecT>
inline
void
DetUpdate(MatA& Minv, 
    VecT& newrow, 
    VecT& rvec, 
    VecT& rvecinv, 
    int rowchanged,
    typename MatA::value_type c_ratio) {

  int ncols=Minv.cols();
  typename MatA::value_type ratio_inv=1.0/c_ratio;
  for(int j=0; j<ncols; j++) {
    if(j == rowchanged) continue;
    typename MatA::value_type temp = 0.0;
    for(int k=0; k<ncols; k++) temp += newrow[k]*Minv(j,k);
    temp *= -ratio_inv;
    for(int k=0; k<ncols; k++) Minv(j,k) += temp*Minv(rowchanged,k);
  }
  for(int k=0; k<ncols; k++) Minv(rowchanged,k) *= ratio_inv;
  /*
     for(int j=0; j<ncols; j++) {
     rvec[j] = DetRatio(Minv,newrow.begin(),j);
     rvecinv[j] = 1.0/rvec[j];
     }

     for(int j=0; j<ncols; j++) {
     Minv(rowchanged, j) *= rvecinv[rowchanged];
     }

     for(int i=0; i<rowchanged; i++) {
     for(int j=0; j<ncols; j++) {
     Minv(i,j) -= Minv(rowchanged,j)*rvec[i];
     }
     }
     for(int i=rowchanged+1; i<Minv.rows(); i++) {
     for(int j=0; j<ncols; j++) {
     Minv(i,j) -= Minv(rowchanged,j)*rvec[i];
     }
     }
     */
}

template<class MatA, class VecT>
inline
void
DetUpdateTranspose(MatA& Minv, 
    VecT& newcol, 
    VecT& rvec, 
    VecT& rvecinv, 
    int colchanged,
    typename MatA::value_type c_ratio) {

  int nrows=Minv.rows();
  typename MatA::value_type ratio_inv=1.0/c_ratio;
  for(int i=0; i<nrows; i++) {
    if(i == colchanged) continue;
    typename MatA::value_type temp = 0.0;
    for(int k=0; k<nrows; k++) temp += newcol[k]*Minv(k,i);
    temp *= -ratio_inv;
    for(int k=0; k<nrows; k++) Minv(k,i) += temp*Minv(k,colchanged);
  }
  for(int k=0; k<nrows; k++) Minv(k,colchanged) *= ratio_inv;
}

/** inline dot product 
 * @param a starting address of an array of type T
 * @param b starting address of an array of type TinyVector<T,D>
 * @param n size
 * @return \f$ {\bf v} = \sum_i a[i] {\bf b}[i]\f$
 */
template<class T, unsigned D>
inline TinyVector<T,D>
dot(const T* a, const TinyVector<T,D>* b, int n) {
  TinyVector<T,D> res;
  for(int i=0; i<n; i++) res += a[i]*b[i];
  return res;
}

/** inline dot product: replacement of BLAS::dot
 * @param a starting address of an array of type T
 * @param b starting address of an array of type T
 * @param n size
 * @return  \f$ res = \sum_i a[i] b[i]\f$
 *
 * The arguments of this inline function are different from BLAS::dot
 * This is more efficient than BLAS::dot due to the function overhead,
 * if a compiler knows how to inline.
 */
template<class T>
inline T 
dot(const T* restrict a, const T* restrict b, int n) {
  T res = 0.0;
  for(int i=0; i<n; i++) res += a[i]*b[i];
  return res;
}

template<class T>
inline std::complex<T> 
dot(const std::complex<T>* restrict a, const T* restrict b, int n) {
  std::complex<T> res = 0.0;
  for(int i=0; i<n; i++) res += a[i]*b[i];
  return res;
}

template<class T>
inline std::complex<T> 
dot(const T* restrict a, const std::complex<T>* restrict b, int n) {
  std::complex<T> res = 0.0;
  for(int i=0; i<n; i++) res += a[i]*b[i];
  return res;
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
