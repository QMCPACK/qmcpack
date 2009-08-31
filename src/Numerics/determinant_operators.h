//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DETERMINANT_OPERATORS_FAST_H
#define QMCPLUSPLUS_DETERMINANT_OPERATORS_FAST_H
#include <complex>
#include <algorithm>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <Numerics/Blasf.h>
#include <Numerics/DeterminantOperators.h>
extern "C"
{
  void dger_(const int& m, const int& n, const double& alpha
      , const double* x, const int& incx, const double* y, const int& incy
      , double* a, const int& lda);
}
namespace qmcplusplus
{

  template<typename T>
    inline void det_row_update_ger(T* restrict pinv,  const T* restrict tv, int m, int rowchanged, T c_ratio)
    {
      const T ratio_inv(1.0/c_ratio);
      T temp[m], rcopy[m];
      dgemv('T', m, m, ratio_inv, pinv, m, tv, 1, T(), temp, 1);
      temp[rowchanged]=1.0-ratio_inv;
      memcpy(rcopy,pinv+m*rowchanged,m*sizeof(T));
      dger_(m,m,-1.0,rcopy,1,temp,1,pinv,m);
    }

  template<typename T>
    inline void det_col_update_ger(T* restrict pinv,  const T* restrict tv, int m, int colchanged, T c_ratio)
    {
      const T ratio_inv(1.0/c_ratio);
      T temp[m], rcopy[m];
      dgemv('N', m, m, ratio_inv, pinv, m, tv, 1, T(), temp, 1);
      temp[colchanged]=1.0-ratio_inv;
      BLAS::copy(m,pinv+colchanged,m,rcopy,1);
      dger_(m,m,-1.0,temp,1,rcopy,1,pinv,m);
    }



  inline double const_one(double c)
  {
    return 1.0;
  }

  inline std::complex<double> const_one(std::complex<double>& c)
  {
    return std::complex<double>(1.0,0.0);
  }

  template<typename T>
    inline void det_row_update(T* restrict pinv,  const T* restrict tv, int m, int rowchanged, T c_ratio)
    {
      const char transa = 'T';
      const T zero=T();
      T ratio_inv=1.0/c_ratio;
      double temp[m], rcopy[m];
      dgemv(transa, m, m, ratio_inv, pinv, m, tv, 1, zero, temp, 1);
      int roffset=rowchanged*m;
      temp[rowchanged]=1.0-ratio_inv;
      std::copy(pinv+roffset,pinv+roffset+m,rcopy);
      for(int i=0,ij=0;i<m;++i)
      {
        T t=temp[i];
        for(int j=0; j<m; ++j) pinv[ij++] -= t*rcopy[j];
      }
    }

  template<typename T>
    inline void multidet_row_update(const T* restrict pinv,  const T* restrict tm, const T* ratios
        , T* restrict new_invs, int m, int rowchanged, int howmany)
    {
      const char transa = 'T';
      const T zero=T();
      double temp[m], rcopy[m];
      int roffset=rowchanged*m;
      int m2=m*m;
      std::copy(pinv+roffset,pinv+roffset+m,rcopy);
      for(int r=0; r<howmany; ++r)
      {
        T ratio_inv=1.0/ratios[r];
        dgemv(transa, m, m, ratio_inv, pinv, m, tm+r*m, 1, zero, temp, 1);
        temp[rowchanged]=1.0-ratio_inv;
        for(int i=0,ij=0,rtot=r*m2;i<m;++i)
        {
          T t=temp[i];
          for(int j=0; j<m; ++j) new_invs[rtot++] = pinv[ij++] - t*rcopy[j];
        }
      }
    }

  template<typename T, typename INDARRAY>
    inline void multidet_row_update(const T* restrict pinv,  const T* restrict tm, const T* ratios
        , T* restrict new_invs, int m, int rowchanged, const INDARRAY& ind)
    {
      const char transa = 'T';
      const T zero=T();
      double temp[m], rcopy[m];
      int roffset=rowchanged*m;
      int m2=m*m;
      std::copy(pinv+roffset,pinv+roffset+m,rcopy);
      for(int k=0; k<ind.size(); ++k)
      {
        int r=ind[k];
        T ratio_inv=1.0/ratios[k];
        dgemv(transa, m, m, ratio_inv, pinv, m, tm+r*m, 1, zero, temp, 1);
        temp[rowchanged]=1.0-ratio_inv;
        for(int i=0,ij=0,rtot=r*m2;i<m;++i)
          for(int j=0; j<m; ++j) new_invs[rtot++] = pinv[ij++] - temp[i]*rcopy[j];
      }
    }

  template<typename MAT, typename VV, typename INDARRAY>
    inline void multidet_row_update(const MAT& pinv,  const MAT& tm, const VV& ratios
        , vector<MAT*> new_invs, int m, int rowchanged, const INDARRAY& ind)
    {
      const char transa = 'T';
      typedef typename MAT::value_type value_type;
      const value_type zero=value_type();
      value_type temp[m], rcopy[m];
      int roffset=rowchanged*m;
      std::copy(pinv.data()+roffset,pinv.data()+roffset+m,rcopy);
      for(int k=0; k<ind.size(); ++k)
      {
        value_type ratio_inv=1.0/ratios[k];
        dgemv(transa, m, m, ratio_inv, pinv.data(), m, tm.data()+(ind[k]+m)*m, 1, zero, temp, 1);
        temp[rowchanged]=1.0-ratio_inv;
        value_type* restrict t=new_invs[k]->data();
        for(int i=0,ij=0;i<m;++i)
          for(int j=0; j<m; ++j,++ij) t[ij] = pinv(ij) - temp[i]*rcopy[j];
      }
    }



  template<typename T>
    inline void det_col_update(T* pinv,  const T* tv, int m, int colchanged, T c_ratio)
    {
      T ratio_inv=1.0/c_ratio;
      for(int i=0; i<m; ++i)
      {
        if(i==colchanged) continue;
        T temp=T();
        for(int k=0; k<m; ++k) temp += tv[k]*pinv[k*m+i];
        temp *= -ratio_inv;
        for(int k=0; k<m; ++k) pinv[k*m+i]+=temp*pinv[k*m+colchanged];
      }
      for(int k=0; k<m; ++k) pinv[k*m+colchanged] *= ratio_inv;
    }


  template<typename T>
    inline void getRatiosByRowSubstitution(const T* restrict tm_new, const T* restrict r_replaced, T* restrict ratios
        , int m, int howmany)
    {
      const char transa = 'T';
      const T one=const_one(T());
      const T zero=T();
      dgemv(transa, m, howmany, one, tm_new, m, r_replaced, 1, zero, ratios, 1);
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
      MAT::value_type res0=BLAS::dot(m,refinv.data()+colchanged,m,tcm.data(),1);
      for(int i=0; i<ind.size(); ++i) ratios[i]=res0-pinned+oldrc*tcm[ind[i]+m];
      return res0;
    }

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
