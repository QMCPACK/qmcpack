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
/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
using namespace qmcplusplus;

inline double const_one(double c)
{
  return 1.0;
}

inline std::complex<double> const_one(std::complex<double>& c)
{
  return complex<double>(1.0,0.0);
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
    for(int j=0; j<m; ++j) pinv[ij++] -= temp[i]*rcopy[j];
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
      for(int j=0; j<m; ++j) new_invs[rtot++] = pinv[ij++] - temp[i]*rcopy[j];
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
inline void det_col_update(T* pinv,  T* tv, int m, int colchanged, T c_ratio)
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



int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo Welcome(argc,argv,OHMMS::Controller->rank());

  int M=64;
  int nc=16;
  int ic=1;
  int unocc=-1;
  int niters=10000;
  while(ic<argc)
  {
    string c(argv[ic]);
    if(c=="-v")//number of valence states
      M=atoi(argv[++ic]);
    else if(c== "-c")//number of conduction states
      nc=atoi(argv[++ic]);
    else if(c=="-u")//valence state from the top, starting zero
      unocc=M-atoi(argv[++ic])-1;
    else if(c=="-i")//number of iterations
      niters=atoi(argv[++ic]);
    ++ic;
  }
  if(unocc<0) unocc=M-1;
  int bigM=M+nc;
  Matrix<double>  psi_0(M,M),psi_1(M,M),psi_2(M,M),psi_saved(M,M),Identity(M,M);
  Matrix<double>  psi_big(bigM,M);

  Vector<double> psi_v(M);
  for(int i=0; i<psi_big.size();++i) psi_big(i)=Random();
  for(int i=0; i<M;++i) psi_v(i)=Random();
  std::copy(psi_big[0],psi_big[M],psi_0.data());

  psi_saved=psi_0;

  {
    psi_1=psi_0;

    double det_0=invert_matrix(psi_0,true);
    Vector<double> newcols(bigM);
    for(int i=0; i<bigM;++i) newcols(i)=Random();

    //replace a column
    //double rc=DetRatioTranspose(psi_0,newcols.begin(),0);
    //double rc=getRatioByColSubstitution(psi_0.data(),newcols.data(),M,0);
    double rc=getRatioByColSubstitution(psi_0.data(),newcols.data(),M,0);

    
    //replace this
    for(int i=0; i<M; ++i) psi_1(0,i)=newcols(i);
    double det_1=invert_matrix(psi_1,true);

    cout << "Checking column substitution= " << rc
      << " " << det_1/det_0 << endl;
  }

  //save the original matrices
  psi_1=psi_saved;
  double phase=0;
  double logdet_0=InvertWithLog(psi_0.data(),M,M,phase);

  cout.setf(std::ios::scientific, std::ios::floatfield);
  MatrixOperators::product(psi_0,psi_saved,Identity);
  for(int i=0; i<M; ++i) Identity(i,i)-=1.0;
  cout << "Checking identity " << BLAS::norm2(Identity.size(),Identity.data())<< endl;

  Vector<double> ratios(nc), ratios_0(nc);
  vector<int> imap(nc);
  for(int i=0; i<nc; ++i) imap[i]=i;

  Timer myclock;
  for(int iter=0;iter<niters; ++iter)
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios_0.data(),M,imap);
  double t_naive=myclock.elapsed();

  myclock.restart();
  for(int iter=0;iter<niters; ++iter)
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios.data(),M,nc);

  double t_gemv=myclock.elapsed();
  cout << "Timing of ratios gemv= " << t_gemv << " naive=" << t_naive << " better=" << t_naive/t_gemv<< endl;
  ratios_0 -= ratios;
  cout << "Error in the ratios= " <<  BLAS::norm2(ratios_0.size(),ratios_0.data())<< endl;

  //Matrix<double> newinv(nc,M*M);
  vector<Matrix<double>*> newinv_v(nc);
  for(int i=0; i<nc; ++i) newinv_v[i]=new Matrix<double>(M,M);

  double t_multi=0.0, t_dir=0.0;
  for(int iter=0; iter<niters;++iter)
  {
    for(int i=0; i<psi_big.size();++i) psi_big(i)=Random();

    myclock.restart();
    getRatiosByRowSubstitution(psi_big[M],psi_0[unocc],ratios.data(),M,imap);
    //multidet_row_update(psi_0.data(),psi_big[M],ratios.data(),newinv.data(),M,unocc,nc);
    //multidet_row_update(psi_0.data(),psi_big[M],ratios.data(),newinv.data(),M,unocc,imap);
    multidet_row_update(psi_0,psi_big,ratios,newinv_v,M,unocc,imap);
    t_multi+=myclock.elapsed();

    for(int k=0; k<imap.size(); ++k)
    {
      int i=imap[k]+M;
      //std::copy(newinv[k],newinv[k+1],psi_2.data());
      std::copy(newinv_v[k]->begin(),newinv_v[k]->end(),psi_2.data());

      psi_1=psi_saved;
      myclock.restart();
      for(int j=0; j<M; ++j) psi_1(j,unocc)=psi_big(i,j);
      double newphase;
      double logdet_1=InvertWithLog(psi_1.data(),M,M,newphase);
      t_dir+=myclock.elapsed();
      if(iter==4)
      {
        psi_2-=psi_1;
        cout << "Error in Inverse matrix = " << BLAS::norm2(psi_2.size(),psi_2.data()) << endl;
        cout << "Inverse matrix norm2 = " <<BLAS::norm2(psi_1.size(),psi_1.data()) << endl; 
        if(phase==newphase)//too lazy
        cout << "DetRatioTranspose = " << ratios[k] << " error=" << ratios[k]-std::exp(logdet_1-logdet_0) << endl;
      }
    }
  }

  cout << M << " " << nc << " " << t_dir/niters << " " << t_multi/niters << " " << t_dir/t_multi
    << endl;

  for(int i=0; i<nc; ++i) delete newinv_v[i];
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
