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
    
    


#ifndef QMCPLUSPLUS_TEST_SPLINE_H
#define QMCPLUSPLUS_TEST_SPLINE_H
namespace qmcplusplus
{
//template<typename T>
//  inline T std::abs(const TinyVector<T,3>& a, const TitnyVector<T,3>& b)
//  {
//    return std::abs(a[0]-b[0])+std::abs(a[1]-b[1])+std::abs(a[2]-b[2]);
//  }

template<typename T1, typename T2>
bool is_same(int n, const T1* restrict a, const T1* restrict b, T2 eps)
{
  bool yes=true;
  T2 diff=0.0;
  for(int i=0; i<n; i++)
  {
    //T2 x=std::abs((a[i]-b[i])/a[i]);
    T2 x=std::abs(a[i]-b[i]);
    diff=std::max(diff,x);
    if(x>eps)
      yes=false;
    //diff=std::max(diff,std::abs((a[i]-b[i])/a[i]));
    ////if(std::abs(a[i]-b[i])>eps) yes=false;
    //if(std::abs(1-b[i]/a[i])>eps) yes=false;
  }
  for(int i=0; i< std::min(n,9); i++)
    std::cout << i << " " << a[i] << " " << b[i] << " " << a[i]-b[i] << std::endl;
  //cout << "Relative max diff = " << diff << std::endl;
  std::cout << "Absolute max diff = " << diff << std::endl;
  return yes;
}

// template<typename T>
// bool is_same(int n, const TinyVector<T,3>* restrict a, const TinyVector<T,3>* restrict b, T eps)
// {
//   bool yes=true;
//   for(int i=0; i<n; i++)
//   {
//     if(std::abs(1-b[i][0]/a[i][0])+std::abs(1-b[i][1]/a[i][1])+std::abs(1-b[i][2]/a[i][2])>eps) yes=false;
//     //if(std::abs(a[i][0]-b[i][0])+std::abs(a[i][1]-b[i][1])+std::abs(a[i][2]-b[i][2])>eps) yes=false;
//     std::cout << a[i] << " " << a[i]-b[i] << std::endl;
//   }
//   return yes;
// }

//template<typename T>
//bool is_same(int n, const std::complex<T>* restrict a, const std::complex<T>* restrict b, T eps)
//{
//  bool yes=true;
//  T diff_r=0.0, diff_i;
//  for(int i=0; i<n; i++)
//  {
//    diff_r=std::max(diff_r,std::abs(1-b[i].real()/a[i].real()));
//    diff_i=std::max(diff_i,std::abs(1-b[i].real()/a[i].real()));
//    //if(std::abs(a[i]-b[i])>eps) yes=false;
//    if(std::abs(1-b[i].real()/a[i].real())>eps || std::abs(1-b[i].imag()/a[i].imag())>eps)
//    {
//      yes=false;
//    }
//  }

//  for(int i=0; i< std::min(n,8); i++)
//    std::cout << i << " " << a[i] << " " << a[i]-b[i] << std::endl;

//  std::cout << "Absolute max diff = " << diff_r << " " << diff_i << std::endl;
//  return yes;
//}

template<typename SPE1, typename SPE2>
void test_bspline(ParticleSet& TargetPtcl, SPE1& a, SPE2& b)
{
  int N=a.OrbitalSetSize;
  SPOSetBase::RealType eps=static_cast<SPOSetBase::RealType>(numeric_limits<float>::epsilon( ));
  //SPOSetBase::RealType eps=1e-6;
  SPOSetBase::ValueVector_t psi_0(N);
  SPOSetBase::ValueVector_t psi_1(N);
  SPOSetBase::GradVector_t  dpsi_0(N);
  SPOSetBase::GradVector_t  dpsi_1(N);
  SPOSetBase::ValueVector_t d2psi_0(N);
  SPOSetBase::ValueVector_t d2psi_1(N);
  a.evaluate(TargetPtcl,0,psi_0);
  b.evaluate(TargetPtcl,0,psi_1);
  std::cout << "Check values " << std::endl;
  if(is_same(N,psi_0.data(),psi_1.data(),eps))
    std::cout << "Value evaluation Success" << std::endl;
  else
    std::cout << "Value evaluation Failed" << std::endl;
  std::cout << std::endl << "Check VGL " << std::endl;
  a.evaluate(TargetPtcl,0,psi_0,dpsi_0,d2psi_0);
  b.evaluate(TargetPtcl,0,psi_1,dpsi_1,d2psi_1);
  if(is_same(N,psi_0.data(),psi_1.data(),eps))
    std::cout << "VGL Value evaluation Success" << std::endl;
  else
    std::cout << "VGL Value evaluation Failed" << std::endl;
  if(is_same(N*3,&(dpsi_0[0][0]),&(dpsi_1[0][0]),eps))
    std::cout << "VGL Grad evaluation Success" << std::endl;
  else
    std::cout << "VGL Grad evaluation Failed" << std::endl;
  if(is_same(N,d2psi_0.data(),d2psi_1.data(),eps))
    std::cout << "VGL Lap evaluation Success" << std::endl;
  else
    std::cout << "VGL Lap evaluation Failed" << std::endl;
  SPOSetBase::ValueMatrix_t psiM_0(N,N);
  SPOSetBase::ValueMatrix_t psiM_1(N,N);
  SPOSetBase::GradMatrix_t  dpsiM_0(N,N);
  SPOSetBase::GradMatrix_t  dpsiM_1(N,N);
  SPOSetBase::ValueMatrix_t d2psiM_0(N,N);
  SPOSetBase::ValueMatrix_t d2psiM_1(N,N);
  std::cout << std::endl << " evaluate_notranspose " << std::endl;
  a.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  b.evaluate_notranspose(TargetPtcl,0,N,psiM_1,dpsiM_1,d2psiM_1);
  if(is_same(N*N,psiM_0.data(),psiM_1.data(),eps))
    std::cout << "Psi Success " << std::endl;
  else
    std::cout << "Psi Failed!!! " << std::endl;
  //if(is_same(N*N,dpsiM_0.data(),dpsiM_1.data(),eps))
  //  std::cout << "dPsi Success " << std::endl;
  //else
  //  std::cout << "dPsi Failed!!! " << std::endl;
  if(is_same(N*N,d2psiM_0.data(),d2psiM_1.data(),eps))
    std::cout << "d2Psi Success " << std::endl;
  else
    std::cout << "d2Psi Failed!!! " << std::endl;
  Timer t;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      a.evaluate(TargetPtcl,j,psi_0);
  std::cout << "ELAPSED VALUE DOUBLE = " << t.elapsed() << std::endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      b.evaluate(TargetPtcl,j,psi_0);
  std::cout << "ELAPSED VALUE NEW = " << t.elapsed() << std::endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      a.evaluate(TargetPtcl,j,psi_0,dpsi_0,d2psi_0);
  std::cout << "ELAPSED VGL DOUBLE = " << t.elapsed() << std::endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      b.evaluate(TargetPtcl,j,psi_1,dpsi_1,d2psi_1);
  std::cout << "ELAPSED VGL NEW = " << t.elapsed() << std::endl;
  t.restart();
  for(int l=0; l<100; ++l)
    a.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  std::cout << "ELAPSED NOTRANSPOSE = " << t.elapsed() << std::endl;
  t.restart();
  for(int l=0; l<100; ++l)
    b.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  std::cout << "ELAPSED NOTRANSPOSE NEW = " << t.elapsed() << std::endl;
}
}
#endif
