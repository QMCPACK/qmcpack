#ifndef QMCPLUSPLUS_TEST_SPLINE_H
#define QMCPLUSPLUS_TEST_SPLINE_H
namespace qmcplusplus
{
//template<typename T>
//  inline T abs(const TinyVector<T,3>& a, const TitnyVector<T,3>& b)
//  {
//    return abs(a[0]-b[0])+abs(a[1]-b[1])+abs(a[2]-b[2]);
//  }

template<typename T1, typename T2>
bool is_same(int n, const T1* restrict a, const T1* restrict b, T2 eps)
{
  bool yes=true;
  T2 diff=0.0;
  for(int i=0; i<n; i++)
  {
    //T2 x=abs((a[i]-b[i])/a[i]);
    T2 x=abs(a[i]-b[i]);
    diff=std::max(diff,x);
    if(x>eps)
      yes=false;
    //diff=std::max(diff,abs((a[i]-b[i])/a[i]));
    ////if(abs(a[i]-b[i])>eps) yes=false;
    //if(abs(1-b[i]/a[i])>eps) yes=false;
  }
  for(int i=0; i< std::min(n,9); i++)
    cout << i << " " << a[i] << " " << b[i] << " " << a[i]-b[i] << endl;
  //cout << "Relative max diff = " << diff << endl;
  cout << "Absolute max diff = " << diff << endl;
  return yes;
}

// template<typename T>
// bool is_same(int n, const TinyVector<T,3>* restrict a, const TinyVector<T,3>* restrict b, T eps)
// {
//   bool yes=true;
//   for(int i=0; i<n; i++)
//   {
//     if(abs(1-b[i][0]/a[i][0])+abs(1-b[i][1]/a[i][1])+abs(1-b[i][2]/a[i][2])>eps) yes=false;
//     //if(abs(a[i][0]-b[i][0])+abs(a[i][1]-b[i][1])+abs(a[i][2]-b[i][2])>eps) yes=false;
//     cout << a[i] << " " << a[i]-b[i] << endl;
//   }
//   return yes;
// }

//template<typename T>
//bool is_same(int n, const complex<T>* restrict a, const complex<T>* restrict b, T eps)
//{
//  bool yes=true;
//  T diff_r=0.0, diff_i;
//  for(int i=0; i<n; i++)
//  {
//    diff_r=std::max(diff_r,abs(1-b[i].real()/a[i].real()));
//    diff_i=std::max(diff_i,abs(1-b[i].real()/a[i].real()));
//    //if(abs(a[i]-b[i])>eps) yes=false;
//    if(abs(1-b[i].real()/a[i].real())>eps || abs(1-b[i].imag()/a[i].imag())>eps)
//    {
//      yes=false;
//    }
//  }

//  for(int i=0; i< std::min(n,8); i++)
//    cout << i << " " << a[i] << " " << a[i]-b[i] << endl;

//  cout << "Absolute max diff = " << diff_r << " " << diff_i << endl;
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
  cout << "Check values " << endl;
  if(is_same(N,psi_0.data(),psi_1.data(),eps))
    cout << "Value evaluation Success" << endl;
  else
    cout << "Value evaluation Failed" << endl;
  cout << endl << "Check VGL " << endl;
  a.evaluate(TargetPtcl,0,psi_0,dpsi_0,d2psi_0);
  b.evaluate(TargetPtcl,0,psi_1,dpsi_1,d2psi_1);
  if(is_same(N,psi_0.data(),psi_1.data(),eps))
    cout << "VGL Value evaluation Success" << endl;
  else
    cout << "VGL Value evaluation Failed" << endl;
  if(is_same(N*3,&(dpsi_0[0][0]),&(dpsi_1[0][0]),eps))
    cout << "VGL Grad evaluation Success" << endl;
  else
    cout << "VGL Grad evaluation Failed" << endl;
  if(is_same(N,d2psi_0.data(),d2psi_1.data(),eps))
    cout << "VGL Lap evaluation Success" << endl;
  else
    cout << "VGL Lap evaluation Failed" << endl;
  SPOSetBase::ValueMatrix_t psiM_0(N,N);
  SPOSetBase::ValueMatrix_t psiM_1(N,N);
  SPOSetBase::GradMatrix_t  dpsiM_0(N,N);
  SPOSetBase::GradMatrix_t  dpsiM_1(N,N);
  SPOSetBase::ValueMatrix_t d2psiM_0(N,N);
  SPOSetBase::ValueMatrix_t d2psiM_1(N,N);
  cout << endl << " evaluate_notranspose " << endl;
  a.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  b.evaluate_notranspose(TargetPtcl,0,N,psiM_1,dpsiM_1,d2psiM_1);
  if(is_same(N*N,psiM_0.data(),psiM_1.data(),eps))
    cout << "Psi Success " << endl;
  else
    cout << "Psi Failed!!! " << endl;
  //if(is_same(N*N,dpsiM_0.data(),dpsiM_1.data(),eps))
  //  cout << "dPsi Success " << endl;
  //else
  //  cout << "dPsi Failed!!! " << endl;
  if(is_same(N*N,d2psiM_0.data(),d2psiM_1.data(),eps))
    cout << "d2Psi Success " << endl;
  else
    cout << "d2Psi Failed!!! " << endl;
  Timer t;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      a.evaluate(TargetPtcl,j,psi_0);
  cout << "ELAPSED VALUE DOUBLE = " << t.elapsed() << endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      b.evaluate(TargetPtcl,j,psi_0);
  cout << "ELAPSED VALUE NEW = " << t.elapsed() << endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      a.evaluate(TargetPtcl,j,psi_0,dpsi_0,d2psi_0);
  cout << "ELAPSED VGL DOUBLE = " << t.elapsed() << endl;
  t.restart();
  for(int l=0; l<100; ++l)
    for(int j=0; j<TargetPtcl.getTotalNum(); ++j)
      b.evaluate(TargetPtcl,j,psi_1,dpsi_1,d2psi_1);
  cout << "ELAPSED VGL NEW = " << t.elapsed() << endl;
  t.restart();
  for(int l=0; l<100; ++l)
    a.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  cout << "ELAPSED NOTRANSPOSE = " << t.elapsed() << endl;
  t.restart();
  for(int l=0; l<100; ++l)
    b.evaluate_notranspose(TargetPtcl,0,N,psiM_0,dpsiM_0,d2psiM_0);
  cout << "ELAPSED NOTRANSPOSE NEW = " << t.elapsed() << endl;
}
}
#endif
