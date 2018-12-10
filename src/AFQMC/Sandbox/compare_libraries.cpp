
#include "AFQMC/config.h"
#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#include<random>

#include "AFQMC/Sandbox/compare_libraries.h"
//#include "AFQMC/Numerics/DenseMatrixOperations.h"
//#include "AFQMC/Numerics/SparseMatrixOperations.h"

#if defined(HAVE_MKL)
#include "mkl.h"
#include "mkl_service.h"
#endif

namespace qmcplusplus
{

typedef std::vector<s2D<ValueType> >::iterator  s2Dit;

void compare_libraries(int NMO, int NAEA, int NAEB, std::vector<s2D<ComplexType> >& Propg_H1, std::vector<IndexType>& Propg_H1_indx, std::vector<s2D<ValueType> >& Vuv, std::vector<s2D<ComplexType> >& vn, std::vector<IndexType>& vn_indx) {

#if defined(HAVE_MKL)
  int ntimes=10;
  std::vector<ComplexType> sigma(vn_indx.size()-1);
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,1.0);
  for (int i=0; i<sigma.size(); ++i)
    sigma[i] = static_cast<ComplexType>(distribution(generator));

  // full matrices  
  ComplexMatrix  P0(2*NMO,NMO), V0(2*NMO*NMO,2*NMO*NMO);
  for(int i=Propg_H1_indx[0]; i<Propg_H1_indx[1]+Propg_H1_indx[0]; i++) {
    IndexType a,b;
    ComplexType v; 
    std::tie(a,b,v) = Propg_H1[i];
    P0(a,b)=v;
  }

/*  
  for(int i=0; i<2*NMO*NMO; i++)
  for(int j=0; j<2*NMO*NMO; j++) V0(i,j)=ComplexType(0.0); 
  for(int i=0; i<Vuv.size(); i++) {
    IndexType a,b;
    ComplexType v;
    std::tie(a,b,v) = Vuv[i];
    V0(a,b)+=v;
  }
*/
    
/*
  for(int i=0; i<Vuv.size(); i++) {
    IndexType ii,kk,jj,ll;
    ii = a/NMO;
    kk = a%NMO;
    jj = b/NMO;
    ll = b%NMO;
    a = (ii*NMO)+ll;
    b = (jj*NMO)+kk;
    V0(a,b)-=v; 
    V2uv.push_back(std::forward_as_tuple(a,b,v));
  } 
*/

 
  ComplexMatrix S0(2*NMO,NAEA), D0(2*NMO,NMO), vHS(2*NMO,NMO);
  for(int i=0; i<NAEA; i++) S0(i,i)=ComplexType(1.0);   
  for(int i=0; i<NAEA; i++) D0(i,i)=ComplexType(1.0);   
  for(int i=0; i<NAEB; i++) S0(i+NMO,i)=ComplexType(1.0);   
  for(int i=0; i<NAEB; i++) D0(i+NMO,i)=ComplexType(1.0);   

  ComplexType one(1.0,0.0); 
  ComplexType zero(0.0,0.0); 
  std::vector<ComplexType> D0_mkl(2*NMO*NMO);
  for(int i=0; i<2*NMO*NMO; i++) D0_mkl[i] = zero; 
  for(int i=0; i<NAEA; i++) D0_mkl[i*NMO+i] = one; 
  for(int i=0; i<NAEB; i++) D0_mkl[(i+NMO)*NMO+i] = one; 

  ComplexMatrix::iterator itG = D0.begin();
  ComplexType epot_manual = 0;
  s2Dit end2 = Vuv.end();
  for(s2Dit it = Vuv.begin(); it != end2; it++) 
    epot_manual += (*(itG + std::get<0>(*it))) * (*(itG + std::get<1>(*it))) * std::get<2>(*it);

  std::cout<<std::endl <<std::endl;
  std::cout<<"******************************************************* \n";
  std::cout<<" TESTING PERFORMANCE OF VARIOUS LINEAR ALGEBRA PACKAGES \n";
  std::cout<<"******************************************************* \n";
  std::cout<<std::endl <<std::endl;

  std::cout<<" OhmmsPETE + mySparse(mkl): \n\n";

  ComplexMatrix vHS_manual(2*NMO,NMO); 
  vHS_manual=ComplexType(0.0);
  for(int i=0; i<vn_indx.size()-1; i++) {
    ComplexType scl = sigma[i];
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++)
      vHS_manual( std::get<0>( vn[n] ) , std::get<1>( vn[n] ) ) -= scl*std::get<2>( vn[n] );
  }

  // OhmmsPETE + mySparseMatrix
  std::string str;
  str="vn_manual"; 
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++) {
      ComplexType scl = sigma[i];
      for(int n = vn_indx[i]; n<vn_indx[i+1]; n++)
        vHS( std::get<0>( vn[n] ) , std::get<1>( vn[n] ) ) -= scl*std::get<2>( vn[n] );
    }  
    Timer.stop(str.c_str()); 
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  std::vector<int> vn_1dindx(vn.size()); 
  std::vector<ComplexType> vn_1ddata(vn.size()); 
  for(int i=0; i<vn_indx.size()-1; i++) 
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) { 
      vn_1dindx[n] = std::get<0>(vn[n])*NMO+std::get<1>(vn[n]);
      vn_1ddata[n] = std::get<2>(vn[n]);
    }

  // switching sigma to -sigma
  for(int i=0; i<vn_indx.size()-1; i++) 
    sigma[i] *= ComplexType(-1.0);   

  ComplexMatrix vHS_zaxpyi(2*NMO,NMO);
  vHS_zaxpyi=ComplexType(0.0);
  for(int i=0; i<vn_indx.size()-1; i++)
    if(vn_indx[i+1] > vn_indx[i])
      cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS_zaxpyi.data()); 


  std::cout<<"Magnitude of difference between vHS_manual and vHS_zaxpyi: \n";
  RealType diff = 0.0;
  for(int i=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++)
     diff += std::abs(vHS_manual(i,j)-vHS_zaxpyi(i,j));
  std::cout<<diff <<std::endl;

  mkl_set_dynamic( 0 );
  str="vn_mkl_cblas_zaxpyi 1 thr";
  mkl_set_num_threads(1); 
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++) 
      if(vn_indx[i+1] > vn_indx[i])
        cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS.data());
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

/*
  str="vn_mkl_cblas_zaxpyi 2 thr";
  mkl_set_num_threads(2);
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++)
      if(vn_indx[i+1] > vn_indx[i])
        cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS.data());
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  str="vn_mkl_cblas_zaxpyi 4 thr";
  mkl_set_num_threads(4);
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++)
      if(vn_indx[i+1] > vn_indx[i])
        cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS.data());
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  str="vn_mkl_cblas_zaxpyi 8 thr";
  mkl_set_num_threads(8);
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++)
      if(vn_indx[i+1] > vn_indx[i])
        cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS.data());
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  str="vn_mkl_cblas_zaxpyi 16 thr";
  mkl_set_num_threads(16);
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    vHS=ComplexType(0.0);
    Timer.start(str.c_str());
    for(int i=0; i<vn_indx.size()-1; i++)
      if(vn_indx[i+1] > vn_indx[i])
        cblas_zaxpyi (vn_indx[i+1]-vn_indx[i], &(sigma[i]), &(vn_1ddata[vn_indx[i]]), &(vn_1dindx[vn_indx[i]]), vHS.data());
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;
*/
/*
  ComplexMatrix vn_full(2*NMO*NMO,vn_indx.size()-1); 
  for(int i=0; i<vn_indx.size()-1; i++) 
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++)
       vn_full( std::get<0>( vn[n] )*NMO+std::get<1>( vn[n] ), i ) = std::get<2>( vn[n] );

  std::vector<int> vn_1drows(2*NMO*NMO);
  std::vector<ComplexType> vn_1ddata2(vn.size()); 
  std::vector<ComplexType> sigma2(vn_indx.size()-1);
  std::vector<ComplexType> vHS2(2*NMO*NMO);
  for(int i=0, cnt=0; i<2*NMO*NMO; i++)  { 
   vn_1drows[i]=cnt; 
   for(int j=0; j<vn_indx.size()-1; j++) { 
     if(std::abs(vn_full(i,j)) > 1e-8) {
      vn_1dindx[cnt] = j; 
      vn_1ddata2[cnt].real = vn_full(i,j).real(); 
      vn_1ddata2[cnt].imag = vn_full(i,j).imag(); 
     }
   }
  }
*/
  ComplexMatrix vn_full(NMO*NMO,vn_indx.size()-1); 
  for(int i=0; i<vn_indx.size()-1; i++) 
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++)
       vn_full( std::get<0>( vn[n] )*NMO+std::get<1>( vn[n] ), i ) = std::get<2>( vn[n] );

  diff=0.0;
  for(int i=0; i<NMO; i++)  { 
   for(int k=0; k<NMO; k++)  { 
    ComplexType t;
    int ik = i*NMO+k;
    for(int j=0; j<vn_indx.size()-1; j++) { 
      t += vn_full(ik,j)*sigma[j];
    } 
    diff += std::abs( vHS_manual(i,k) - t );
   } 
  }
  std::cout<<"Magnitude of difference between vHS_manual and vHS_full: \n";
  std::cout<<diff <<std::endl;  
  

  std::vector<int> vn_1drows(NMO*NMO+1);
  std::vector<ComplexType> vn_1ddata2(vn.size()); 
  std::vector<ComplexType> sigma2(vn_indx.size()-1);
  std::vector<ComplexType> vHS2(2*NMO*NMO);
  int cnt=0;
  for(int i=0; i<NMO*NMO; i++)  { 
   vn_1drows[i]=cnt; 
   for(int j=0; j<vn_indx.size()-1; j++) { 
     if(std::abs(vn_full(i,j)) > 1e-10) {
      vn_1dindx[cnt] = j; 
      vn_1ddata2[cnt++] = vn_full(i,j);
     }
   }
  }
  vn_1drows.back() = vn_1ddata2.size();

  for(int j=0; j<vn_indx.size()-1; j++) 
    sigma2[j] = sigma[j];

  int nrows = vn_full.rows();
  int ncols = vn_full.cols();
  char trans = 'N';

  for(int i=0; i<2*NMO*NMO; i++) vHS2[i] = 0.0; 

  std::vector<int> vn_1drows3(NMO*NMO+1);
  std::vector<int> vn_1dindx3(vn.size());
  std::vector<ComplexType> vn_1ddata3(vn.size());

  for(int i=0; i<vn_indx.size()-1; i++) {
    vn_1drows3[i] = vn_indx[i]; 
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) {
       vn_1ddata3[n] = std::get<2>( vn[n] );
       vn_1dindx3[n] = std::get<0>( vn[n] )*NMO + std::get<1>( vn[n] );
    }
  }
  vn_1drows3.back() = vn_1ddata3.size();  

  mkl_set_num_threads(1);
  char matdes[6];
  matdes[0] = 'G'; 
  matdes[3] = 'C'; 
  
  str="vn_mkl_zcsrmv (transposed) 1 thr";
  trans = 'T';
  mkl_zcsrmv( &trans, &ncols, &nrows, &one, matdes, vn_1ddata3.data() , vn_1dindx3.data(),  vn_1drows3.data() ,  &(vn_1drows3[1]), sigma2.data(), &zero, vHS2.data() );

  std::cout<<"Magnitude of difference between vHS_manual and vHS_zcsrmv (in transposed form): \n";
  diff = 0.0;
  for(int i=0, k=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++,k++) 
     diff += std::abs(vHS_manual(i,j)-vHS2[k]);
  std::cout<<diff <<std::endl;

  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    Timer.start(str.c_str());
    mkl_zcsrmv( &trans, &ncols, &nrows, &one, matdes, vn_1ddata3.data() , vn_1dindx3.data(),  vn_1drows3.data() ,  &(vn_1drows3[1]), sigma2.data(), &zero, vHS2.data() );
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  str="vn_mkl_zcsrmv 1 thr";
  trans = 'N';

  mkl_zcsrmv( &trans, &nrows, &ncols, &one, matdes, vn_1ddata2.data() , vn_1dindx.data(),  vn_1drows.data() ,  &(vn_1drows[1]), sigma2.data(), &zero, vHS2.data() );

  std::cout<<"Magnitude of difference between vHS_manual and vHS_zcsrmv: \n";
  diff = 0.0;
  for(int i=0, k=0; i<2*NMO; i++)
   for(int j=0; j<NMO; j++,k++)
     diff += std::abs(vHS_manual(i,j)-vHS2[k]);
  std::cout<<diff <<std::endl;

  //mkl_set_num_threads(1);
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    Timer.start(str.c_str());
    mkl_zcsrmv( &trans, &nrows, &ncols, &one, matdes, vn_1ddata2.data() , vn_1dindx.data(),  vn_1drows.data() ,  &(vn_1drows[1]), sigma2.data(), &zero, vHS2.data() );
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  std::cout<<"\n\n\n**********************************************************\n";
  std::cout<<"     Testing local energy   \n";
  std::cout<<"**********************************************************\n";


  str="epot_manual";
  itG = D0.begin();
  end2 = Vuv.end();
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    epot_manual = 0;
    Timer.start(str.c_str());  
    for(s2Dit it = Vuv.begin(); it != end2; it++)
      epot_manual += (*(itG + std::get<0>(*it))) * (*(itG + std::get<1>(*it))) * std::get<2>(*it);
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  vn_1ddata2.resize(Vuv.size()); 
  vn_1dindx.resize(Vuv.size());
  vn_1drows.resize(2*NMO*NMO+1);
  vHS2.resize(2*NMO*NMO);
  int curr=-1;
  for(int n=0; n<Vuv.size(); n++) {
    if( std::get<0>(Vuv[n]) != curr ) { 
      int old = curr;  
      curr = std::get<0>(Vuv[n]);
      for(int i=old+1; i<=curr; i++) vn_1drows[i] = n;
    }
    vn_1ddata2[n] = std::get<2>(Vuv[n]);
    vn_1dindx[n] = std::get<1>(Vuv[n]); 
  }
  for(int i=curr+1; i<vn_1drows.size(); i++)
    vn_1drows[i] = vn_1ddata2.size();

  int one_int = 1;
  trans = 'N';
  nrows = ncols = 2*NMO*NMO;
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmv( &trans, &nrows, &ncols, &one, matdes, vn_1ddata2.data() , vn_1dindx.data(),  vn_1drows.data() ,  &(vn_1drows[1]), D0_mkl.data(), &zero, vHS2.data() );
  ComplexType epot_mkl = 0; 
  for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i];
  std::cout<<"Difference between epot_manual and epot_mkl: " <<epot_manual-epot_mkl <<std::endl;

  str="epot_mkl_zcsrmv";
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    epot_mkl = 0;
    Timer.start(str.c_str());
    mkl_zcsrmv( &trans, &nrows, &ncols, &one, matdes, vn_1ddata2.data() , vn_1dindx.data(),  vn_1drows.data() ,  &(vn_1drows[1]), D0_mkl.data(), &zero, vHS2.data() );
    for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i];
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

  mkl_cspblas_zcsrgemv (&trans, &nrows, vn_1ddata2.data() , vn_1drows.data(), vn_1dindx.data(),  D0_mkl.data(), vHS2.data());
  epot_mkl = 0;
  for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i]; 
  std::cout<<"Difference between epot_manual and epot_cspblas_zcsrgemv: " <<epot_manual-epot_mkl <<std::endl;

  str="epot_mkl_cspblas_zcsrgemv";
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    epot_mkl = 0;
    Timer.start(str.c_str());
    mkl_cspblas_zcsrgemv (&trans, &nrows, vn_1ddata2.data() , vn_1drows.data(), vn_1dindx.data(),  D0_mkl.data(), vHS2.data());
    for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i]; 
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

/*
  std::vector<std::complex<float> > vn_1ddata4(vn_1ddata2.size());
  std::vector<std::complex<float> > vHS4(vHS2.size());
  std::vector<std::complex<float> > D0_mkl_float(D0_mkl.size());
  for(int i=0; i<vn_1ddata4.size(); i++) vn_1ddata4[i] = static_cast<std::complex<float> >(vn_1ddata2[i]);
  for(int i=0; i<vHS4.size(); i++) vHS4[i] = 0;
  for(int i=0; i<D0_mkl_float.size(); i++) D0_mkl_float[i] = static_cast<std::complex<float> >(D0_mkl[i]);
  mkl_cspblas_ccsrgemv (&trans, &nrows, vn_1ddata4.data() , vn_1drows.data(), vn_1dindx.data(),  D0_mkl_float.data(), vHS4.data());
  std::complex<float> epot_mkl_float = 0;
  for(int i=0; i<2*NMO*NMO; i++) epot_mkl_float += vHS4[i] * D0_mkl_float[i];
  std::cout<<"Difference between epot_manual and epot_cspblas_zcsrgemv(float): " <<epot_manual-epot_mkl_float <<std::endl;

  str="epot_mkl_cspblas_zcsrgemv(float)";
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    epot_mkl_float = 0;
    Timer.start(str.c_str());
    //mkl_cspblas_ccsrgemv (&trans, &nrows, vn_1ddata4.data() , vn_1drows.data(), vn_1dindx.data(),  D0_mkl_float.data(), vHS4.data());
    for(int i=0; i<2*NMO*NMO; i++) epot_mkl_float += vHS4[i] * D0_mkl_float[i];
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;
 */ 

  std::cout<<"\n\n\n**********************************************************\n";
  std::cout<<"     Testing SparseMatrix class   \n";
  std::cout<<"**********************************************************\n";

  std::cout<<"\n\n\n**********************************************************\n";
  std::cout<<"     Testing  potential energy   \n";
  std::cout<<"**********************************************************\n";

  ComplexSpMat SpVuv(2*NMO*NMO); 
  SpVuv.reserve(Vuv.size());
  for(int n=0; n<Vuv.size(); n++) {
    ComplexType t(std::get<2>(Vuv[n]));
    SpVuv.add( std::get<0>(Vuv[n]) , std::get<1>(Vuv[n]), t); 
  }
  SpVuv.compress();

  mkl_cspblas_zcsrgemv (&trans, &nrows, SpVuv.values() , SpVuv.row_index(), SpVuv.column_data(),  D0_mkl.data(), vHS2.data());
  epot_mkl = 0;
  for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i];
  std::cout<<"Difference between epot_manual and epot_spmat_cspblas_zcsrgemv: " <<epot_manual-epot_mkl <<std::endl;

  str="epot_spmat_mkl_cspblas_zcsrgemv";
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    epot_mkl = 0;
    Timer.start(str.c_str());
    mkl_cspblas_zcsrgemv (&trans, &nrows, SpVuv.values() , SpVuv.row_index(), SpVuv.column_data(),  D0_mkl.data(), vHS2.data());
    for(int i=0; i<2*NMO*NMO; i++) epot_mkl += vHS2[i] * D0_mkl[i];
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;

/*  
  std::cout<<std::endl <<std::endl <<" EIGEN: \n\n";
  std::vector<Eigen::SparseMatrix<ComplexType> > vn_eigen;
  Eigen::Matrix<ComplexType,Eigen::Dynamic, Eigen::Dynamic> eig_vHS(2*NMO,NMO);   
  vn_eigen.resize(vn_indx.size()-1);  
  for(int i=0; i<vn_indx.size()-1; i++) {
    vn_eigen[i].resize(2*NMO,NMO);
    vn_eigen[i].reserve(vn_indx[i+1]-vn_indx[i]);
    for(int n = vn_indx[i]; n<vn_indx[i+1]; n++) 
      vn_eigen[i].insert(std::get<0>( vn[n] ),std::get<1>( vn[n] )) = std::get<2>( vn[n] );
    vn_eigen[i].makeCompressed();
  }

  str="vn_Eigen";
  Timer.reset(str.c_str());
  for(int nt=0; nt<ntimes; nt++) {
    eig_vHS.setZero();
    Timer.start(str.c_str());
    for(int i=0; i<vn_eigen.size(); i++)
      eig_vHS += sigma[i]*vn_eigen[i]; 
    Timer.stop(str.c_str());
  }
  std::cout<<str <<":  " <<Timer.average(str.c_str()) <<std::endl <<std::endl;  
*/

#endif

}

}
