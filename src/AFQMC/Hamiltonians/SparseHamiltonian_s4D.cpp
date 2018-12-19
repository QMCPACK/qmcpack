#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include <Platforms/sysutil.h>
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/SimpleParser.h"
#include "Configuration.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D.h"
#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Matrix/hdf5_readers.hpp"
#include "AFQMC/Matrix/array_partition.hpp"

#include "AFQMC/Matrix/csr_matrix.hpp"
//#include "AFQMC/Matrix/spma_communications.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
   SparseHamiltonian_s4D::shm_csr_matrix SparseHamiltonian_s4D::calculateHSPotentials(RealType cut, ComplexMatrix& vn0, TaskGroup& TGprop) {

     return shm_csr_matrix({0,0},0,boost::mpi3::intranode::allocator<SPComplexType>(TG.Node()));

   }
*/
/*
  void SparseHamiltonian_s4D::calculateHSPotentials(RealType cut, RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool parallel)
  {

    if(skip_V2)
      APP_ABORT("Error: Calling SparseHamiltonian_s4D routines with skip_V2=yes. \n\n\n");

    if(distribute_Ham && !parallel) {
        APP_ABORT("Error: Distributed hamiltonian requires parallel Cholesky factorization. \n\n\n");
    }

    int rnk=0;
    rnk = TG.getGlobalRank();
    //cholesky_residuals.clear();
    //cholesky_residuals.reserve(2*NMO*NMO);

      // 0. vn0 = -0.5* sum_{i,l,sigma} (sum_j <i_sigma,j_sigma|j_sigma,l_sigma> ) c+i_sigma cl_sigma 
      assert(vn0.rows() >= NMO);
      assert(vn0.cols() == NMO);
      std::vector<s4D<ValueType> > v2sym; 
      vn0 = ValueType(0);
      {
        long nt=0,rk=0,npr=1; 
        if(distribute_Ham) {
          npr = static_cast<long>(TG.getTGSize());
          rk = static_cast<long>(TG.getTGRank());
        } else if(parallel) { 
          npr = static_cast<long>(TG.getGlobalSize());
          rk = static_cast<long>(TG.getGlobalRank());
        }
        for(s4Dit it = V2.begin(); it != V2.end(); it++, nt++) {
          if( nt%npr != rk) continue;
          // dumb and slow, but easy for now
          find_equivalent_OneBar_for_integral_list(*it,v2sym); 
          for(int n=0; n<v2sym.size(); n++) {   
            IndexType i,j,k,l;
            ValueType V;
            std::tie (i,j,k,l,V) = v2sym[n];   
            int sector = getSpinSector(NMO,i,j,k,l);
            // <i,j | j,l> -> v(i,l)  
            if( (sector==0 || sector==3) && j==k) 
              vn0(i,Index2Col(NMO,l)) -= 0.5*V;
          }
        }  
        if(parallel) {
          std::vector<ComplexType> g(NMO*NMO);
          std::copy(vn0.begin(),vn0.begin()+NMO*NMO,g.begin());
          TG.Global().all_reduce(g.begin(),g.end(),vn0.begin());
        }
      }
*/
      /********************************************************************
      *               Calculate Cholesky decomposition 
      *
      *   1. The mapping of the 2-el repulsion integrals to a 2-D matrix
      *      is done as follows:
      *         V(i,j,k,l) -->  V( i*NMO+k,  l*NMO+j )
      ********************************************************************/
/*
     Timer.reset("Generic");
     Timer.start("Generic");

     // used to split (i,k) space over all processors for the construction and storage of cholesky vectors
     int npr = TG.getGlobalSize(), rk = TG.getGlobalRank(); 
     int ik0=0,ik1=NMO*NMO;
     std::vector<int> ik_partition;
     // used for split of (i,k) space over the processors in a TG during the evaluation of H(i,kmax,k,imax) piece
     // in case the Hamiltonian is distributed
     int tg_npr = TG.getTGSize(), tg_rk = TG.getTGRank(); 
     int tg_ik0=0,tg_ik1=NMO*NMO;
     std::vector<int> tgik_partition;
     std::vector<int> cnts;
     if(parallel) {
       FairDivide(NMO*NMO,npr,ik_partition); 
       ik0 = ik_partition[rk];
       ik1 = ik_partition[rk+1];
       cnts.resize(npr);
       for(int i=0; i<npr; i++)
         cnts[i] = ik_partition[i+1]-ik_partition[i];
     } else {
       cnts.resize(1);
       cnts[0] = ik1-ik0; 
     }
     if(distribute_Ham) {
       FairDivide(NMO*NMO,tg_npr,tgik_partition);
       tg_ik0 = tgik_partition[tg_rk];
       tg_ik1 = tgik_partition[tg_rk+1];       
     }
     int nterms = ik1-ik0; 
     int maxnterms = *std::max_element(cnts.begin(),cnts.end()); 
     if(nterms < 1) {
       APP_ABORT("Error: Too many processors in parallel calculation of HS potential. Try reducing the number of cores, calculating in serial or reading from a file. \n\n\n ");
     }

     // will store full Cholesky std::vectors now and keep sparse versions only at the end.
     // want to avoid possible numerical issues from truncation
     std::vector< std::vector<ValueType> > L;
     L.reserve(NMO*NMO);

     std::vector<ValueType> Lnmax(NMO*NMO);
     std::vector<s2D<ComplexType> > IKLmax(npr); 
     s2D<ComplexType> mymax;
     int maxloc=0;

     std::vector<ValueType> Lcomm, Lcomm2;
     if(distribute_Ham) Lcomm.resize(NMO*NMO);

     // to store diagonal elements to avoid search, since they are used often 
     std::vector<ValueType> Duv(nterms);
     if(distribute_Ham) {
       std::fill(Lcomm.begin(),Lcomm.end(),ValueType(0));
       for(IndexType i=0, nt=0; i<NMO; i++) {
         for(IndexType k=0; k<NMO; k++,nt++) {
           if(nt<tg_ik0 || nt>=tg_ik1) continue;
           // <i,k|k,i> 
           if( k<i ) {
#if defined(QMC_COMPLEX)
             // <k,i|i,k> 
             if(k>=min_i && k<max_i) Lcomm[nt] = H(k,i,i,k);
#else
             // <k,k|i,i> 
             if(k>=min_i && k<max_i) Lcomm[nt] = H(k,k,i,i);
#endif
           } else {
#if defined(QMC_COMPLEX)
             // <i,k|k,i> 
             if(i>=min_i && i<max_i) Lcomm[nt] = H(i,k,k,i);
#else
             // <i,i|k,k> 
             if(i>=min_i && i<max_i) Lcomm[nt] = H(i,i,k,k);
#endif
           }
#if defined(QMC_COMPLEX)
           if(Lcomm[nt].imag() > 1e-10 || Lcomm[nt].real() < RealType(0)) {
             app_error()<<" WARNING: Found negative/complex Duv: " <<i <<" " <<k <<" " <<Lcomm[nt] <<std::endl;
             if(zero_bad_diag_2eints)
               Lcomm[nt] = ValueType(0);
           }
#else
           if(Lcomm[nt] < ValueType(0)) {
             app_error()<<" WARNING: Found negative Duv: " <<i <<" " <<k <<" " <<Lcomm[nt] <<std::endl;
             if(zero_bad_diag_2eints)
               Lcomm[nt] = ValueType(0);
           }
#endif
         }
         if(nt>=tg_ik1) break;
       }
       {
         Lcomm2 = Lcomm;
         TG.Global().all_reduce(Lcomm2.begin(),Lcomm2.end(),Lcomm.begin()); 
       }  
       std::copy( Lcomm.begin()+ik0, Lcomm.begin()+ik1, Duv.begin() );
     } else {
       for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
         for(IndexType k=0; k<NMO; k++,nt++) {
           if(nt<ik0 || nt>=ik1) continue;
           Duv[ik] = H(i,k,k,i);  
#ifndef QMC_COMPLEX
           if(Duv[ik] < ValueType(0)) {
             app_error()<<" WARNING: Found negative Duv: " <<i <<" " <<k <<" " <<Duv[ik] <<std::endl;  
             if(zero_bad_diag_2eints) 
               Duv[ik] = ValueType(0);
           }
#else
           if(Duv[ik].imag() > 1e-10 || Duv[ik].real() < RealType(0)) {
             app_error()<<" WARNING: Found negative/complex Duv: " <<i <<" " <<k <<" " <<Duv[ik] <<std::endl;
             if(zero_bad_diag_2eints)
               Duv[ik] = ValueType(0);
           }
#endif
           ik++;
         }
         if(nt>=ik1) break;
       }
     }

     // D(ik,lj) = H(i,j,k,l) - sum_p Lp(ik) Lp*(lj)
     // Diagonal:  D(ik,ik) = H(i,k,k,i) - sum_p Lp(ik) Lp*(ik) 
     RealType max=0;
     IndexType ii=-1,kk=-1;
     mymax = std::make_tuple(-1,-1,ComplexType(0));
     for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
      for(IndexType k=0; k<NMO; k++,nt++) {
        if(nt<ik0 || nt>=ik1) continue;
        if( std::abs(Duv[ik]) > max) {
          //max = std::get<2>(mymax) =std::abs(Duv[ik]);  
          max = std::abs(Duv[ik]);  
          std::get<2>(mymax) = Duv[ik];
          ii=std::get<0>(mymax)=i;
          kk=std::get<1>(mymax)=k;
        } 
        ik++;
      }
      if(nt>=ik1) break;
     }
     if(ii<0 || kk<0) {
      app_error()<<"Problems with Cholesky decomposition. \n";
      APP_ABORT("Problems with Cholesky decomposition. \n");   
     }

     if(parallel) {
       // test this call directly through mpi3 later
       MPI_Allgather(reinterpret_cast<char*>(&mymax),sizeof(s2D<ComplexType>),MPI_CHAR,reinterpret_cast<char*>(IKLmax.data()),sizeof(s2D<ComplexType>),MPI_CHAR,TG.Global().impl_);
       ii=-1;kk=-1;
       max=0;
       for(int i=0; i<npr; i++) {
         if(std::abs(std::get<2>(IKLmax[i]))>max) {
           ii=std::get<0>(IKLmax[i]);
           kk=std::get<1>(IKLmax[i]);
           max=std::abs(std::get<2>(IKLmax[i]));
           maxloc=i;
         } 
       }
     }

     if(printEig) {
      app_log()<<"Residuals of Cholesky factorization at each iteration: \n";
      app_log()<<L.size() <<" " <<std::abs(max) <<"\n";
     }

     Timer.reset("Generic1");
     Timer.reset("Generic2");

     if(test_2eint && !distribute_Ham) {
*/
       /* <ij|kl> <--> M_(ik)_(lj) where M must be positive definite.
        * --> |M_(ik)_(lj)| <= sqrt( |M_(ik)_(ik)| |M_(lj)_(lj)| )
        * --> |<ij|kl>| <= sqrt( |<ik|ki>| |<lj|jl| )   
        *  
        *  MAM: Being careful to check influence of cutoff in check!
        */
/*
       for(s4Dit it = V2.begin(); it != V2.end(); it++) {
         IndexType i,j,k,l;
         ValueType w1,w2,w3;
         std::tie (i,j,k,l,w1) = *it;  
 
         w2 = H(i,k,k,i);
         w3 = H(l,j,j,l);
         ValueType w4 = H(j,l,l,j);
         if( std::abs(w1) > std::sqrt(std::abs(w2*w3)) ) {
           
           if(std::abs(w2) == 0.0 && std::abs(w3) > 0.0) {
             
              // <ik|ki> was truncated, so positive-definiteness is broken
              if( std::abs(w1)*std::abs(w1)/std::abs(w3) < cutoff1bar ) continue;

           } else if(std::abs(w3) == 0.0 && std::abs(w2) > 0.0) { 

              // <lj|jl> was truncated, so positive-definiteness is broken
              if( std::abs(w1)*std::abs(w1)/std::abs(w2) < cutoff1bar ) continue;

           } 

           app_log()<<" Problems with positive-definiteness: " 
                    <<i <<" " <<j <<" " <<k <<" " <<l <<" " 
                    <<w1 <<" " <<w2 <<" " <<w3 <<" " <<w4 <<std::endl; 
         }

       } 

     }

     int cnt_energy_increases=0;
     RealType max_old;
     while(max > cutoff_cholesky) {

       Timer.start("Generic1");
       RealType oneOverMax = 1/std::sqrt(std::abs(max));

       //cholesky_residuals.push_back(std::abs(max));

       // calculate new cholesky std::vector based on (ii,kk)
       L.push_back(std::vector<ValueType>(maxnterms));  
       std::vector<ValueType>& Ln = L.back();
       std::vector<ValueType>::iterator it = Ln.begin();

       if(rk==maxloc) {
         for(int n=0; n<L.size()-1; n++)
           Lnmax[n] = L[n][ii*NMO+kk-ik0]; 
       }
       if(parallel && L.size()>1) 
         TG.Global().broadcast_n(Lnmax.begin(),L.size()-1,maxloc);

       if(distribute_Ham) {
         std::fill(Lcomm.begin(),Lcomm.end(),ValueType(0));
         s4D<ValueType> s;
         for(IndexType i=0, nt=0; i<NMO; i++) {
           for(IndexType k=0; k<NMO; k++, nt++) {
             if(nt<tg_ik0 || nt>=tg_ik1) continue;
             s = std::make_tuple(i,kk,k,ii,ValueType(0));
             bool cjgt = find_smallest_permutation(s);
             if( std::get<0>(s) >= min_i && std::get<0>(s) < max_i ) {
               s4Dit it_ = std::lower_bound( V2.begin(), V2.end(), s, 
                 [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs)
                 {
                   return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
                 }
               );
               if (it_ != V2.end() &&  std::get<0>(*it_)==std::get<0>(s) &&  std::get<1>(*it_)==std::get<1>(s) && std::get<2>(*it_)==std::get<2>(s) && std::get<3>(*it_)==std::get<3>(s) ) {
#if defined(QMC_COMPLEX)
                 if(cjgt)
                   Lcomm[nt] = std::conj(std::get<4>(*it_));
                 else
#endif
                   Lcomm[nt] = std::get<4>(*it_);
               }
             } 
           }
           if(nt>=tg_ik1) break;
         }
         {
           Lcomm2 = Lcomm;
           TG.Global().all_reduce(Lcomm2.begin(),Lcomm2.end(),Lcomm.begin()); 
         }  
         std::copy( Lcomm.begin()+ik0, Lcomm.begin()+ik1, it );         
       } else { 
         for(IndexType i=0, nt=0; i<NMO; i++) {
           for(IndexType k=0; k<NMO; k++, nt++) {
             if(nt<ik0 || nt>=ik1) continue;
             *(it++) = H(i,kk,k,ii);
           }
           if(nt>=ik1) break;
         }
       }

       for(int n=0; n<L.size()-1; n++) {
         //ValueType scl = myconj(L[n][ii*NMO+kk]); 
         ValueType scl = myconj(Lnmax[n]); 
         std::vector<ValueType>::iterator it1 = L[n].begin();
         it = Ln.begin();
         for(IndexType i=0; i<nterms; i++) 
           *(it++) -= *(it1++)*scl; 
       }

       it = Ln.begin();
       for(IndexType i=0; i<nterms; i++)
         *(it++) *= oneOverMax;
       Timer.stop("Generic1");
       
       Timer.start("Generic2");
       max_old = max;
       IndexType ii0=ii,kk0=kk;
       max=0;
       ii=-1;
       kk=-1;
       mymax = std::make_tuple(-1,-1,ComplexType(0));
       for(IndexType i=0,ik=0,nt=0; i<NMO; i++) {
        for(IndexType k=0; k<NMO; k++,nt++) {
         if(nt<ik0 || nt>=ik1) continue;
         Duv[ik] -= Ln[ik]*myconj(Ln[ik]);  
         if(zero_bad_diag_2eints) {
           if( std::abs(Duv[ik]) > max && toComplex(Duv[ik]).real() > 0) {
             //max = std::get<2>(mymax) =std::abs(Duv[ik]);  
             max = std::abs(Duv[ik]);  
             std::get<2>(mymax) = Duv[ik];
             ii=std::get<0>(mymax)=i;
             kk=std::get<1>(mymax)=k;
           }
         } else {
           if( std::abs(Duv[ik]) > max) {
             //max = std::get<2>(mymax) =std::abs(Duv[ik]);  
             max = std::abs(Duv[ik]);  
             std::get<2>(mymax) = Duv[ik];
             ii=std::get<0>(mymax)=i;
             kk=std::get<1>(mymax)=k;
           }
         }
         ik++;
        }
        if(nt>=ik1) break;
       }
       if(parallel) {
         MPI_Allgather(reinterpret_cast<char*>(&mymax),sizeof(s2D<ComplexType>),MPI_CHAR,reinterpret_cast<char*>(IKLmax.data()),sizeof(s2D<ComplexType>),MPI_CHAR,TG.Global().impl_);
         ii=-1;kk=-1;
         max=0;
         for(int i=0; i<npr; i++) {
           if(std::abs(std::get<2>(IKLmax[i]))>max) {
             ii=std::get<0>(IKLmax[i]);
             kk=std::get<1>(IKLmax[i]);
             max=std::abs(std::get<2>(IKLmax[i]));
             maxloc=i;
           }
         }
       }
       if(ii<0 || kk<0) {
        app_error()<<"Problems with Cholesky decomposition. \n";
        APP_ABORT("Problems with Cholesky decomposition. \n");
       }
       Timer.stop("Generic2");
       if( TG.getGlobalRank()==0 && toComplex(std::get<2>(IKLmax[maxloc])).real() < 0.0 ) {
         app_error()<<"\n Warning: Error matrix is not positive definite in Cholesky factorization" <<std::get<2>(IKLmax[maxloc])  <<"\n"
                    <<" Factorization is likely to fail. \n"
                    <<" This likely happens because: \n"
                    <<"      1) Cutoff of 2-electron integrals is too agresive, \n"
                    <<"      2) 2-electron integrals are incorrect. \n" 
                    <<" Consider using fix_2eint=true to remove orbital pairs that lead to break down of positive definiteness. \n"  
                    <<" Otherwise, 1) reduce cutoff_1bar (preferred) or 2) increase cutoff_cholesky. \n\n"; 
       }  
       if( TG.getGlobalRank()==0 && printEig)
         app_log()<<L.size() <<" " <<ii <<" " <<kk <<" " <<std::abs(max) <<" " <<Timer.total("Generic1") <<" " <<Timer.total("Generic2")   <<"\n";
       if(max > max_old) {
         cnt_energy_increases++;
         if(cnt_energy_increases == 3) { 
           app_error()<<"ERROR: Problems with convergence of Cholesky decomposition. \n" 
             <<"Number of std::vectors found so far: " <<L.size() <<"\n"
             <<"Current value of truncation error: " <<max_old <<" " <<max <<std::endl;  
             APP_ABORT("Problems with convergence of Cholesky decomposition.\n"); 
         }
       }
     }
     app_log()<<" Found: " <<L.size() <<" Cholesky std::vectors with a cutoff of: " <<cutoff_cholesky <<"\n";   

     Timer.stop("Generic");
     app_log()<<" -- Time to generate Cholesky factorization: " <<Timer.average("Generic") <<std::endl;

     if(test_breakup && !parallel && !distribute_Ham) {

      app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0,nt=0,ik=0; i<NMO; i++)
       for(IndexType j=0; j<NMO; j++) 
        for(IndexType k=0; k<NMO; k++)
         for(IndexType l=0; l<NMO; l++,nt++) {     
           if(nt<ik0||nt>=ik1) continue;
           ValueType v2 = H(i,j,k,l);
           ValueType v2c = 0.0;
           // is it L*L or LL*???
           for(int n=0; n<L.size(); n++) v2c += L[n][i*NMO+k]*myconj(L[n][l*NMO+j]);
           s+=std::abs(v2-v2c);
           if( max < std::abs(v2-v2c) ) max = std::abs(v2-v2c); 
           if( std::abs(v2-v2c) > 10*cutoff_cholesky ) {
             app_error()<<" Problems with Cholesky decomposition, i,j,k,l,H2,H2c: "
                       <<i <<" "
                       <<j <<" "
                       <<k <<" "
                       <<l <<" "
                       <<v2 <<" "
                       <<v2c <<std::endl;
           }
           ik++;
         }
      app_log()<<"\n ********************************************\n Average error due to truncated Cholesky factorization (in units of cutoff), maximum error   : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO <<"  " <<max <<" \n********************************************\n"<<std::endl; 

       Timer.stop("Generic");
       app_log()<<" -- Time to test Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     }
*/
      /********************************************************************
      *  You get 2 potentials per Cholesky std::vector   
      *
      *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
      ********************************************************************/
/*
      Timer.reset("Generic");
      Timer.start("Generic");

      ValueType sqrtdt = std::sqrt(dt)*0.5;

      std::vector<int> cnt_per_vec(2*L.size());
      std::vector<int> ik2padded;
      if(parallel) { 
        Lcomm.resize(npr*maxnterms); 
        ik2padded.resize(npr*maxnterms);
        for(int i=0; i<ik2padded.size(); i++) ik2padded[i]=i;
        // to make it independent of behavior of FairDivide
        const int tag = npr*maxnterms+10000;
        std::vector<int>::iterator it = ik2padded.begin(), itc = cnts.begin();
        for(; itc!=cnts.end(); itc++, it+=maxnterms) 
          std::fill(it+(*itc), it+maxnterms,tag); 
        std::stable_partition( ik2padded.begin(), ik2padded.end(), 
            [tag] (const int& i) { return i<tag; }
              ); 
      } else {
        ik2padded.resize(NMO*NMO);
        for(int i=0; i<NMO*NMO; i++) ik2padded[i]=i;
      }  

      Timer.reset("Generic2");
      Timer.start("Generic2");
      if(!sparse) cut=1e-8;
      else if(cut < 1e-12) cut=1e-12;
      int cnt=0, cntn=0;
      // generate sparse version
      int nvecs=L.size();
      for(int n=0; n<L.size(); n++) { 
        ValueType* Ls; 
        if(parallel) {
          TG.Global().gather_n(L[n].begin(),maxnterms,Lcomm.begin(),0);
          if(TG.getGlobalRank()==0) Ls = Lcomm.data();
        } else {
          Ls = L[n].data();
        } 
        if(TG.getGlobalRank()==0) {
          int np=0, nm=0;
          for(IndexType i=0; i<NMO; i++) 
           for(IndexType k=0; k<NMO; k++) { 
             // v+
             if(std::abs( (Ls[ik2padded[i*NMO+k]] + myconj(Ls[ik2padded[k*NMO+i]])) ) > cut) 
               np++;
             // v-
             if(std::abs( (Ls[ik2padded[i*NMO+k]] - myconj(Ls[ik2padded[k*NMO+i]])) ) > cut)  
               nm++;
           }
          cnt_per_vec[2*n] = np;
          cnt_per_vec[2*n+1] = nm;
          if(nm>0) cntn++;
        }
        if(n>0 && (n%100==0))
          TG.global_barrier();
      } 

      int nnodes = TGprop.getNNodesPerTG();
      int cv0=0,cvN=2*L.size();
      std::vector<int> sets;
      std::vector<int> sz_per_node(nnodes);
      int node_number = TGprop.getLocalNodeNumber();
      if(parallel) {
        TG.Global().broadcast(cnt_per_vec.begin(),cnt_per_vec.end());
        int nvec = std::count_if(cnt_per_vec.begin(),cnt_per_vec.end(),
                 [] (int i) { return i>0; } ); 
        if(sparse) Spvn.setDims(NMO*NMO,nvec);

        nvec_per_node.resize(nnodes);
        if(nnodes==1) {
          cv0=0;
          cvN=2*L.size();
          cnt = std::accumulate(cnt_per_vec.begin(),cnt_per_vec.end(),0);
          if(sparse) Spvn.reserve(cnt);
          else {
            Dvn.allocate(NMO*NMO*nvec);
            Dvn.resize(NMO*NMO*nvec);
          }
          nvec_per_node[0] = nvec; 
          sz_per_node[0] = cnt;
        } else {
          sets.resize(nnodes+1);
          if(TG.getGlobalRank()==0) {
            // partition std::vectors over nodes in TG
            std::vector<int> blocks(2*L.size()+1); 
            blocks[0]=0;
            cnt=0;
            for(int i=0; i<2*L.size(); i++) {
              if(sparse) cnt+=cnt_per_vec[i];
              else cnt+= (cnt_per_vec[i]>0)?1:0;
              blocks[i+1] = cnt; 
            }
            balance_partition_ordered_set(2*L.size(),blocks.data(),sets);
            TG.Cores().broadcast(sets.begin(),sets.end());
            cv0 = sets[node_number];
            cvN = sets[node_number+1];
          
            // since many std::vectors might have zero size and will be discarded below,
            // count only non-zero
            for(int i=0; i<nnodes; i++) 
              nvec_per_node[i] = std::count_if(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],
                 [] (int i) { return i>0; } );      
            for(int i=0; i<nnodes; i++) 
              sz_per_node[i] = std::accumulate(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],0);
          } else if(TG.getCoreID()==0) {
            TG.Cores().broadcast(sets.begin(),sets.end());
            cv0 = sets[node_number];
            cvN = sets[node_number+1];
          }
          TG.Global().broadcast(nvec_per_node.begin(),nvec_per_node.end());  
          TG.Global().broadcast(sz_per_node.begin(),sz_per_node.end());
          if(sparse) Spvn.reserve(sz_per_node[node_number]); 
          else {
            Dvn.allocate(NMO*NMO*nvec_per_node[node_number]);
            Dvn.resize(NMO*NMO*nvec_per_node[node_number]);
          }
        }
      } else {
        int nvec=0;
        cnt=0; 
        for (int i : cnt_per_vec ) 
          if(i>0) { 
            cnt+=i;
            nvec++;
          } 
        nvec_per_node[0] = nvec;
        if(sparse) {
          Spvn.setDims(NMO*NMO,nvec);
          Spvn.allocate_serial(cnt);
        } else {
          Dvn.allocate_serial(NMO*NMO*nvec);
          Dvn.resize_serial(NMO*NMO*nvec);
        }
      } 

      int ncols = nvec_per_node[node_number];
      if(parallel) TG.global_barrier(); 

      Timer.stop("Generic2");
      app_log()<<"     -- setup: " <<Timer.average("Generic2") <<"\n";

      Timer.reset("Generic2");
      Timer.reset("Generic3");

#ifndef QMC_COMPLEX
      if(cntn>0)
        APP_ABORT("Found real Cholesky vectors with real integrals. This is not allowed with dense cholesky vectors. Run with sparse vectors or compile with complex integrals. \n\n\n");
#endif

      cnt=std::accumulate(nvec_per_node.begin(),nvec_per_node.begin()+node_number,0);
      for(int n=0; n<L.size(); n++) { 
       if( cnt_per_vec[2*n]==0 && cnt_per_vec[2*n+1]==0 ) continue;
       ValueType* Ls;
       Timer.start("Generic2");
       if(parallel) {
#if defined(QMC_COMPLEX)
         TG.Global().gather_n(L[n].begin(),maxnterms,Lcomm.begin(),0);
         if(TG.getCoreID()==0)
           TG.Cores().broadcast(Lcomm.begin(),Lcomm.end());  
#else
         TG.Global().all_gather_n(L[n].begin(),maxnterms,Lcomm.begin(),0);
#endif
         Ls = Lcomm.data();
       } else {
         Ls = L[n].data();
       }
       Timer.stop("Generic2");
       Timer.start("Generic3");
       if(TG.getCoreID()==0 && 2*n>=cv0 && 2*n<cvN && cnt_per_vec[2*n]>0) {
         int np=0;
         // v+
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) { 
           ValueType V = (Ls[ik2padded[i*NMO+k]] + myconj(Ls[ik2padded[k*NMO+i]])); 
           if(std::abs(V) > cut) { 
             V*=sqrtdt;
             if(sparse) Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(V));
             else Dvn[(i*NMO+k)*ncols + cnt]=static_cast<SPValueType>(V);
             ++np;
           } 
         }
         ++cnt;
         if(np==0) 
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n"); 
       }
       if(TG.getCoreID()==0 && (2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) {
#if defined(QMC_COMPLEX)
         int np=0;
         // v-
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) { 
           ValueType V = (Ls[ik2padded[i*NMO+k]] - myconj(Ls[ik2padded[k*NMO+i]]));
           if(std::abs(V) > cut) {
             V*=ComplexType(0.0,1.0)*sqrtdt;
             if(sparse) Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(V));
             else Dvn[(i*NMO+k)*ncols + cnt]=static_cast<SPValueType>(V);
             ++np;
           }
          }
         ++cnt;
         if(np==0) 
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n"); 
#else
       APP_ABORT("Error: This should not happen. Found negative cholesky vector. \n"); 
#endif
       }
       // necessary to avoid the avalanche of messages to the root from cores that are not head_of_nodes
       if(n>0 && (n%100==0))
         TG.global_barrier();
       Timer.stop("Generic3");
      }
      app_log()<<"     -- av comm time: " <<Timer.average("Generic2") <<"\n";
      app_log()<<"     -- av insert time: " <<Timer.average("Generic3") <<"\n";

      Timer.stop("Generic");
      app_log()<<" -- Time to assemble Cholesky Matrix: " <<Timer.average("Generic") <<"\n";

      if(TG.getGlobalRank()==0 && nnodes>1 && parallel) { 
        app_log()<<" Partition of Cholesky Vectors: 0 ";
        cnt=0;
        for(int i=0; i<nnodes; i++) { 
          cnt+=nvec_per_node[i]; 
          app_log()<<cnt <<" ";  
        }
        app_log()<<std::endl;
        if(sparse) {
          app_log()<<" Number of terms in Spvn per node in TG: ";
          for(int i : sz_per_node ) app_log()<<i <<" "; 
          app_log()<<std::endl;
        }
      }

      app_log()<<"Number of HS potentials: " <<ncols <<std::endl;
      if(sparse) {

        app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;
        app_log()<<"Compressing Spvn. \n";

        Timer.reset("Generic");
        Timer.start("Generic");

        if(parallel) Spvn.compress(TG.Node().impl_);
        else if(TG.getCoreID()==0) Spvn.compress();

        Timer.stop("Generic");
        app_log()<<"Done Compressing Spvn. \n";
        if(rnk==0) app_log()<<" -- Time to Compress Cholesky Matrix: " <<Timer.average("Generic") <<"\n";

      }

      if(parallel) TG.global_barrier(); 

     if(test_breakup) {

      if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      int* cols = Spvn.column_data();
      int* rows = Spvn.row_data();
      int* indx = Spvn.row_index();
      SPValueType* vals = Spvn.values();

      int NMO2 = NMO*NMO;

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<NMO; i++)
       for(IndexType j=0; j<NMO; j++)
        for(IndexType k=0; k<NMO; k++)
         for(IndexType l=0; l<NMO; l++) {
           ValueType v2 = H(i,j,k,l);
           int ik = i*NMO+k;
           int jl = j*NMO+l;
           ValueType v2c = SparseMatrixOperators::product_SpVSpV<ValueType>(indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[jl+1]-indx[jl],cols+indx[jl],vals+indx[jl]) / dt;
           s+=std::abs(v2-v2c);
           if( max < std::abs(v2-v2c) ) max = std::abs(v2-v2c);
           if( std::abs(v2-v2c) > 10*cutoff_cholesky ) {
             app_error()<<" Problems with H2 decomposition, i,j,k,l,H2,H2c: "
                       <<i <<" "
                       <<j <<" "
                       <<k <<" "
                       <<l <<" "
                       <<v2 <<" "
                       <<v2c <<std::endl;
           }
         }
      double scl = 1.0;
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/scl <<" " <<max <<" \n ********************************************\n"<<std::endl;

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }

  }

  bool SparseHamiltonian_s4D::createHamiltonianForPureDeterminant(int walker_type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b , std::vector<s1D<ValueType> >& hij, SPValueSMSpMat& Vijkl, const RealType cut)
  {

    // used to identify equal index sets (value is not compared)  
    _myEqv_snD_ myEqv;

    // teporary until mpi3 is fully integrated
    TaskGroup_ TGham(TG,"DummyHS",1,TG.getNCoresPerTG());

    createHij(walker_type,NMO,occ_a,occ_b,H1,hij,cut);

    long cnt2=0, number_of_terms=0;


      Timer.reset("Generic");
      Timer.start("Generic");

      std::vector<s4D<ValueType> > vs4D;
      vs4D.reserve(48);
      long N = NMO;

// right now, the algorithm will add all the terms associated with a given quartet (ijkl)
// at once. For the given ijkl, 3 (6) possible combinations can appear in the list for real (complex)  
//  Only do the smallest of the 3 (6) possible combinations to avoid duplicated 

#if defined(QMC_COMPLEX)
      std::vector<s4D<ValueType>> ineq_ijkl(6);
      std::vector<bool> setJ(6);
#else
      std::vector<s4D<ValueType>> ineq_ijkl(3);
      std::vector<bool> setJ(3);
#endif

      auto search_in_V2_IJ = [] (s4Dit it1, s4Dit it2,s4D<ValueType>& ijkl) {
          s4Dit first = std::lower_bound(it1,it2,ijkl,
             [] (const s4D<ValueType>& a, const s4D<ValueType>& b)
             {return (std::get<2>(a)<std::get<2>(b)) ||
                     (!(std::get<2>(b)<std::get<2>(a))&&(std::get<3>(a)<std::get<3>(b)));} );
          if (first!=it2 && (std::get<2>(ijkl)==std::get<2>(*first)) && (std::get<3>(ijkl)==std::get<3>(*first)))
            return std::make_tuple(std::get<4>(*first),true);
          return std::make_tuple(ValueType(0),false);
      };

      ValueType zero = ValueType(0);
      cnt2=0;
      long npr = TG.getGlobalSize(), rk = TG.getGlobalRank();
      OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
      ValueType J1,J2,J3,J1a=zero,J2a=zero,J3a=zero,fct;
      long p_min=0, p_max=IJ.size()-1;
// if integrals are distributed, 
//   distribute work over min_i/max_i sector among processors in the Hamiltonian TG
      if(distribute_Ham) {
        p_min = mapUT(static_cast<long>(min_i),static_cast<long>(min_i),N);
        if( max_i != NMO)   // FIX FIX FIX with SPIN_RESTRICTED 
          p_max = mapUT(static_cast<long>(max_i),static_cast<long>(max_i),N);
        // in this case, npr and rk are based on the extended TG (the one which includes all cores within a node)
        npr = TG.getTGSize();
        rk = TG.getTGRank();
      }
      for(long p=p_min, nt=0; p<p_max; p++) {
        // from n->m, I have all non-zero (k,l) for a given (i,j), with i<=j  
        long n = IJ[p];
        long m = IJ[p+1];
        if(n==m) continue;
        nt++;  // found good term, increase counter
        if( nt%npr != rk ) continue;
        s4Dit end = V2.begin()+m;
        for(s4Dit it = V2.begin()+n; it != end; it++) {
          // J1 = <ij|kl>   
          // J2 = <ij|lk> or <ik|lj>   
          // J3 = <ik|jl> or <il|jk>  
          std::tie (i,j,k,l,J1) = *it;

          IndexType occi,occj,occk,occl;
          occi = (occ_a[i]||occ_b[i])?1:0;
          occj = (occ_a[j]||occ_b[j])?1:0;
          occk = (occ_a[k]||occ_b[k])?1:0;
          occl = (occ_a[l]||occ_b[l])?1:0;
          if( occi+occj+occk+occl < 2) continue;

// the smallest permutation will satisfy  i<=j<=k<=l
// so need to generte smallest of the other sectors, as long as the smallest one is 
// in the list. You need to make sure that the smallest one in the list is processed
// Algorithm:
//  1. If (i,j,k,l) is the smallest permutation, keep going
//  2. If it is not, make a ordered list of all inequivalent permutations and record your position in the list 
//      a. For every term in the list before yours, 
//          i1. If the element is in the list, do nothing and continue to next element in V2.
//          i2. If not in the list, set the value of the appropriate Jx to zero and test next element in permutation list  
//  At the end you either do nothing because there is a permutationally inequivalent term smaller than you in the list, 
//  or you process the current term with all previous Js set to zero to avoid recalculating them.         


          if( i<=j && i<=k && i<=l && j<=k && j<=l && k<=l) {
            ineq_ijkl[0] = *it;
            setJ[0]=true;
            std::fill(setJ.begin()+1,setJ.end(),false);
#if defined(QMC_COMPLEX)
            ineq_ijkl[1] = std::forward_as_tuple(i,j,l,k,zero);
            ineq_ijkl[2] = std::forward_as_tuple(i,k,j,l,zero);
            ineq_ijkl[3] = std::forward_as_tuple(i,k,l,j,zero);
            ineq_ijkl[4] = std::forward_as_tuple(i,l,j,k,zero);
            ineq_ijkl[5] = std::forward_as_tuple(i,l,k,j,zero);
#else
            ineq_ijkl[1] = (j<=k)?(std::forward_as_tuple(i,j,l,k,zero)):(std::forward_as_tuple(i,k,l,j,zero));
            ineq_ijkl[2] = (k<=l)?(std::forward_as_tuple(i,k,j,l,zero)):(std::forward_as_tuple(i,l,j,k,zero));
#endif
          } else {

            // make sure there is a smaller one on the list
            // 1. generate smallest permutation
            OrbitalType i_=i,j_=j,k_=k,l_=l;
            if(j_<i_) std::swap(i_,j_);
            if(k_<i_) std::swap(i_,k_);
            if(l_<i_) std::swap(i_,l_);
            if(k_<j_) std::swap(j_,k_);
            if(l_<j_) std::swap(j_,l_);
            if(l_<k_) std::swap(k_,l_);
            std::fill(setJ.begin(),setJ.end(),false);
            ineq_ijkl[0] = std::forward_as_tuple(i_,j_,k_,l_,zero);
#if defined(QMC_COMPLEX)
            ineq_ijkl[1] = std::forward_as_tuple(i_,j_,l_,k_,zero);
            ineq_ijkl[2] = std::forward_as_tuple(i_,k_,j_,l_,zero);
            ineq_ijkl[3] = std::forward_as_tuple(i_,k_,l_,j_,zero);
            ineq_ijkl[4] = std::forward_as_tuple(i_,l_,j_,k_,zero);
            ineq_ijkl[5] = std::forward_as_tuple(i_,l_,k_,j_,zero);
#else
            ineq_ijkl[1] = (j_<=k_)?(std::forward_as_tuple(i_,j_,l_,k_,zero)):(std::forward_as_tuple(i_,k_,l_,j_,zero));
            ineq_ijkl[2] = (k_<=l_)?(std::forward_as_tuple(i_,k_,j_,l_,zero)):(std::forward_as_tuple(i_,l_,j_,k_,zero));
#endif
            bool process=false;
            for(int i=0; i<ineq_ijkl.size(); i++) {
              if( myEqv(ineq_ijkl[i],*it) ) {
                std::get<4>(ineq_ijkl[i])=J1;
                setJ[i]=true;
                process=true;
                break;
              } else {
                long p0 = mapUT(std::get<0>(ineq_ijkl[i]),std::get<1>(ineq_ijkl[i]),N);
                if(std::get<1>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[i]))) {
                  process=false;
                  break;
                } else {
                  setJ[i]=true;
                  std::get<4>(ineq_ijkl[i])=zero;
                }
              }
            }
            if(!process) continue;
          }

          // at this point, I know that:
          //    1. this is a term I must process 
          //    2. the smallest permutation is in ineq_ijkl[0] with J1=std::get<4>(ineq_ijkl[0])
          //    3. some values of Js might already be calculated based on setJ[k] 
          std::tie (i,j,k,l,J1) = ineq_ijkl[0];

          // look for <ij|lk>
          if(setJ[1]) {
            J2 = std::get<4>(ineq_ijkl[1]);
          } else if(i==j || l==k) {
            J2=J1;
          } else {
            long p0 = mapUT(i,std::get<1>(ineq_ijkl[1]),N);
            J2 = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[1]));
          }

          // look for <ik|jl>
          if(setJ[2]) {
            J3 = std::get<4>(ineq_ijkl[2]);
          } else if(j==k) {
            J3=J1;
          } else if(i==l) {
            J3 = myconj(J1);
          } else {
            long p0 = mapUT(i,std::get<1>(ineq_ijkl[2]),N);
            J3 = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[2]));
          }

#if defined(QMC_COMPLEX)
            J1a=J2a=J3a=zero;

            //  J2a = <ik|lj>
            if(setJ[3]) {
              J2a = std::get<4>(ineq_ijkl[3]);
            } else if(l==j) {
              J2a=J3;
            } else if(i==k) {
              J2a=J3;
            } else if(k==j) {
              J2a=J2;
            } else if(i==l) {
              J2a=std::conj(J2);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[3]),N);
              J2a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[3]));
            }

            //  J3a = <il|jk> 
            if(setJ[4]) {
              J3a = std::get<4>(ineq_ijkl[4]);
            } else if(l==j) {
              J3a=J2;
            } else if(i==k) {
              J3a=std::conj(J2);
            } else if(k==l) {
              J3a=J3;
            } else if(i==j) {
              J3a=std::conj(J3);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[4]),N);
              J3a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[4]));
            }

            //  For complex, there are 3 extra non-symmetric terms:
            //  J1a = <il|kj>
            if(setJ[5]) {
              J1a = std::get<4>(ineq_ijkl[5]);
            } else if(k==l) {
              J1a=J2a;
            } else if(i==j) {
              J1a=std::conj(J2a);
            } else if(j==k) {
              J1a=J3a;
            } else if(i==l) {
              J1a=J3a;
            } else if(l==j) {
              J1a=J1;
            } else if(i==k) {
              J1a=std::conj(J1);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[5]),N);
              J1a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[5]));
            }

#endif

          vs4D.clear();
          if(walker_type==0) {
            find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==1) {
            find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==2) {
              find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else {
            APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
          }
          cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);

        }
      }


      cnt2 = TG.Global().all_reduce_value(cnt2);
      Vijkl.reserve(cnt2+1000);
      cnt2=0;

      for(long p=p_min, nt=0; p<p_max; p++) {
        // from n->m, I have all non-zero (k,l) for a given (i,j), with i<=j  
        long n = IJ[p];
        long m = IJ[p+1];
        if(n==m) continue;
        nt++;  // found good term, increase counter
        if( nt%npr != rk ) continue;
        s4Dit end = V2.begin()+m;
        for(s4Dit it = V2.begin()+n; it != end; it++) {

          // for details see above
          std::tie (i,j,k,l,J1) = *it;

          IndexType occi,occj,occk,occl;
          occi = (occ_a[i]||occ_b[i])?1:0;
          occj = (occ_a[j]||occ_b[j])?1:0;
          occk = (occ_a[k]||occ_b[k])?1:0;
          occl = (occ_a[l]||occ_b[l])?1:0;
          if( occi+occj+occk+occl < 2) continue;

          if( i<=j && i<=k && i<=l && j<=k && j<=l && k<=l) {
            ineq_ijkl[0] = *it;
            setJ[0]=true;
            std::fill(setJ.begin()+1,setJ.end(),false);
#if defined(QMC_COMPLEX)
            ineq_ijkl[1] = std::forward_as_tuple(i,j,l,k,zero);
            ineq_ijkl[2] = std::forward_as_tuple(i,k,j,l,zero);
            ineq_ijkl[3] = std::forward_as_tuple(i,k,l,j,zero);
            ineq_ijkl[4] = std::forward_as_tuple(i,l,j,k,zero);
            ineq_ijkl[5] = std::forward_as_tuple(i,l,k,j,zero);
#else
            ineq_ijkl[1] = (j<=k)?(std::forward_as_tuple(i,j,l,k,zero)):(std::forward_as_tuple(i,k,l,j,zero));
            ineq_ijkl[2] = (k<=l)?(std::forward_as_tuple(i,k,j,l,zero)):(std::forward_as_tuple(i,l,j,k,zero));
#endif
          } else {

            OrbitalType i_=i,j_=j,k_=k,l_=l;
            if(j_<i_) std::swap(i_,j_);
            if(k_<i_) std::swap(i_,k_);
            if(l_<i_) std::swap(i_,l_);
            if(k_<j_) std::swap(j_,k_);
            if(l_<j_) std::swap(j_,l_);
            if(l_<k_) std::swap(k_,l_);
            std::fill(setJ.begin(),setJ.end(),false);
            ineq_ijkl[0] = std::forward_as_tuple(i_,j_,k_,l_,zero);
#if defined(QMC_COMPLEX)
            ineq_ijkl[1] = std::forward_as_tuple(i_,j_,l_,k_,zero);
            ineq_ijkl[2] = std::forward_as_tuple(i_,k_,j_,l_,zero);
            ineq_ijkl[3] = std::forward_as_tuple(i_,k_,l_,j_,zero);
            ineq_ijkl[4] = std::forward_as_tuple(i_,l_,j_,k_,zero);
            ineq_ijkl[5] = std::forward_as_tuple(i_,l_,k_,j_,zero);
#else
            ineq_ijkl[1] = (j_<=k_)?(std::forward_as_tuple(i_,j_,l_,k_,zero)):(std::forward_as_tuple(i_,k_,l_,j_,zero));
            ineq_ijkl[2] = (k_<=l_)?(std::forward_as_tuple(i_,k_,j_,l_,zero)):(std::forward_as_tuple(i_,l_,j_,k_,zero));
#endif
            bool process=false;
            for(int i=0; i<ineq_ijkl.size(); i++) {
              if( myEqv(ineq_ijkl[i],*it) ) {
                std::get<4>(ineq_ijkl[i])=J1;
                setJ[i]=true;
                process=true;
                break;
              } else {
                long p0 = mapUT(std::get<0>(ineq_ijkl[i]),std::get<1>(ineq_ijkl[i]),N);
                if(std::get<1>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[i]))) {
                  process=false;
                  break;
                } else {
                  setJ[i]=true;
                  std::get<4>(ineq_ijkl[i])=zero;
                }
              }
            }
            if(!process) continue;
          }

          std::tie (i,j,k,l,J1) = ineq_ijkl[0];

          // look for <ij|lk>
          if(setJ[1]) {
            J2 = std::get<4>(ineq_ijkl[1]);
          } else if(i==j || l==k) {
            J2=J1;
          } else {
            long p0 = mapUT(i,std::get<1>(ineq_ijkl[1]),N);
            J2 = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[1]));
          }

          // look for <ik|jl>
          if(setJ[2]) {
            J3 = std::get<4>(ineq_ijkl[2]);
          } else if(j==k) {
            J3=J1;
          } else if(i==l) {
            J3 = myconj(J1);
          } else {
            long p0 = mapUT(i,std::get<1>(ineq_ijkl[2]),N);
            J3 = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[2]));
          }

#if defined(QMC_COMPLEX)
            J1a=J2a=J3a=zero;

            //  J2a = <ik|lj>
            if(setJ[3]) {
              J2a = std::get<4>(ineq_ijkl[3]);
            } else if(l==j) {
              J2a=J3;
            } else if(i==k) {
              J2a=J3;
            } else if(k==j) {
              J2a=J2;
            } else if(i==l) {
              J2a=std::conj(J2);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[3]),N);
              J2a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[3]));
            }

            //  J3a = <il|jk> 
            if(setJ[4]) {
              J3a = std::get<4>(ineq_ijkl[4]);
            } else if(l==j) {
              J3a=J2;
            } else if(i==k) {
              J3a=std::conj(J2);
            } else if(k==l) {
              J3a=J3;
            } else if(i==j) {
              J3a=std::conj(J3);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[4]),N);
              J3a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[4]));
            }

            //  For complex, there are 3 extra non-symmetric terms:
            //  J1a = <il|kj>
            if(setJ[5]) {
              J1a = std::get<4>(ineq_ijkl[5]);
            } else if(k==l) {
              J1a=J2a;
            } else if(i==j) {
              J1a=std::conj(J2a);
            } else if(j==k) {
              J1a=J3a;
            } else if(i==l) {
              J1a=J3a;
            } else if(l==j) {
              J1a=J1;
            } else if(i==k) {
              J1a=std::conj(J1);
            } else {
              long p0 = mapUT(i,std::get<1>(ineq_ijkl[5]),N);
              J1a = std::get<0>(search_in_V2_IJ(V2.begin()+IJ[p0], V2.begin()+IJ[p0+1], ineq_ijkl[5]));
            }
#endif

          vs4D.clear();
          if(walker_type==0) {
            find_all_contributions_to_hamiltonian_closed_shell(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==1) {
            find_all_contributions_to_hamiltonian_collinear(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==2) {
              find_all_contributions_to_hamiltonian_ghf(NMO,aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else {
            APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
          }
          cnt2+=add_allowed_terms(NMO,vs4D,occ_a,occ_b, Vijkl, true, walker_type==2);
        }
      }

    TG.global_barrier();
    Timer.stop("Generic");
    app_log()<<"Time to generate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");
    redistribute_sparse_matrix(TGham,Vijkl);
    Timer.stop("Generic");
    app_log()<<"Time to communicate Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");

    if(!Vijkl.remove_repeated_and_compress(TG.Node().impl_)) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

    Timer.stop("Generic");
    app_log()<<"Time to remove_repeated_and_compress Hamiltonian: " <<Timer.total("Generic") <<std::endl;
    
    return true;

  }

 
  bool SparseHamiltonian_s4D::createHamiltonianForGeneralDeterminant(int walker_type, const ComplexMatrix& A,std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vijkl, const RealType cut)
  {
    if(skip_V2)
      APP_ABORT("Error: Calling SparseGeneralHamiltonian routines with skip_V2=yes. \n\n\n");

    // teporary until mpi3 is fully integrated
    TaskGroup_ TGham(TG,"DummyHS",1,TG.getNCoresPerTG());

    ComplexMatrix M,N;
    const ComplexType one = ComplexType(1.0);
    const ComplexType zero = ComplexType(0.0);
    int npr = TG.getGlobalSize(), rk = TG.getGlobalRank();

    if(walker_type == 0) {
      app_log()<<" Generating rotated hamiltonian matrices for RHF walker type. \n";
      if(A.rows() < NMO || A.cols() != NAEA) {
        app_error()<<" Error: Incorrect dimensions in Slater Matrix in  SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant: " <<A.rows() <<" " <<A.cols() <<std::endl;
        return false;
      }
    } else if(walker_type == 1) {
      app_log()<<" Generating rotated hamiltonian matrices for ROHF/UHF walker type. \n";
      // both UHF/GHF wavefunctions allowed
      if(A.rows() != 2*NMO || (A.cols() != NAEA && A.cols() != NAEA+NAEB)) {
        app_error()<<" Error: Incorrect dimensions in Slater Matrix in  SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant: " <<A.rows() <<" " <<A.cols() <<std::endl;
        return false;
      }
    } else if(walker_type==2) {
      app_log()<<" Generating rotated hamiltonian matrices for GHF walker type. \n";
      // both UHF/GHF wavefunctions allowed
      if(A.rows() != 2*NMO || (A.cols() != NAEA && A.cols() != NAEA+NAEB)) {
        app_error()<<" Error: Incorrect dimensions in Slater Matrix in  SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant: " <<A.rows() <<" " <<A.cols() <<std::endl;
        return false;
      }
      app_error()<<" Hamiltonian rotation not yet implemented for GHF type of walker. \n";
      return false;
    } else {
      app_error()<<" Error: Unacceptable walker_type in SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant: " <<walker_type <<std::endl;
      return false;
    }

    // half-rotate Hij to A basis
    // hij = sum_k A*(k,i) * H1(k,j) 
    rotateHij(walker_type,NMO,NAEA,NAEB,A,H1,hij,cut);

    if(!v2_full_init)
      V2_full.setup(TG.getCoreID()==0,std::string("SparseGeneralHamiltonian_V2_full"),TG.Node().impl_);
    rotateHijkl_s4D_integrals(distribute_Ham,walker_type,TG,NMO,NAEA,NAEB,A,V2,V2_full,KL,Vijkl,cut);

    Timer.reset("Generic");
    Timer.start("Generic");

    redistribute_sparse_matrix(TGham,Vijkl);
    TG.global_barrier();

    Timer.stop("Generic");
    app_log()<<"Time to communicate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;

#ifdef AFQMC_DEBUG
    app_log()<<" Done generating sparse hamiltonians. " <<std::endl;
    app_log()<<" Compressing sparse hamiltonians. " <<std::endl;
#endif

    if(!Vijkl.remove_repeated_and_compress(TG.Node().impl_)) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

#ifdef AFQMC_DEBUG
    app_log()<<" Done compressing sparse hamiltonians. " <<std::endl;
#endif

    TG.global_barrier();

    return true;
  }

*/

}

}

