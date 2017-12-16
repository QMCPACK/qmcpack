#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/SimpleParser.h"
#include "Configuration.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/ProjectorBase.h"
#include "AFQMC/Hamiltonians/DDProjector.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

namespace qmcplusplus
{

  bool DDProjector::initFromGuess()
  {
  
    // start by just setting to zero
    Pmat.resize(2*NMO,2*NMO); 
    Pmat=0;
  
    return true;
  }

  bool DDProjector::initFromHDF5(const std::string& fileName)
  {

    return true;

/*
    hdf_archive dump(myComm);
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      app_error()<<" Error opening integral file in SparseGeneralHamiltonian. \n";
      return false;
    }

    std::string path = "/Hamiltonian/SparseGeneralHamiltonian";
    if(!dump.is_group( path )) {
      app_error()<<" ERROR: H5Group /Hamiltonian/SparseGeneralHamiltonian does not exists in restart file. \n";
      return false;
    }

    if(!dump.push("Hamiltonian",false)) return false;
    if(!dump.push("SparseGeneralHamiltonian",false)) return false;

    std::vector<int> Idata(7);
    if(!dump.read(Idata,"dims")) return false;

    H1.resize(Idata[0]);
    V2.resize(Idata[1]);
    V2_2bar.resize(Idata[2]);

    if(NMO < 0) NMO = Idata[3];
    if(NAEA < 0) NAEA = Idata[4];
    if(NAEB < 0) NAEB = Idata[5];
    if(Idata[3] != NMO) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n"; 
      return false;
    }
    if(Idata[4] != NAEA) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n"; 
      return false;
    }
    if(Idata[5] != NAEB) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n"; 
      return false;
    }
    spinRestricted = (Idata[6]==0)?(true):(false);

    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);
    Idata.resize(NAEA+NAEB);
    if(!dump.read(Idata,"occups")) return false;
    for(int i=0; i<NAEA; i++) occup_alpha[i] = Idata[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = Idata[i];

    std::vector<double> Rdata(2);
    if(!dump.read(Rdata,"Energies")) return false;
    NuclearCoulombEnergy = Rdata[0];
    FrozenCoreEnergy = Rdata[0];

    int sz = std::max( 2*H1.size(), std::max( 4*V2.size(), 4*V2_2bar.size() ) );
    std::vector<IndexType> ivec;    
    ivec.reserve(sz);

    ivec.resize(2*H1.size());
    if(!dump.read(ivec,"H1_indx")) return false;
    for(int i=0, j=0; i<H1.size(); i++, j+=2)        
      H1[i] = std::make_tuple(ivec[j],ivec[j+1],0);  

    ivec.clear();
    ivec.resize(4*V2.size());
    if(!dump.read(ivec,"V2_indx")) return false;
    for(int i=0, j=0; i<V2.size(); i++, j+=4)
      V2[i] = std::make_tuple(ivec[j],ivec[j+1],ivec[j+2],ivec[j+3],0);  

    ivec.clear();
    ivec.resize(4*V2_2bar.size());
    if(!dump.read(ivec,"V2_2bar_indx")) return false;
    for(int i=0, j=0; i<V2_2bar.size(); i++, j+=4)
      V2_2bar[i] = std::make_tuple(ivec[j],ivec[j+1],ivec[j+2],ivec[j+3],0);  

    std::vector<IndexType>().swap(ivec);

    sz = std::max( H1.size(), std::max( V2.size(), V2_2bar.size() ) );
    std::vector<ValueType> vvec;
    vvec.reserve(sz); 

    vvec.resize(H1.size());
    if(!dump.read(vvec,"H1")) return false;
    for(int i=0; i<H1.size(); i++)           
      std::get<2>(H1[i]) = vvec[i]; 

    vvec.clear();
    vvec.resize(V2.size());
    if(!dump.read(vvec,"V2")) return false;
    for(int i=0; i<V2.size(); i++)
      std::get<4>(V2[i]) = vvec[i]; 

    vvec.clear();
    vvec.resize(V2_2bar.size());
    if(!dump.read(vvec,"V2_2bar")) return false;
    for(int i=0; i<V2_2bar.size(); i++)
      std::get<4>(V2_2bar[i]) = vvec[i]; 

    std::vector<ValueType>().swap(vvec);   

    dump.pop();
    dump.pop();

    dump.close();

    return true;
*/
  }

  void DDProjector::hdf_write() {
/*
    if(hdf_write_file == std::string("")) return;

    hdf_archive dump(myComm);
    if(!dump.create(hdf_write_file)) {
      app_error()<<" Error opening restart file in SparseGeneralHamiltonian. \n";
      return;
    }

    std::string path = "/Hamiltonian/SparseGeneralHamiltonian";
    if(dump.is_group( path )) {
      app_error()<<" ERROR: H5Group /Hamiltonian/SparseGeneralHamiltonian already exists in restart file. Not over-writing data in file. \n"; 
      return;
    }

    dump.push("Hamiltonian");
    dump.push("SparseGeneralHamiltonian");

    std::vector<int> Idata(7);
    Idata[0]=H1.size();
    Idata[1]=V2.size();
    Idata[2]=V2_2bar.size();
    Idata[3]=NMO;
    Idata[4]=NAEA;
    Idata[5]=NAEB;
    Idata[6]=spinRestricted?(0):(1);
    dump.write(Idata,"dims");
    
    Idata.resize(NAEA+NAEB);
    for(int i=0; i<NAEA; i++) Idata[i] = occup_alpha[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) Idata[i] = occup_beta[j];
    dump.write(Idata,"occups");

    std::vector<double> Rdata(2);
    Rdata[0] = NuclearCoulombEnergy;
    Rdata[1] = FrozenCoreEnergy; 
    dump.write(Rdata,"Energies");

    int sz = std::max( 2*H1.size(), std::max( 4*V2.size(), 4*V2_2bar.size() ) );
    std::vector<IndexType> ivec; 
    ivec.reserve(sz);

    ivec.resize(2*H1.size());
    for(int i=0, j=0; i<H1.size(); i++, j+=2) 
      std::tie (ivec[j],ivec[j+1],std::ignore) = H1[i];   
    dump.write(ivec,"H1_indx");

    ivec.clear();
    ivec.resize(4*V2.size());
    for(int i=0, j=0; i<V2.size(); i++, j+=4) 
      std::tie (ivec[j],ivec[j+1],ivec[j+2],ivec[j+3],std::ignore) = V2[i];   
    dump.write(ivec,"V2_indx");

    ivec.clear();
    ivec.resize(4*V2_2bar.size());
    for(int i=0, j=0; i<V2_2bar.size(); i++, j+=4) 
      std::tie (ivec[j],ivec[j+1],ivec[j+2],ivec[j+3],std::ignore) = V2_2bar[i];   
    dump.write(ivec,"V2_2bar_indx");

    std::vector<IndexType>().swap(ivec); 

    sz = std::max( H1.size(), std::max( V2.size(), V2_2bar.size() ) );
    std::vector<ValueType> vvec;    
    vvec.reserve(sz);    

    vvec.resize(H1.size());
    for(int i=0; i<H1.size(); i++)       
      std::tie (std::ignore,std::ignore,vvec[i]) = H1[i];
    dump.write(vvec,"H1");

    vvec.clear();
    vvec.resize(V2.size());
    for(int i=0; i<V2.size(); i++)
      std::tie (std::ignore,std::ignore,std::ignore,std::ignore,vvec[i]) = V2[i];
    dump.write(vvec,"V2");

    vvec.clear();
    vvec.resize(V2_2bar.size());
    for(int i=0; i<V2_2bar.size(); i++)
      std::tie (std::ignore,std::ignore,std::ignore,std::ignore,vvec[i]) = V2_2bar[i];
    dump.write(vvec,"V2_2bar");    

    std::vector<ValueType>().swap(vvec); 

    dump.pop();
    dump.pop();

    dump.flush();
    dump.close();
*/
  }

  void DDProjector::calculateHSPotentials_Diagonalization(ComplexSMSpMat& Spvn)
  {

    int rnk=0;
#if defined(USE_MPI)
    rnk = rank();
#endif

      int NMO2 = 2*NMO; 

      Timer.reset("Generic");       
      Timer.start("Generic");       

      ValueMatrix eigVec(NMO2);
      RealVector eigVal(NMO2);
      if(!DenseMatrixOperators::symEigenSysAll(NMO2,Pmat.data(),NMO2,eigVal.data(),eigVec.data(),NMO2) ) {
        app_error()<<"Problems with eigenvalue/eigenvector calculation in DDProjector::calculateHSPotentials_Diagonalization.\n";
        APP_ABORT("Problems with eigenvalue/eigenvector calculation in DDProjector::calculateHSPotentials_Diagonalization.\n");
      }

      Timer.stop("Generic");       
      if(rnk==0) app_log()<<" -- Time to solve eigenvalue problem: " <<Timer.average("Generic") <<"\n";

      for(int i=0; i<NMO2; i++) { 
        RealType scl = std::sqrt( std::max(0.0,eigVal[i]) )  ;
        for(int j=0; j<NMO2; j++)  
          eigVec(j,i) *= scl; 
      }        

      int cnt1=0;
      int cnt2=0;
      for(int i=0; i<NMO2; i++) { 
       if(eigVal[i] > std::abs(eigcut)) { 
         int cnt3=0;
         for(int j=0; j<NMO2; j++) 
            if(std::abs(eigVec(j,i)) > cutoff_sparse) cnt3++;
         if(cnt3 > 0) {
           cnt1++;
           cnt2 += cnt3;
         }
       }
      }

     // later on, instead of doing all ~M^4 terms, choose a few thousand randomly
     if(test_breakup) {

      if(rnk==0) app_log()<<" -- Testing Projector factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<2*NMO; i++) {
       for(IndexType j=0; j<2*NMO; j++) { 
           ValueType v2 = Pmat(i,j);
           ValueType v2c = 0.0;
           for(int n=0; n<NMO2; n++) v2c += eigVec(i,n)*myconj(eigVec(j,n));
           s+=std::abs(v2-v2c);
           if( max < std::abs(v2-v2c) ) max = std::abs(v2-v2c);
           if( std::abs(v2-v2c) > 10*eigcut ) {
             app_error()<<" Problems with Projector decomposition, i,j,P,Pc: "
                       <<i <<" "
                       <<j <<" "
                       <<v2 <<" "
                       <<v2c <<std::endl;
           }
         }
      }
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/eigcut/NMO/NMO/4.0 <<" " <<max <<" \n ********************************************\n"<<std::endl; 

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }

      Spvn.setDims(NMO2,cnt1);
      Spvn.allocate_serial(cnt2);

      ComplexType ifac = ComplexType(0.0,1.0);

      cnt1=0;
      for(int i=0; i<NMO2; i++) {
       if(eigVal[i] > std::abs(eigcut)) {
         int cnt3=0;
         for(int j=0; j<2*NMO; j++) {
           if(std::abs(ifac*eigVec(j,i)) > cutoff_sparse) {
             cnt3++;
             int jk = j*NMO+j;
             if(j>=NMO) jk-=NMO;
             Spvn.add(jk,cnt1,ifac*eigVec(j,i));
           }
         }
         if(cnt3 > 0) 
           cnt1++;
       }
      }

      app_log()<<"Number of HS potentials in DDProjector: " <<Spvn.cols() <<std::endl;
      app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;
      app_log()<<"Compressing Spvn. \n";
      Spvn.compress();
      app_log()<<"Done Compressing Spvn. \n";

  }

  void DDProjector::calculateHSPotentials(ComplexSMSpMat& Spvn)
  {
    //if(use_eig)
    calculateHSPotentials_Diagonalization(Spvn);
  }


  bool DDProjector::parse(xmlNodePtr cur)
  {

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur; 
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    std::string bkp("no"); 
    std::string use("no"); 
    ParameterSet m_param;
    m_param.add(eigcut,"cutoff_decomp","double");
    m_param.add(eigcut,"cutoff_decomposition","double");
    m_param.add(eigcut,"cutoff_factorization","double");
    m_param.add(eigcut,"cutoff_cholesky","double");
    m_param.add(cutoff_sparse,"cutoff_sparse","double");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(filename,"filename","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    m_param.add(bkp,"test_breakup","std::string");
    m_param.add(use,"useCholesky","std::string");

    std::string par("no");
    m_param.add(par,"paral_fac","std::string");
    m_param.add(par,"parallel_fac","std::string");
    m_param.add(par,"parallel_factorization","std::string");

    m_param.put(cur);

    std::transform(par.begin(),par.end(),par.begin(),(int (*)(int)) tolower);
    if(par == "yes" || par == "true") parallel_factorization = true;

    use_eig=true;
    std::transform(use.begin(),use.end(),use.begin(),(int (*)(int)) tolower);
    if(use == "true" || use == "yes") use_eig = false;
    std::transform(filetype.begin(),filetype.end(),filetype.begin(),(int (*)(int))tolower);
    std::transform(bkp.begin(),bkp.end(),bkp.begin(),(int (*)(int))tolower);
    if(bkp == "yes" || bkp == "true") test_breakup = true;  

    if(use_eig)
      app_log()<<"Calculating factorization of 2 body interaction with direct diagonalization.\n";
    else
      app_log()<<"Calculating factorization of 2 body interaction with Cholesky method.\n";

    std::transform(par.begin(),par.end(),par.begin(),(int (*)(int)) tolower);
    if(par == "yes" || par == "true") parallel_factorization = true;

    if(parallel_factorization)
      app_log()<<"Calculating factorization of 2-bofy hamiltonian in parallel. \n";
   
    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;
  } 

  // do a general check of parameters for consistency
  // make sure object was build consistently and correctly 
  bool DDProjector::checkObject() 
  {
    return true;
  }

}

