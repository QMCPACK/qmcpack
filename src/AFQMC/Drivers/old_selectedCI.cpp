#include<tuple>
#include<map>
#include<string>
#include<iomanip>

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include<Message/MPIObjectBase.h>
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"
#include <qmc_common.h>

#include "AFQMC/config.h"
#include "AFQMC/Drivers/selectedCI.h"

#include "AFQMC/Utilities/Utils.h"

namespace qmcplusplus {

bool selectedCI::run()
{

  app_log()<<" Running selectedCI. \n";

  int ne = NAEA+NAEB;
  std::vector<IndexType> intm;
  intm.reserve(ne*100000);
  std::vector<RealType> eigVal(1);
  std::vector<IndexType> vira,virb;
  std::vector<IndexType> indx(ne);
  vira.reserve(NMO-NAEA);
  virb.reserve(NMO-NAEA);
  ComplexMatrix eigVec;

  occ_orbs.reserve(ne*100000);
  ci.reserve(100000);
  for(int i=0; i<NAEA; i++)
    occ_orbs.push_back(i);
  for(int i=NMO; i<NMO+NAEB; i++)
    occ_orbs.push_back(i);
  ci.push_back(1.0);
  int nn = occ_orbs.size()/ne;
  for(int i=0; i<nn; i++)
    std::sort(occ_orbs.begin()+i*ne,occ_orbs.begin()+(i+1)*ne);

// ideas for speed-ups
// 1. Use the fact that H is stored in sparse ordered form to
//    only consider terms that have hmat*ci > cut   
//    To do this, make a routine taht given (i,j) it gets all
//    (k,l) such that h(i,j,k,l)*ci > cut     
// 2. Keep a different list for those determinants added in the last
//    step (those new to the list that haven't been "excited" before)
//    This is a reasonable approximation since the ci coeffs do not change
//    much. Then only excite from this list.  
// 3. For parallelization, we need to solve both hamiltonian storage (which will be bigger now)
//    and fare share of work. For now, assume hamiltonian is  
// 4. use (real) parallel? sparse eigenvalue solver that produces only a given # of states 

  for(int i=0; i<maxit; i++) {
    intm.clear();
    int nci = occ_orbs.size()/ne;
    int nterms=0;
    Timer.reset("Generic");
    Timer.start("Generic");
    std::vector<IndexType>::iterator it=occ_orbs.begin();
    for(int nc = 0; nc<nci; nc++, it+=ne) {
      vira.clear();
      virb.clear();
      for(int iv=0; iv<NMO; iv++)
        if(!binary_search (it, it+NAEA, iv)) vira.push_back(iv);
      for(int iv=NMO; iv<2*NMO; iv++)
        if(!binary_search (it+NAEA, it+ne, iv)) virb.push_back(iv);
      if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
      if(std::search(intm.begin(),intm.end(),it,it+ne)==intm.end()) {
        for(int ii=0; ii<ne; ii++) intm.push_back(*(it+ii));
        nterms++;
      }

      // double excitations 
      for(int ia=0; ia<NAEA; ia++) {
        int a = *(it+ia);
        // aa
        for(int ib=ia+1; ib<NAEA; ib++) {
          int b = *(it+ib);
          for(int ic=0; ic<NMO-NAEA; ic++) {
            int c = vira[ic];
            for(int id=ic+1; id<NMO-NAEA; id++) {
              int d = vira[id];
              ValueType Hij = sHam->H(a,b,c,d) - sHam->H(a,b,d,c);
              if(std::abs(Hij*ci[nc]) > cutoff_list) {
                if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                int sz = intm.size();
                indx.clear();
                for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                indx[ia] = c;
                indx[ib] = d;
                std::sort(indx.begin(),indx.end());
                if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                  for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                  nterms++;
                }
                if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                sz = intm.size();
                indx.clear();
                for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                indx[ia] = d;
                indx[ib] = c;
                std::sort(indx.begin(),indx.end());
                if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                  for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                  nterms++;
                }
              }
            }
          }
        }

        // ab
        for(int ib=NAEA; ib<ne; ib++) {
          int b = *(it+ib);
          for(int ic=0; ic<NMO-NAEA; ic++) {
           int c = vira[ic];
            for(int id=0; id<NMO-NAEB; id++) {
              int d = virb[id];
              ValueType Hij = sHam->H(a,b,c,d);
              if(std::abs(Hij*ci[nc]) > cutoff_list) {
                if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                int sz = intm.size();
                indx.clear();
                for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                indx[ia] = c;
                indx[ib] = d;
                std::sort(indx.begin(),indx.end());
                if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                  for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                  nterms++;
                }
              }
            }
          }
        }
      } 
      for(int ia=NAEA; ia<ne; ia++) {
        int a = *(it+ia);
        for(int ib=ia+1; ib<ne; ib++) {
          int b = *(it+ib);
          for(int ic=0; ic<NMO-NAEB; ic++) {
            int c = virb[ic];
            for(int id=ic+1; id<NMO-NAEB; id++) {
              int d = virb[id];
              ValueType Hij = sHam->H(a,b,c,d) - sHam->H(a,b,d,c);
              if(std::abs(Hij*ci[nc]) > cutoff_list) {
                if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                int sz = intm.size();
                indx.clear();
                for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                indx[ia] = c;
                indx[ib] = d;
                std::sort(indx.begin(),indx.end());
                if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                  for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                  nterms++;
                }
                if(nterms*ne == intm.capacity()) intm.reserve(ne*(nterms+10000));
                sz = intm.size();
                indx.clear();
                for(int ii=0; ii<ne; ii++) indx.push_back(*(it+ii));
                indx[ia] = d;
                indx[ib] = c;
                std::sort(indx.begin(),indx.end());
                if(std::search(intm.begin(),intm.end(),indx.begin(),indx.end())==intm.end()) {
                  for(int ii=0; ii<ne; ii++) intm.push_back(indx[ii]);
                  nterms++;
                }
              }  
            }  
          }    
        }
      }
    }  // states in occ_orbs
    Timer.stop("Generic");

    app_log()<<" Iteration: " <<i <<std::endl;
    app_log()<<" Intermediate list has " <<nterms <<" terms" <<std::endl;

    Timer.reset("Generic1");
    Timer.start("Generic1");

    bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,intm,nterms);
    if(!sucess) {
      app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
      return false;
    }
    app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic3") <<std::endl;
    app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic4") <<std::endl;
    occ_orbs.reserve(intm.size());
    occ_orbs.clear();
    ci.clear();
    ci.reserve(nterms);
    for(int ii=0,nt=0; ii<nterms; ii++) {
      //app_log()<<"ci " <<ii <<" " <<eigVec(0,ii) <<std::endl; 
      if(std::abs(eigVec(0,ii)) > cutoff_diag) {
        ci.push_back(eigVec(0,ii));
        for(int j=0; j<ne; j++)
         occ_orbs.push_back(intm[ii*ne+j]);
        nt++;
      }   
    }
    Timer.stop("Generic1");

    std::ofstream out("iterativeCI.dat");
    if(out.fail()) {
      app_error()<<" Problems opening iterativeCI.dat \n";
      return false;
    }
    {

      std::vector<std::tuple<double,int> > dets(ci.size());
      for(int i=0; i<ci.size(); i++) dets[i] = std::make_tuple( std::abs(ci[i]),i);
      std::sort( dets.begin(), dets.end(),
      [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
               {return (std::get<0>(a)>std::get<0>(b));} );
      out<<" &FCI \n NCI = " <<ci.size() <<" \n /\n";
      for(int i=0; i<ci.size(); i++) {
        out<<ci[std::get<1>(dets[i])] <<" ";
        for(int j=0; j<ne; j++) out<<occ_orbs[ std::get<1>(dets[i])*ne+j]+1 <<" ";
        out<<"\n";
      }
      out.close();
    }

    app_log()<<" Energy: " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
    app_log()<<" Number of determinants after truncation: " <<occ_orbs.size()/ne <<std::endl;
    app_log()<<" Timings: " <<Timer.total("Generic") <<" " <<Timer.total("Generic1") <<std::endl;

  } // iteration 

  if(diag_in_steps>0) {
    app_log()<<"\n***********************************************\n  #Determinants        Energy: " <<"\n";
    std::vector<std::tuple<double,int> > dets(ci.size());
    for(int i=0; i<ci.size(); i++) dets[i] = std::make_tuple( std::abs(ci[i]),i);
    std::sort( dets.begin(), dets.end(),
      [] (const std::tuple<double,int>& a, const std::tuple<double,int>& b)
               {return (std::get<0>(a)>std::get<0>(b));} );
    for(int i=1; i<ci.size(); i+=diag_in_steps) {
      intm.clear();
      for(int ki=0; ki<i; ki++) {
       int kk = std::get<1>(dets[ki]);
       for(int kj=0; kj<ne; kj++) intm.push_back(occ_orbs[ kk*ne+kj]);
      }
      bool sucess = diagonalizeTrialWavefunction(eigVal,eigVec,intm,i,false);
      if(!sucess) {
        app_error()<<" Error: Problems with diagonalizeTrialWavefunction. \n";
        return false;
      }
      app_log()<<i <<" " <<eigVal[0]+NuclearCoulombEnergy <<std::endl;
    }
    app_log()<<"***********************************************" <<std::endl <<std::endl;
  }

  return true; 
}

bool selectedCI::diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ComplexMatrix& eigVec, std::vector<IndexType>& occv, int nci, bool eigV )
{
  ComplexType one = ComplexType(1.0,0.0);
  ComplexType zero = ComplexType(0.0,0.0);
  bool sucess;

    for(int i=0; i<nci; i++)
      std::sort(occv.begin()+i*(NAEA+NAEB),occv.begin()+(i+1)*(NAEA+NAEB));

    if(myComm->rank()==0) {

      Timer.reset("Generic3");
      Timer.start("Generic3");
      ComplexMatrix hm(nci,nci);
      ComplexMatrix ov(nci,nci);
      ComplexType let;
      RealType sg;
      std::vector<IndexType> occ(NAEA+NAEB);
      IndexType n0,n1,n2,n3;
      std::vector<IndexType> DL(NAEA+NAEB);
      std::vector<IndexType> DR(NAEA+NAEB);

      // don't rely on H2_2bar in case it is removed
      for(int ki=0; ki<nci; ki++) {
        // i==j
        let=zero;
        std::vector<IndexType>::iterator it = occv.begin()+ki*(NAEA+NAEB);
        for(int i=0; i<NAEA+NAEB; i++)
        {
          let += sHam->H(*(it+i),*(it+i));
          for(int j=i+1; j<NAEA+NAEB; j++) {
            let += sHam->H(*(it+i),*(it+j),*(it+i),*(it+j)) - sHam->H(*(it+i),*(it+j),*(it+j),*(it+i));
          }
        }

        ov(ki,ki) = one;
        hm(ki,ki) = let;
        for(int kj=ki+1; kj<nci; kj++) {

          std::vector<IndexType>::iterator jt = occv.begin()+kj*(NAEA+NAEB);
          std::copy(it,it+NAEA+NAEB,DL.begin());
          std::copy(jt,jt+NAEA+NAEB,DR.begin());
          int cnt = cntExcitations(NAEA,NAEB,DL,DR,n0,n1,n2,n3,occ,sg);

          if(cnt==0) {
            app_error()<<" Error: Found repeated determinant in trial wave function in MultiPureSingleDeterminant \n";
            return false;
          } else if(cnt==2) {
            int nterms = NAEA+NAEB-1;
            let=sHam->H(n0,n1);
            for(int i=0; i<nterms; i++)
              let+=sHam->H(n0,occ[i],n1,occ[i]) - sHam->H(n0,occ[i],occ[i],n1);
            hm(ki,kj)=let*sg;
          } else if(cnt==4) {
            hm(ki,kj) = sg*(sHam->H(n0,n1,n2,n3) - sHam->H(n0,n1,n3,n2));
          } else {
            hm(ki,kj) =  zero;
          }

          ov(ki,kj) = ov(kj,ki) = zero;
          hm(kj,ki) = myconj(hm(ki,kj));
        }
      }
      Timer.stop("Generic3");
      //app_log()<<" Time to generate hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

      Timer.reset("Generic4");
      Timer.start("Generic4");

      eigVal.resize(1);
      eigVec.resize(1,nci);
      std::vector<int> ifail(nci);
      sucess = DenseMatrixOperators::genHermitianEigenSysSelect(nci,hm.data(),nci,ov.data(),nci,1,eigVal.data(),eigV,eigVec.data(),eigVec.size2(),ifail.data());

      Timer.stop("Generic4");
      //app_log()<<" Time to diagonalize hamiltonian in diagonalizeTrialWavefunction: " <<Timer.total("Generic2") <<std::endl;

    } else {
      eigVal.resize(1);
      eigVec.resize(1,nci);
    }

    myComm->bcast(sucess);
    myComm->bcast(eigVal.data(),eigVal.size(),0,myComm->getMPI());
    myComm->bcast(eigVec.data(),eigVec.size1()*eigVec.size2(),0,myComm->getMPI());

    return sucess;

}

bool selectedCI::parse(xmlNodePtr cur)
{
  if(cur==NULL) return false;

  if(cur == NULL)
    return false;

  xmlNodePtr curRoot=cur;

  maxit=0;
  cutoff_list=cutoff_diag=0;
  build_full_hamiltonian = true;

  std::string str("yes");
  ParameterSet m_param;
  m_param.add(str,"full_hamiltonian","std::string");
  m_param.add(output_filename,"output_filename","std::string");
  m_param.add(output_filename,"output","std::string");
  m_param.add(cutoff_list,"cutoff_list","double");
  m_param.add(cutoff_diag,"cutoff_diag","double");
  m_param.add(maxit,"maxit","int");
  m_param.add(diag_in_steps,"diag_steps","int");
  m_param.put(cur);

  std::transform(str.begin(),str.end(),str.begin(),(int (*)(int)) tolower);
  if(str == "no" || str == "false") build_full_hamiltonian = false; 

  return true;
}


bool selectedCI::setup(HamPtr h0, WSetPtr w0, PropPtr p0, WfnPtr wf0)
{

  if(h0==NULL) {
    app_error()<<" Error: Null Hamiltonian pointer in selectedCI::setup(). \n";
    return false; 
  }

  ham0=h0;
  wlkBucket=NULL;
  prop0=NULL;
  wfn0=NULL;

  sHam = dynamic_cast<SparseGeneralHamiltonian*>(ham0);
  if(!sHam) {
    app_error()<<" Error in MultiPureSingleDeterminant::getHamiltonian. \n"
               <<" Hamiltonian associated with MultiPureSingleDeterminant must of the \n"
               <<" type SparseGeneralHamiltonian. \n";
    APP_ABORT("");
  }

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Beginning Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  app_log()<<" Using " <<ncores_per_TG <<" cores per node in a TaskGroup. \n";
  // right now this TG is not used. It is needed for setup purposes and to  
  // get a unique TG number for every group of cores on a node (used in the WalkerSet)
  TG.setup(ncores_per_TG,1,false);  
  std::vector<int> TGdata(5);
  TG.getSetupInfo(TGdata); 
  CommBuffer.setup(TG.getCoreRank()==0,std::string("COMMBuffer_")+std::to_string(myComm->rank()),ncores_per_TG);
  TG.setBuffer(&CommBuffer);

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Hamiltonian \n"
           <<"****************************************************\n"
           <<std::endl;

  // hamiltonian
  if(!ham0->init(TGdata,&CommBuffer)) {
    app_error()<<"Error initializing Hamiltonian in selectedCI::setup" <<std::endl; 
    return false; 
  }   

  NuclearCoulombEnergy = static_cast<ValueType>(sHam->NuclearCoulombEnergy);

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Finished Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  return true;
}

// writes checkpoint file
bool selectedCI::checkpoint(int block, int step) 
{
  return true;
}

// sets up restart archive and reads  
bool selectedCI::restart(hdf_archive& read)
{
  return true;
}

bool selectedCI::clear()
{
  return true;
}

}

