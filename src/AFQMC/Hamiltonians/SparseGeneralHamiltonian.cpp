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
#include "Message/CommOperators.h"

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/exceptions.hpp>

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Hamiltonians/SparseGeneralHamiltonian.h"
#include "AFQMC/Utilities/readHeader.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"
#include "AFQMC/Utilities/Utils.h"

// use lambdas FIX FIX FIX
 struct my_pair_sorter {
  bool operator() (std::pair<qmcplusplus::RealType,int> a, std::pair<qmcplusplus::RealType,int> b) {return a.first<b.first;}
  bool operator() (std::pair<std::complex<qmcplusplus::RealType>,int> a, std::pair<std::complex<qmcplusplus::RealType>,int> b) {return a.first.real()<b.first.real();}
} myobject;

namespace qmcplusplus
{


  bool SparseGeneralHamiltonian::readFCIDUMP(const std::string& fileName, bool minimizeIO) //,
//      std::vector<s2D<ValueType> >& hij, std::vector<s4D<ValueType> >& Vijkl)
  {

    int rnk=0;
#if defined(USE_MPI)
    rnk = rank();
#endif

    std::vector<s2D<ValueType> > hij_core;
    std::vector<s4D<ValueType> > Vijkl_core, Vijkl_mixed;
    ValueSMSpMat V2_fact_c, V2_fact_m;
    
     std::ifstream in;
#ifdef _BGQ_     
     // trying to speed up I/O in BGQ
     // this might cause a small memory leak, fix later
     //int bsize = 100000;   
     //char *mybuffer = new char[bsize]; 
     //in.rdbuf()->pubsetbuf(mybuffer,bsize);
#endif
     in.open(fileName.c_str());
     if(in.fail()) {
        app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
        return false;
     }

     if(!readHeader(in,NMAX,NMO_FULL,NETOT,NAEA,NAEB,NCA,NCB,MS2,spinRestricted,ISYM,occup_alpha,occup_beta,orbSymm,occupPerSymm_alpha,occupPerSymm_beta,orderStates,factorizedHamiltonian)) {
       app_error()<<" Error: Problem with header section of file. \n";
       return false;
     }

     NMO = NMO_FULL-NCA;

     // no hamiltonian distribution yet
     min_i = 0;
     max_i = NMO;
     if(!spinRestricted) max_i *= 2;  // or 4?

     if(orderStates && NMAX != NMO_FULL) {
       app_error()<<"Error in SparseGeneralHamiltonian::readFCIDUMP. orderStates can not be used in combination with NMAX!=NMO_FULL. \n ";
       return false;
     }

     if(NCA != NCB) {
       app_error()<<"Error in SparseGeneralHamiltonian::readFCIDUMP. NCA!=NCB. Not sure how to implement this! \n";
       return false;
     }

     // MMORALES: It is possible to define variables (e.g. NMO_FULL, NAEA, ...) in the input xml file and these might be different from those in this file. But after the file header is read, all variables must be in an initialized state  
     if(!checkAFQMCInfoState()) {
       app_error()<<" ERROR: Problem with basic parameters during Hamiltonian initialization. \n"
                  <<"        State of Hamiltonian::AFQMCInfo after reading header section of ASCII integral file. \n" <<std::endl;
       printAFQMCInfoState(app_error());  
       return false;
     } 
     
     if(!spinRestricted) app_log()<<"Found UHF orbitals in Molpro integral file. Running Spin-Unrestricted calculation. \n";
     else 
       app_log()<<"Found RHF/ROHF orbitals in Molpro integral file. Running Spin-Restricted calculation. \n";

     if(occup_alpha.size() == 0 && occup_beta.size() == 0 ) {
       app_log()<<" WARNING: OCCUP tag not found on integral file. Assuming ground-state occupation. \n";
       occup_alpha.resize(NAEA);
       occup_beta.resize(NAEB);
       for(int i=0; i<NAEA; i++) occup_alpha[i]=i;
       for(int i=0; i<NAEB; i++) occup_beta[i]=i+NMO;
     }

     if(occup_alpha.size() != NAEA) {
       app_error() <<"Error: size of OCCUP_ALPHA must be equal to NAEA. \n" <<std::endl;
       return false;
     }
     if(occup_beta.size() != NAEB) {
       app_error() <<"Error: size of OCCUP_BETA must be equal to NAEB. \n" <<std::endl;
       return false;
     }


     std::streampos start = in.tellg();
     int nOne_core,nOne,nTwo,nTwo_core,nTwo_mixed;
     int nThree,nThree_core,nThree_mixed;

     std::map<IndexType,IndexType> orbMapA; 
     std::map<IndexType,IndexType> orbMapB; 
     orbMapA[0]=0;
     orbMapB[0]=0;
     for(int i=1; i<=NMO_FULL; i++) orbMapA[i]=i;   
     for(int i=1; i<=NMO_FULL; i++) orbMapB[i]=i+NMO_FULL;   

     // setup SM objects 
     V2.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2"),TG.getNodeCommLocal());
     V2_2bar.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2_2bar"),TG.getNodeCommLocal());
     V2_full.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2_full"),TG.getNodeCommLocal());

     V2_fact.setup(head_of_nodes,"V2_fact",TG.getNodeCommLocal());
     V2_fact_c.setup(head_of_nodes,"V2_fact_c",TG.getNodeCommLocal());
     V2_fact_m.setup(head_of_nodes,"V2_fact_m",TG.getNodeCommLocal());
     int n3Vecs=0;

     if(orderStates) { 
     // hacking these routines a bit to allow FC rotation 
    
     // app_error()<<" Error: reordering of states has benn temporarily disabled!!! " <<std::endl;
     // return false;
 
      int NCA_tmp=NCA, NCB_tmp=NCB;
      NCA=NCB=0;
      NMO=NMO_FULL;

    app_log()<<" Free memory before count() : " <<freemem() <<" MB. \n";
    
      if(!countElementsFromFCIDUMP(in,nOne,nOne_core,nTwo,nTwo_core,nTwo_mixed,nThree,nThree_core,nThree_mixed,orbMapA,orbMapB,n3Vecs)) {
        app_error()<<"Error in readFCIDUMP. Problem counting elements. \n" <<std::endl;
        return false;
      }

    app_log()<<" Free memory after count() : " <<freemem() <<" MB. \n";

      if(nThree > 0)  {
        if(nTwo > 0) {
          app_error()<<"Found both 3-Index and 4-Index terms in FCIDUMP. Only one form is allowed. " <<std::endl; 
          return false;
        }
        factorizedHamiltonian = true;
      } else
        factorizedHamiltonian = false;

      if(nOne_core > 0 || nTwo_core > 0 || nTwo_mixed > 0 || nThree_core > 0 || nThree_mixed > 0) {
        app_error()<<"Error in readFCIDUMP. Problem counting elements. #core > 0 in orderStates. \n" <<std::endl;
        return false;
      }
      if(factorizedHamiltonian) {
        int NMO2 = NMO*NMO;
        if(!spinRestricted) NMO2 *= 2; 
        V2_fact.setDims(NMO2,n3Vecs);
        V2_fact_c.setDims(NMO2,n3Vecs);
        V2_fact_m.setDims(NMO2,n3Vecs);
        V2_fact.reserve(nThree);    // V2_fact must be reserved, since elements are added with "add" 
      } else {
        if(n3Vecs > 0) {
          app_error()<<"Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!" <<std::endl;
          APP_ABORT("Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!");
        }
        V2.resize(nTwo);
      }

      app_log()<<" Number of one- and two- electron integrals read from file (assuming no core states for rotation)  : " <<nOne <<" " <<nTwo <<std::endl;

      in.clear();
      in.seekg(start);

      Timer.reset("Generic");
      Timer.start("Generic");

      H1.resize(nOne); 
      hij_core.resize(nOne_core); 

      Vijkl_core.resize(nTwo_core);  
      Vijkl_mixed.resize(nTwo_mixed);  

      if(!readElementsFromFCIDUMP(in,H1,hij_core,V2,Vijkl_core,Vijkl_mixed,V2_fact,V2_fact_c,V2_fact_m,orbMapA,orbMapB)) {
        app_error()<<"Error in return from readElementsFromFCIDUMP in readFCIDUMP.\n" <<std::endl; 
        return false;
      }

      if(factorizedHamiltonian)
        if(H1.size() != nOne || hij_core.size() != nOne_core) {
          app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl;
          return false;
        }
      else
        if(H1.size() != nOne || hij_core.size() != nOne_core  || V2.size()!=nTwo ||  Vijkl_core.size()!=nTwo_core ||  Vijkl_mixed.size()!=nTwo_mixed) {
          app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl;
          return false;
        }

      Timer.stop("Generic");
      if(rnk==0) app_log()<<" -- Time to read ASCII integral file: " <<Timer.average("Generic") <<"\n";


      if(occupPerSymm_alpha.size() != occupPerSymm_beta.size()) {
        app_error()<<" Error: Different number of symmetry groups between alpha/beta. \n";
        return false; 
      } 
      int nSym=occupPerSymm_alpha.size();
      std::vector< std::vector<int> > OrbsPerSymm(nSym);
      for(int i=0; i<orbSymm.size(); i++) 
        OrbsPerSymm[orbSymm[i]-1].push_back(i);

      // sort integral tables
      std::sort (H1.begin(), H1.end(),mySort);
      if(factorizedHamiltonian)
        V2_fact.compress();
      else
        if(head_of_nodes) std::sort (V2.begin(), V2.end(),mySort);
      myComm->barrier();
    
      orbMapA.clear(); 
      orbMapA[0]=0;
      orbMapB.clear(); 
      orbMapB[0]=0;

      if(spinRestricted) {
        std::vector< std::pair<ValueType,int> > orbEnergies(NMO_FULL);

        for(IndexType i=0; i<NMO_FULL; i++) {
          orbEnergies[i].second = i+1; 
          orbEnergies[i].first = H(i,i); 
          // loop through symmetry groups
          for(int j=0; j<occupPerSymm_alpha.size(); j++) {
            // loop through all orbitals in symmetry group j
            for(int k=0; k<occupPerSymm_alpha[j]; k++) {
              int b = OrbsPerSymm[j][k];
              if(i != b) orbEnergies[i].first += (H(i,b,i,b)-H(i,b,b,i));
            }
          }
          for(int j=0; j<occupPerSymm_beta.size(); j++) {
            for(int k=0; k<occupPerSymm_beta[j]; k++) {
              int b = OrbsPerSymm[j][k]+NMO_FULL;
              orbEnergies[i].first += H(i,b,i,b);
            }
          }
        }

        std::sort(orbEnergies.begin(),orbEnergies.end(),myobject);
        for(int i=0; i<orbEnergies.size(); i++) orbMapA[orbEnergies[i].second] = i+1;
        for(int i=0; i<orbEnergies.size(); i++) orbMapB[orbEnergies[i].second] = i+1+NMO_FULL;

        if(rnk==0) {
          app_log()<<"Re-ordering orbitals to order them by energy. \n"
              <<"New Index      Old Index      Orbital Energy: " <<"\n";
          for(int i=0; i<NMO_FULL; i++) 
            app_log()<<i+1     <<"    "    <<orbEnergies[i].second    <<"     "   <<orbEnergies[i].first <<"\n"; 
        }

      } else {

        std::vector< std::pair<ValueType,int> > orbEnergies(NMO_FULL*2);
        for(int i=0; i<NMO_FULL; i++) {
          orbEnergies[i].second = i+1;
          orbEnergies[i].first = H(i,i); 
          orbEnergies[i+NMO_FULL].second = i+1;
          orbEnergies[i+NMO_FULL].first = H(NMO_FULL+i,NMO_FULL+i); 
          // loop through symmetry groups
          for(int j=0; j<occupPerSymm_alpha.size(); j++) {
            for(int k=0; k<occupPerSymm_alpha[j]; k++) {
              int b = OrbsPerSymm[j][k];
              if(i != b) orbEnergies[i].first += (H(i,b,i,b)-H(i,b,b,i));
              orbEnergies[i+NMO_FULL].first += (H(i+NMO_FULL,b,i+NMO_FULL,b)-H(i+NMO_FULL,b,b,i+NMO_FULL));
            }
          }
          for(int j=0; j<occupPerSymm_beta.size(); j++) {
            for(int k=0; k<occupPerSymm_beta[j]; k++) {
              int b = OrbsPerSymm[j][k]+NMO_FULL;
              orbEnergies[i].first += (H(i,b,i,b)-H(i,b,b,i));
              if(i+NMO_FULL != b) orbEnergies[i+NMO_FULL].first += (H(i+NMO_FULL,b,i+NMO_FULL,b)-H(i+NMO_FULL,b,b,i+NMO_FULL));
            }
          }
        }

        std::sort(orbEnergies.begin(),orbEnergies.begin()+NMO_FULL,myobject);
        for(int i=0; i<NMO_FULL; i++) orbMapA[orbEnergies[i].second] = i+1;
        std::sort(orbEnergies.begin()+NMO_FULL,orbEnergies.end(),myobject);
        for(int i=0; i<NMO_FULL; i++) orbMapB[orbEnergies[i+NMO_FULL].second] = i+1+NMO_FULL;

        if(rnk==0) {
          app_log()<<"Re-ordering orbitals to order them by energy. \n"
              <<"New Index      Old Index      Orbital Energy   : " <<"\n";
          for(int i=0; i<NMO_FULL; i++) {
            app_log()<<i+1     <<" alpha   "    <<orbEnergies[i].second    <<"     "   <<orbEnergies[i].first <<"\n";
            app_log()<<"  "     <<" beta   "    <<orbEnergies[i+NMO_FULL].second    <<"     "   <<orbEnergies[i+NMO_FULL].first <<"\n";
          }
        }

      } 

       in.clear();
       in.seekg(start);

       // restore values
       NCA=NCA_tmp;
       NCB=NCB_tmp;
       NMO=NMO_FULL-NCA;

       V2_fact.clear();

       myComm->barrier();
     } 

     if(!countElementsFromFCIDUMP(in,nOne,nOne_core,nTwo,nTwo_core,nTwo_mixed,nThree,nThree_core,nThree_mixed,orbMapA,orbMapB,n3Vecs)) {
       app_error()<<"Error in readFCIDUMP. Problem counting elements. \n" <<std::endl;
       return false;
     }

     // need to fix the problem wirh H(i,j,k,l,Vijkl) with factorized hamiltonians
     // right now it doesn't use the correct data structure. Need a new routine for V2_fact_mixed 
     if(nThree_core > 0) {
       app_error()<<" Error: Core states have been temporarily disabled with factorized hamiltonians." <<std::endl;
       return false;
     }

     H1.resize(nOne);
     hij_core.resize(nOne_core);

     if(nThree > 0)  {
       if(nTwo > 0) {
         app_error()<<"Found both 3-Index and 4-Index terms in FCIDUMP. Only one form is allowed. " <<std::endl;
         return false;
       }
       app_log()<<"Found integral file with Cholesky decomposed integrals.\n";
       app_log()<<"# Chol Vec, # terms:" <<n3Vecs <<" " <<nThree <<std::endl;
       factorizedHamiltonian = true;
       int NMO2 = NMO*NMO;
       if(!spinRestricted) NMO2 *= 2; 
       V2_fact.setDims(NMO2,n3Vecs);
       V2_fact_c.setDims(NMO2,n3Vecs);
       V2_fact_m.setDims(NMO2,n3Vecs);
       // V2_fact must be reserved, since elements are added with "add"
       V2_fact.reserve(nThree);
       V2_fact_c.reserve(nThree_core);
       V2_fact_m.reserve(nThree_mixed);
     } else {
       if(n3Vecs > 0) {
         app_error()<<"Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!" <<std::endl; 
         APP_ABORT("Found three index terms in FCIDUMP. (Only allowed with factorized hamiltonian. Check!!!"); 
       } 
       factorizedHamiltonian = false;

       V2.resize(nTwo,true);  // allow reduction of size if necessary
       Vijkl_core.resize(nTwo_core);
       Vijkl_mixed.resize(nTwo_mixed);
     }
 
     in.clear();
     in.seekg(start);

     Timer.reset("Generic");
     Timer.start("Generic");

     if(!readElementsFromFCIDUMP(in,H1,hij_core,V2,Vijkl_core,Vijkl_mixed,V2_fact,V2_fact_c,V2_fact_m,orbMapA,orbMapB)) {
       app_error()<<"Error in return from readElementsFromFCIDUMP in readFCIDUMP.\n" <<std::endl; 
       return false;
     }

     if(factorizedHamiltonian) {
       if(H1.size() != nOne || hij_core.size() != nOne_core ) {
         app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl; 
         return false;
       }
     } else {
       if(H1.size() != nOne || hij_core.size() != nOne_core  || V2.size()!=nTwo ||  Vijkl_core.size()!=nTwo_core ||  Vijkl_mixed.size()!=nTwo_mixed) {
         app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl; 
         return false;
       }
     }
     myComm->barrier(); 

     Timer.stop("Generic");
     if(rnk==0) app_log()<<" -- Time to read ASCII integral file (2nd time after reorder): " <<Timer.average("Generic") <<"\n";

//     delete mybuffer; 
     in.close();

     if(minimizeIO) {

       //MPI_Bcast (&(oneE[0]), nq, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

       //MPI_Bcast (&(twoEaa[0]), nq, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

       //MPI_Bcast (&(twoEbb[0]), nq, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

     }

#ifdef AFQMC_DEBUG
    app_log()<<" Finished reading FCIDUMP file." <<std::endl; 
    app_log()<<" Sorting integral tables." <<std::endl; 
#endif

    Timer.reset("Generic");
    Timer.start("Generic");

    // sort integral tables
    // remove repeated, needed in the (beta/alpha) and (alpha/beta) blocks
    std::sort (H1.begin(), H1.end(),mySort);
    s2Dit ith = std::unique(H1.begin(),H1.end(),myEqv);
    H1.resize( std::distance(H1.begin(),ith) );

    std::sort (hij_core.begin(), hij_core.end(),mySort);
    ith = std::unique(hij_core.begin(),hij_core.end(),myEqv);
    hij_core.resize( std::distance(hij_core.begin(),ith) );

    if(factorizedHamiltonian) {
      V2_fact.compress();
      V2_fact_m.compress();
      V2_fact_c.compress();
      myComm->barrier();
    } else {
      if(head_of_nodes) {
        std::sort (V2.begin(), V2.end(),mySort);
        s4Dit itV = std::unique(V2.begin(),V2.end(),myEqv);
        // keeping extra space because array might be quite large to reallocate 
        V2.resize_serial( std::distance(V2.begin(),itV));
      }

      std::sort (Vijkl_core.begin(), Vijkl_core.end(),mySort);
      std::vector<s4D<ValueType>>::iterator itV = std::unique(Vijkl_core.begin(),Vijkl_core.end(),myEqv);
      Vijkl_core.resize( std::distance(Vijkl_core.begin(),itV) );

      std::sort (Vijkl_mixed.begin(), Vijkl_mixed.end(),mySort);
      itV = std::unique(Vijkl_mixed.begin(),Vijkl_mixed.end(),myEqv);
      Vijkl_mixed.resize( std::distance(Vijkl_mixed.begin(),itV) );
    }

    Timer.stop("Generic");
    if(rnk==0) app_log()<<" -- Time to sort sparse integral tables: " <<Timer.average("Generic") <<"\n";

    // calculate list of virtuals 
    isOcc_alpha.clear();
    isOcc_beta.clear();
    for(IndexType i=0; i<2*NMO; i++) isOcc_alpha[i]=false;
    for(IndexType i=0; i<2*NMO; i++) isOcc_beta[i]=false;
    for(int i=0; i<occup_alpha.size(); i++) isOcc_alpha[ occup_alpha[i] ]=true;
    for(int i=0; i<occup_beta.size(); i++) isOcc_beta[ occup_beta[i] ]=true;
    virtual_alpha.clear();
    virtual_beta.clear();
    virtual_alpha.reserve(NMO-NAEA);
    virtual_beta.reserve(NMO-NAEB);
    for(IndexType i=0; i<NMO; i++)
      if(!isOcc_alpha[i]) virtual_alpha.push_back(i);
    for(IndexType i=NMO; i<2*NMO; i++)
      if(!isOcc_beta[i]) virtual_beta.push_back(i);
    
    // define close_shell
    if(NCA==NCB && NAEA == NAEB && occup_alpha.size() == occup_beta.size() ) { 
      close_shell = true;
      for(int i=0; i<occup_alpha.size(); i++) 
        if( occup_alpha[i] != occup_beta[i]-NMO )
          close_shell = false;
    } else
      close_shell = false;

#ifdef AFQMC_DEBUG
    app_log()<<" Finished sorting integral tables." <<std::endl; 
#endif

    // FZ core rotation
    if(NCA > 0 || NCB > 0) {
#ifdef AFQMC_DEBUG
    app_log()<<" Performing Frozen-Core manipulations. \n"; 
#endif   

    Timer.reset("Generic");
    Timer.start("Generic");

    FrozenCoreEnergy = ValueType(0.0);

    for(int i=0; i<hij_core.size(); i++) 
      if( std::get<0>( hij_core[i] ) == std::get<1>( hij_core[i] )  ) {
        FrozenCoreEnergy += std::get<2>( hij_core[i] );
        if(spinRestricted) FrozenCoreEnergy += std::get<2>( hij_core[i] );
      }
    app_log()<<" One-Body Core Energy: " <<FrozenCoreEnergy <<std::endl;     

    for(int i=0; i<NCA; i++)
    {
      for(int j=i+1; j<NCA; j++)
        FrozenCoreEnergy += (H(i,j,i,j,Vijkl_core,NMO_FULL) - H(i,j,j,i,Vijkl_core,NMO_FULL) );
      for(int j=NMO_FULL; j<NMO_FULL+NCB; j++)
        FrozenCoreEnergy += H(i,j,i,j,Vijkl_core,NMO_FULL);
    }
    for(int i=NMO_FULL; i<NMO_FULL+NCB; i++)
    {
      for(int j=i+1; j<NMO_FULL+NCB; j++)
        FrozenCoreEnergy += (H(i,j,i,j,Vijkl_core,NMO_FULL) - H(i,j,j,i,Vijkl_core,NMO_FULL));
    }    
    app_log()<<" Frozen core energy: " <<FrozenCoreEnergy <<std::endl;
    NuclearCoulombEnergy += FrozenCoreEnergy;

    // lazy for now
    std::vector<s2D<ValueType> > h1_new;
    int cnt=0;
    for(int i=0; i<NMO; i++) {
      for(int j=i; j<NMO; j++) {
        ValueType V = H(i,j); 
        for(int k=0; k<NCA; k++) 
          V += H(k,NCA+i,k,NCA+j,Vijkl_mixed,NMO_FULL) - H(k,NCA+i,NCA+j,k,Vijkl_mixed,NMO_FULL); 
        for(int k=NMO_FULL; k<NMO_FULL+NCB; k++) 
          V += H(k,NCA+i,k,NCA+j,Vijkl_mixed,NMO_FULL);
        if(std::abs(V) > cutoff1bar ) cnt++;
      }
    }
    if(!spinRestricted) {
      for(int i=NMO; i<2*NMO; i++) {
        for(int j=i; j<2*NMO; j++) {
          ValueType V = H(i,j);
          for(int k=0; k<NCA; k++) 
            V += H(k,i-NMO+NMO_FULL+NCB,k,j-NMO+NMO_FULL+NCB,Vijkl_mixed,NMO_FULL);
          for(int k=NMO_FULL; k<NMO_FULL+NCB; k++) 
            V += H(k,i-NMO+NMO_FULL+NCB,k,j-NMO+NMO_FULL+NCB,Vijkl_mixed,NMO_FULL) - H(k,i-NMO+NMO_FULL+NCB,j-NMO+NMO_FULL+NCB,k,Vijkl_mixed,NMO_FULL);
          if(std::abs(V) > cutoff1bar ) cnt++;
        }
      }
    }
    h1_new.reserve(cnt);
    for(int i=0; i<NMO; i++) {
      for(int j=i; j<NMO; j++) {
        ValueType V = H(i,j);
        for(int k=0; k<NCA; k++) 
          V += H(k,NCA+i,k,NCA+j,Vijkl_mixed,NMO_FULL) - H(k,NCA+i,NCA+j,k,Vijkl_mixed,NMO_FULL);
        for(int k=NMO_FULL; k<NMO_FULL+NCB; k++)
          V += H(k,NCA+i,k,NCA+j,Vijkl_mixed,NMO_FULL);
        if(std::abs(V) > cutoff1bar ) h1_new.push_back(std::forward_as_tuple(i,j,V));
      }
    }
    if(!spinRestricted) {
      for(int i=NMO; i<2*NMO; i++) {
        for(int j=i; j<2*NMO; j++) {
          ValueType V = H(i,j);
          for(int k=0; k<NCA; k++) 
            V += H(k,i-NMO+NMO_FULL+NCB,k,j-NMO+NMO_FULL+NCB,Vijkl_mixed,NMO_FULL);
          for(int k=NMO_FULL; k<NMO_FULL+NCB; k++)
            V += H(k,i-NMO+NMO_FULL+NCB,k,j-NMO+NMO_FULL+NCB,Vijkl_mixed,NMO_FULL) - H(k,i-NMO+NMO_FULL+NCB,j-NMO+NMO_FULL+NCB,k,Vijkl_mixed,NMO_FULL);
        if(std::abs(V) > cutoff1bar ) h1_new.push_back(std::forward_as_tuple(i,j,V));
        }
      }
    }  

    H1 = h1_new;
    std::sort (H1.begin(), H1.end(),mySort);

    Timer.stop("Generic");
    if(rnk==0) app_log()<<" -- Time to perform FC manipulations: " <<Timer.average("Generic") <<"\n";

#ifdef AFQMC_DEBUG
    app_log()<<" Finished performing Frozen-Core manipulations. \n"; 
#endif    
    }

    app_log()<<" Calculating eigenvalues: " <<std::endl; 

    ValueType ekin=0.0,epot_bb=0.0,epot_ab=0.0,epot_aa=0.0,epot_coul=0.0,epot=NuclearCoulombEnergy,Emp2=0.0;
    for(std::vector<IndexType>::iterator ita = occup_alpha.begin(); ita!=occup_alpha.end(); ita++)
    {
      ekin += H(*ita,*ita);
      for(std::vector<IndexType>::iterator itb = ita+1; itb<occup_alpha.end(); itb++) 
        epot += H_2bar(*ita,*itb,*ita,*itb);
      for(std::vector<IndexType>::iterator itb = occup_beta.begin(); itb!=occup_beta.end(); itb++) 
        epot += H(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_beta.begin(); ita!=occup_beta.end(); ita++)
    {
      ekin += H(*ita,*ita);
      for(std::vector<IndexType>::iterator itb = ita+1; itb!=occup_beta.end(); itb++) 
        epot += H_2bar(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_alpha.begin(); ita!=occup_alpha.end(); ita++)
    {
      for(std::vector<IndexType>::iterator itb = ita+1; itb<occup_alpha.end(); itb++)
        epot_coul += H(*ita,*itb,*ita,*itb);
      for(std::vector<IndexType>::iterator itb = occup_beta.begin(); itb!=occup_beta.end(); itb++)
        epot_coul += H(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_beta.begin(); ita!=occup_beta.end(); ita++)
    {
      for(std::vector<IndexType>::iterator itb = ita+1; itb!=occup_beta.end(); itb++)
        epot_coul += H(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_alpha.begin(); ita!=occup_alpha.end(); ita++)
    {
      for(std::vector<IndexType>::iterator itb = ita+1; itb<occup_alpha.end(); itb++)
        epot_aa += H_2bar(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_beta.begin(); ita!=occup_beta.end(); ita++)
    {
      for(std::vector<IndexType>::iterator itb = ita+1; itb!=occup_beta.end(); itb++)
        epot_bb += H_2bar(*ita,*itb,*ita,*itb);
    }
    for(std::vector<IndexType>::iterator ita = occup_alpha.begin(); ita!=occup_alpha.end(); ita++)
    {
      for(std::vector<IndexType>::iterator itb = occup_beta.begin(); itb!=occup_beta.end(); itb++)
        epot_ab += H(*ita,*itb,*ita,*itb);
    }

    // for spinRestricted, only alpha sector is populated
    app_log()<<"Eigenvalues: " <<std::endl;   
    eig.resize(2*NMO);
    for(IndexType i=0; i<NMO; i++) {
      eig[i] = H(i,i);
      for(std::vector<IndexType>::iterator it = occup_alpha.begin(); it<occup_alpha.end(); it++) 
        eig[i] += H_2bar(i,*it,i,*it);
      for(std::vector<IndexType>::iterator it = occup_beta.begin(); it<occup_beta.end(); it++) 
        eig[i] += H(i,*it,i,*it);
      app_log()<<i <<"   " <<eig[i] <<std::endl;
    }    
    if(!spinRestricted) {
      for(IndexType i=NMO; i<2*NMO; i++) {
        eig[i] = H(i,i);
        for(std::vector<IndexType>::iterator it = occup_alpha.begin(); it<occup_alpha.end(); it++) 
          eig[i] += H(i,*it,i,*it);
        for(std::vector<IndexType>::iterator it = occup_beta.begin(); it<occup_beta.end(); it++) 
          eig[i] += H_2bar(i,*it,i,*it);
        app_log()<<i <<"   " <<eig[i] <<std::endl;
      }    
    }

    bool printFockMatrix=false; 
    if(printFockMatrix) {
      app_log()<<" Off-diagonal elements of the Fock matrix larger than " <<5*std::max(cutoff2bar,cutoff1bar) <<std::endl;
      for(IndexType i=0; i<NMO; i++) {
       for(IndexType j=i+1; j<NMO; j++) {
        ValueType fij = H(i,j);
        for(std::vector<IndexType>::iterator it = occup_alpha.begin(); it<occup_alpha.end(); it++)
          fij += H_2bar(i,*it,j,*it);
        for(std::vector<IndexType>::iterator it = occup_beta.begin(); it<occup_beta.end(); it++) 
          fij += H(i,*it,j,*it);
        if( std::abs(fij) > 5*std::max(cutoff2bar,cutoff1bar) )
          app_log()<<i <<"  " <<j <<"     " <<fij <<std::endl;
       }
      }
      if(!spinRestricted) {
      }
      app_log()<<" Done printing Fock matrix." <<std::endl;
    } 

    app_log()<<std::setprecision(12) <<std::endl <<" Ekin:     " <<ekin <<std::endl;
    app_log()<<" Epot:     " <<epot <<std::endl;
    app_log()<<" Epot_Coulomb:     " <<epot_coul <<std::endl;
    app_log()<<" Epot_aa:     " <<epot_aa <<std::endl;
    app_log()<<" Epot_bb:     " <<epot_bb <<std::endl;
    app_log()<<" Epot_ab:     " <<epot_ab <<std::endl;
    app_log()<<" Ehf:     " <<epot+ekin <<std::endl;
    if(spinRestricted && NCA == NCB && NAEA==NAEB && close_shell) {
      // right now only for close-shell RHF 
      // using "trick" to get 1 bar terms
      std::vector<IndexType>::iterator iti,itj,ita,itb; 
/*
      for(iti=occup_alpha.begin(); iti<occup_alpha.end(); iti++)
      for(itj=occup_alpha.begin(); itj<occup_alpha.end(); itj++)
      for(ita=virtual_alpha.begin(); ita<virtual_alpha.end(); ita++)
      for(itb=virtual_alpha.begin(); itb<virtual_alpha.end(); itb++)
        Emp2 += H(*iti,*itj,*ita,*itb)*(2.0*H(*ita,*itb,*iti,*itj) - H(*ita,*itb,*itj,*iti))/(eig[*iti]+eig[*itj]-eig[*ita]-eig[*itb]); 
*/
      app_log()<<" Emp2:     " <<Emp2 <<std::endl;
      app_log()<<" Ehf+Emp2: " <<epot+ekin+Emp2 <<std::endl <<std::endl;
    } 

    app_log()<<" Memory used by 2-el integral table: " <<sizeof(s4D<ValueType>)*V2.size()/1024/1024 <<" MB. " <<std::endl;

    if(rnk == 0) {
      hdf_write();
      ascii_write();
    }

    return true;
  }

  bool SparseGeneralHamiltonian::generate_V2_2bar() {

    if(V2_2bar.size() != 0) return true;

    if(factorizedHamiltonian) {
      app_error()<<" Error in SparseGeneralHamiltonian::generate_V2_2bar will not generate 2bar integrals with factorized hamiltonian. \n\n\n";
      return false;
    }

#ifdef AFQMC_DEBUG
    app_log()<<" Generating sparse 2-bar integrals. " <<std::endl; 
#endif

    int cnt=0;
    Timer.reset("Generic");
    Timer.start("Generic");


    // using dumbest algorithm right now, improve later
    if(head_of_nodes) {

      cnt=0; 
      for(IndexType i=0; i<NMO; i++)  {  
      for(IndexType j=i+1; j<NMO; j++)  {  
      for(IndexType k=i; k<NMO; k++)  {  
      for(IndexType l=k+1; l<NMO; l++)  {  
        if( i==j || k==l ) continue;
        s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0)); 
        s4D<ValueType> s = find_smaller_equivalent_TwoBar_for_integral_list(ijkl);
        if(ijkl == s) { // found an element, now see if it is non-zero
          ValueType v = H(i,j,k,l)-H(i,j,l,k);
          if( std::abs(v) > cutoff2bar )
            cnt++;  
        }
      }
      }
      }
      }
      if(!spinRestricted) {
       for(IndexType i=NMO; i<2*NMO; i++)  {  
       for(IndexType j=i+1; j<2*NMO; j++)  {  
       for(IndexType k=i; k<2*NMO; k++)  {  
       for(IndexType l=k+1; l<2*NMO; l++)  {  
        if( i==j || k==l ) continue;
         s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0)); 
         s4D<ValueType> s = find_smaller_equivalent_TwoBar_for_integral_list(ijkl);
         if(ijkl == s) { // found an element, now see if it is non-zero
           ValueType v = H(i,j,k,l)-H(i,j,l,k);
           if( std::abs(v) > cutoff2bar )
             cnt++;  
         }
       }
       }
       }
       }
      }

    }

    // allocate shared memory 
    V2_2bar.reserve(cnt); 

    if(head_of_nodes) {
      for(IndexType i=0; i<NMO; i++)  {  
      for(IndexType j=i+1; j<NMO; j++)  {  
      for(IndexType k=i; k<NMO; k++)  {  
      for(IndexType l=k+1; l<NMO; l++)  {  
        if( i==j || k==l ) continue;
        s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0)); 
        s4D<ValueType> s = find_smaller_equivalent_TwoBar_for_integral_list(ijkl);
        if(ijkl == s) { // found an element, now see if it is non-zero
          ValueType v = H(i,j,k,l)-H(i,j,l,k);
          if( std::abs(v)>static_cast<RealType>(0.0) && std::abs(v) > cutoff2bar ) {
            std::get<4>(ijkl) = v; 
            V2_2bar.push_back(ijkl);
          }
        }
      }
      }
      }
      }
      if(!spinRestricted) {
        for(IndexType i=NMO; i<2*NMO; i++)  {  
        for(IndexType j=i+1; j<2*NMO; j++)  {  
        for(IndexType k=i; k<2*NMO; k++)  {  
        for(IndexType l=k+1; l<2*NMO; l++)  {
          if( i==j || k==l ) continue;
          s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0)); 
          s4D<ValueType> s = find_smaller_equivalent_TwoBar_for_integral_list(ijkl);
          if(ijkl == s) { // found an element, now see if it is non-zero
            ValueType v = H(i,j,k,l)-H(i,j,l,k);
            if( std::abs(v)>static_cast<RealType>(0.0) && std::abs(v) > cutoff2bar ) {
              std::get<4>(ijkl) = v; 
              V2_2bar.push_back(ijkl);
            }
          }
        }
        }
        }
        }
      }
    }
    myComm->barrier();

    Timer.stop("Generic");
    app_log()<<" -- Time to sort generate sparse 2-bar integral tables: " <<Timer.average("Generic") <<"\n";

#ifdef AFQMC_DEBUG
    app_log()<<" 2bar list contains " <<V2_2bar.size() <<" (" <<cnt <<") elements. \n" <<std::endl; 
    app_log()<<" Finished generating sparse 2-bar integrals. " <<std::endl; 
    app_log()<<" Sorting 2-bar integrals. " <<std::endl; 
#endif

    Timer.reset("Generic");
    Timer.start("Generic");

    if(head_of_nodes) std::sort (V2_2bar.begin(), V2_2bar.end(),mySort);
    myComm->barrier();

    Timer.stop("Generic");
    app_log()<<" -- Time to sort sparse 2-bar integral tables: " <<Timer.average("Generic") <<"\n";

#ifdef AFQMC_DEBUG
    app_log()<<" Finished sorting 2-bar integrals. " <<std::endl; 
#endif
    return true;
  }

  void SparseGeneralHamiltonian::get_selCI_excitations(OrbitalType I, OrbitalType J, int spinSector, RealType cutoff, OrbitalType* occs, std::vector<OrbitalType>& KLs ) {

    if(!has_hamiltonian_for_selCI) 
      APP_ABORT("Error: SparseGeneralHamiltonian::get_selCI_excitations() 2eInts not setup. \n\n\n");

    KLs.clear();
    KLs.reserve(2*nmax_KL_selCI);
    SMDenseVector<s2D<ValueType> >::iterator itV, itend;
    OrbitalType NK=0,NL=0;
    if( spinSector == 0 ) {  // aa
      IndexType pos = mapUT_woD(I,J,NMO); 
      IndexType n0=IJ_aa[pos];
      IndexType n1=IJ_aa[pos+1];
      if(n0==n1) return; 
      itV = V2_selCI_aa.begin() + n0;   
      itend = V2_selCI_aa.begin() + n1;   
    } else if(spinSector == 1) { // ab 
      IndexType pos;
      if(spinRestricted) {
        pos = mapUT(I,J-NMO,NMO);
        if(I <= J-NMO) NL=NMO; 
        else NK=NMO; 
      } else 
        pos = (I*NMO)*(J-NMO); 
      IndexType n0=IJ_ab[pos];
      IndexType n1=IJ_ab[pos+1];
      if(n0==n1) return;
      itV = V2_selCI_ab.begin() + n0;   
      itend = V2_selCI_ab.begin() + n1;   
    } else if(spinSector == 2) { // ba
      APP_ABORT(" Error: SparseGeneralHamiltonian::get_selCI_excitations().  Should not be here. \n\n\n");
    } else if(spinSector == 3) { // bb
      IndexType pos = mapUT_woD(I-NMO,J-NMO,NMO);
      if(spinRestricted) {
        IndexType n0=IJ_aa[pos];
        IndexType n1=IJ_aa[pos+1];
        if(n0==n1) return; 
        itV = V2_selCI_aa.begin() + n0;  
        itend = V2_selCI_aa.begin() + n1;
        NK=NL=NMO;
      } else {
        IndexType n0=IJ_bb[pos];
        IndexType n1=IJ_bb[pos+1];
        if(n0==n1) return; 
        itV = V2_selCI_bb.begin() + n0;  
        itend = V2_selCI_bb.begin() + n1;
      }
    }

    register OrbitalType K,L;
    while(itV < itend && std::abs(std::get<2>(*itV)) > cutoff) {
      std::tie(K,L,std::ignore) = *(itV++);
      K+=NK;  // needed for spinRestricted-ab
      L+=NL;  // needed for spinRestricted-ab
      if( !std::binary_search(occs,occs+NAEA+NAEB,K) && !std::binary_search(occs,occs+NAEA+NAEB,L) ) {
        KLs.push_back(K);
        KLs.push_back(L);
      }
    }
  }

  void SparseGeneralHamiltonian::generate_selCI_Ham(double cutoff) {

//  IndexMatrix IJ_aa, IJ_bb, IJ_ab;
//  SMDenseVector<s2D<ValueType> > V2_selCI_aa;
//  SMDenseVector<s2D<ValueType> > V2_selCI_ab;
//  SMDenseVector<s2D<ValueType> > V2_selCI_bb;

    nmax_KL_selCI=0;
    has_hamiltonian_for_selCI = true;
    IJ_aa.resize(mapUT_woD(NMO-2,NMO-1,NMO)+2);
    if(!spinRestricted) IJ_bb.resize( mapUT_woD(NMO-2,NMO-1,NMO)+2 );
    if(spinRestricted)
      IJ_ab.resize( mapUT(NMO-1,NMO-1,NMO)+2 );
    else
      IJ_ab.resize( NMO*NMO+1 );
    V2_selCI_aa.setup(head_of_nodes,std::string("V2_selCI_aa_"),TG.getNodeCommLocal());
    V2_selCI_ab.setup(head_of_nodes,std::string("V2_selCI_ab_"),TG.getNodeCommLocal());
    if(!spinRestricted)
      V2_selCI_bb.setup(head_of_nodes,std::string("V2_selCI_bb_"),TG.getNodeCommLocal());

    if(factorizedHamiltonian) {
      app_error()<<" Error in SparseGeneralHamiltonian::generate_selCI_Ham: not implemented \n"; 
      APP_ABORT(" Error in SparseGeneralHamiltonian::generate_selCI_Ham: not implemented \n"); 
    }

#ifdef AFQMC_DEBUG
    app_log()<<" Generating selCI integrals. " <<std::endl; 
#endif

    int cnt=0;
    Timer.reset("Generic");
    Timer.start("Generic");

    IndexType cntaa=0; 
    IndexType cntab=0; 
    IndexType cntbb=0; 

    // using dumbest algorithm right now, improve later
    if(head_of_nodes) {

// notes:
// 1. For aa/bb: only i<j and k<l are stored
// 2. For spinRestricted-ab only i<j is stored
      int tmp=0;
      for(OrbitalType i=0; i<NMO; i++)  {  
      for(OrbitalType j=i+1; j<NMO; j++)  {  
        if(mapUT_woD(i,j,NMO) != tmp++) 
          APP_ABORT("Error: Problems with mapUT_woD. \n\n\n"); 
      for(OrbitalType k=0; k<NMO; k++)  {  
      for(OrbitalType l=k+1; l<NMO; l++)  {  
        ValueType v = H(i,j,k,l)-H(i,j,l,k);
        if(std::abs(v) > cutoff) 
          cntaa++;  
      }
      }
      }
      }
      if(!spinRestricted) {
       for(OrbitalType i=0; i<NMO; i++)  {
       for(OrbitalType j=NMO; j<2*NMO; j++)  {
       for(OrbitalType k=0; k<NMO; k++)  {
       for(OrbitalType l=NMO; l<2*NMO; l++)  {  
         ValueType v = H(i,j,k,l);
         if(std::abs(v) > cutoff)
           cntab++;
       }
       }
       }
       }
       for(OrbitalType i=NMO; i<2*NMO; i++)  {  
       for(OrbitalType j=i+1; j<2*NMO; j++)  {  
       for(OrbitalType k=NMO; k<2*NMO; k++)  {  
       for(OrbitalType l=k+1; l<2*NMO; l++)  {  
        ValueType v = H(i,j,k,l)-H(i,j,l,k);
        if(std::abs(v) > cutoff)
          cntbb++;
       }
       }
       }
       }
      } else {
       tmp=0;
       for(OrbitalType i=0; i<NMO; i++)  {
       for(OrbitalType j=i; j<NMO; j++)  {
        if(mapUT(i,j,NMO) != tmp++) 
          APP_ABORT("Error: Problems with mapUT. \n\n\n"); 
       for(OrbitalType k=0; k<NMO; k++)  {
       for(OrbitalType l=0; l<NMO; l++)  {
         ValueType v = H(i,j,k,l);
         if(std::abs(v) > cutoff)
           cntab++;
       }
       }
       }
       }
      }

    }


    V2_selCI_aa.resize(cntaa);
    V2_selCI_ab.resize(cntab);
    if(!spinRestricted)
      V2_selCI_aa.resize(cntbb);

    app_log()<<" Memory used by selectedCI integrals: " 
     <<((cntaa+cntbb+cntab)*sizeof(s2D<ValueType>))/1024.0/1024.0 <<" MB. \n";

    if(head_of_nodes) {

      cnt = 0;
      IndexVector::iterator itIJ = IJ_aa.begin();
      SMDenseVector<s2D<ValueType> >::iterator itV = V2_selCI_aa.begin();
      for(OrbitalType i=0; i<NMO; i++)  {
      for(OrbitalType j=i+1; j<NMO; j++)  {
        *(itIJ++) = cnt;
        for(OrbitalType k=0; k<NMO; k++)  {
        for(OrbitalType l=k+1; l<NMO; l++)  {
          ValueType v = H(i,j,k,l)-H(i,j,l,k);
          if(std::abs(v) > cutoff) {
            *(itV++) = std::forward_as_tuple(k,l,v); 
            cnt++;
          }
        }
        }
      }
      }
      *(itIJ) = cnt;

      // sort integrals
      IndexVector::iterator itend = IJ_aa.end()-1;
      for(itIJ = IJ_aa.begin(), itV = V2_selCI_aa.begin(); itIJ < itend; itIJ++) { 
        int nt = *(itIJ+1) - *itIJ;
        if(nt > nmax_KL_selCI) nmax_KL_selCI = nt;
        std::sort(itV, itV+nt,  
               [] (const s2D<ValueType>& a, const s2D<ValueType>& b)
               {return std::abs(std::get<2>(a))>std::abs(std::get<2>(b));} );
        itV+=nt;
      }  

      if(!spinRestricted) {

        cnt=0;
        itIJ = IJ_ab.begin();
        itV = V2_selCI_ab.begin(); 
        for(OrbitalType i=0; i<NMO; i++)  {
        for(OrbitalType j=NMO; j<2*NMO; j++)  {
          *(itIJ++) = cnt;
          for(OrbitalType k=0; k<NMO; k++)  {
          for(OrbitalType l=NMO; l<2*NMO; l++)  {
            ValueType v = H(i,j,k,l);
            if(std::abs(v) > cutoff) {
              *(itV++) = std::forward_as_tuple(k,l,v);
              cnt++;
            }
          }
          }
        }
        }
        *(itIJ) = cnt;

        cnt=0;
        itIJ = IJ_bb.begin();
        itV = V2_selCI_bb.begin(); 
        for(OrbitalType i=NMO; i<2*NMO; i++)  {
        for(OrbitalType j=i+1; j<2*NMO; j++)  {
          *(itIJ++) = cnt;
          for(OrbitalType k=NMO; k<2*NMO; k++)  {
          for(OrbitalType l=k+1; l<2*NMO; l++)  {
            ValueType v = H(i,j,k,l)-H(i,j,l,k);
            if(std::abs(v) > cutoff) {
              *(itV++) = std::forward_as_tuple(k,l,v);
              cnt++;
            }
          }
          }
        }
        }
        *(itIJ) = cnt;

        // sort integrals
        itend = IJ_bb.end()-1;
        for(itIJ = IJ_bb.begin(), itV = V2_selCI_bb.begin(); itIJ < itend; itIJ++) { 
          int nt = *(itIJ+1) - *itIJ;
          if(nt > nmax_KL_selCI) nmax_KL_selCI = nt;
          std::sort(itV, itV+nt,  
                 [] (const s2D<ValueType>& a, const s2D<ValueType>& b)
                 {return std::abs(std::get<2>(a))>std::abs(std::get<2>(b));} );
          itV+=nt;
        }  

        itend = IJ_ab.end()-1;
        for(itIJ = IJ_ab.begin(), itV = V2_selCI_ab.begin(); itIJ < itend; itIJ++) {
          int nt = *(itIJ+1) - *itIJ;
          if(nt > nmax_KL_selCI) nmax_KL_selCI = nt;
          std::sort(itV, itV+nt,  
                 [] (const s2D<ValueType>& a, const s2D<ValueType>& b)
                 {return std::abs(std::get<2>(a))>std::abs(std::get<2>(b));} );
          itV+=nt;
        }  

      } else {

        // to be able to use i<j, I need to map I and J to [0:NMO-1]
        // this means that you need to be careful to map back to original sector
        
        cnt=0;
        itIJ = IJ_ab.begin();
        itV = V2_selCI_ab.begin(); 
        for(OrbitalType i=0; i<NMO; i++)  {
        for(OrbitalType j=i; j<NMO; j++)  {
          *(itIJ++) = cnt;
          for(OrbitalType k=0; k<NMO; k++)  {
          for(OrbitalType l=0; l<NMO; l++)  {
            ValueType v = H(i,j,k,l);
            if(std::abs(v) > cutoff) {
              *(itV++) = std::forward_as_tuple(k,l,v);
              cnt++;
            }
          }
          }
        }
        }
        *(itIJ) = cnt;

        itend = IJ_ab.end()-1;
        for(itIJ = IJ_ab.begin(), itV = V2_selCI_ab.begin(); itIJ < itend; itIJ++) {
          int nt = *(itIJ+1) - *itIJ;
          if(nt > nmax_KL_selCI) nmax_KL_selCI = nt;
          std::sort(itV, itV+nt,  
                 [] (const s2D<ValueType>& a, const s2D<ValueType>& b)
                 {return std::abs(std::get<2>(a))>std::abs(std::get<2>(b));} );
          itV+=nt;
        }  

      }
    }

    myComm->bcast<int>(nmax_KL_selCI);
    myComm->bcast<IndexType>(IJ_aa.data(),IJ_aa.size(),MPI_COMM_HEAD_OF_NODES); 
    myComm->bcast<IndexType>(IJ_ab.data(),IJ_ab.size(),MPI_COMM_HEAD_OF_NODES); 
    if(!spinRestricted)
      myComm->bcast<IndexType>(IJ_bb.data(),IJ_bb.size(),MPI_COMM_HEAD_OF_NODES); 
    myComm->barrier();

    Timer.stop("Generic");
    app_log()<<" -- Time to generate sparse 2-bar integral tables for selectedCI: " <<Timer.average("Generic") <<"\n";

  }

  // This routine operates on the FULL MO set, including CORE, ACTIVE and IGNORED states. 
  bool SparseGeneralHamiltonian::countElementsFromFCIDUMP(std::ifstream& in, int& n1, int& n1_c, int& n2, int& n2_c, int& n2_m, int& n3, int& n3_c, int& n3_m, std::map<IndexType,IndexType>& orbMapA, std::map<IndexType,IndexType>& orbMapB, int& n3_vec )
  {

     IndexType a,b,c,d;
     IndexType ap,bp,cp,dp;
     ValueType val;
     int uhf_block=0;
     n3_vec=n1=n2=n1_c=n2_c=n2_m=n3=n3_c=n3_m=0;

     while(!in.eof()) {
       in>>val >>ap >>bp >>cp >>dp;
       if(in.fail()) break;

       // to reduce problems from numerical roundoff 
       if(std::abs(toComplex(val).imag()) < 1e-8) setImag(val,0.0);

       // FCIDUMP format encodes UHF with uhf_block. All indexes are smaller than NMAX/NMO 
       // THIS SHOULD BE DONE AFTER MAPPING!!!
       // BUT YOU CAN'T CREATE THE MAPPING WITHOUT READING ALL STATES, CORRECT???
       // HOW TO SOLVE THIS??? 
       // READ ALL THE FIRST TIME????
      
       // trying to fix problem in 3-index terms
       auto it = orbMapA.find(cp);
       if (it == orbMapA.end())
        orbMapA[cp] = cp; 
       it = orbMapB.find(cp);
       if (it == orbMapB.end())
        orbMapB[cp] = cp; 

       if(spinRestricted) {
        a=orbMapA[ap];
        b=orbMapA[bp];
        c=orbMapA[cp];
        d=orbMapA[dp];
       } else {
         if(uhf_block == 0) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapA[cp];
           d=orbMapA[dp];
         } else if(uhf_block==1) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==2) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==3) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=3.\n";
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n";
              return false;
           }
           c=d=0;
         } else if(uhf_block==4) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=4.\n";
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n";
              return false;
           }
           c=d=0;
         }
       }

       if(a>0 && b>0 && c>0 && d > 0) {  // 2-electron ME (ab|cd) = <ac|bd>
         if( std::abs(val) > cutoff1bar ) {
           if( !goodSpinSector(a-1,c-1,b-1,d-1,NMO_FULL) ) {
             app_error()<<" Problems in SparseGeneralHamiltonian::countElementsFromFCIDUMP. Inconsistent two body term in integral file: " <<a <<" " <<b <<" " <<c <<" " <<d <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<" " <<val <<std::endl;
             return false;
           }
           int tmp1,tmp2,tmp3,tmp4;
           switch(uhf_block) {
             case 0: 
             {  // AAAA
               tmp1 = tmp2= NCA;
               break;
             }
             case 1: 
             { // BBBB
               tmp1 = tmp2 = NMO_FULL+NCB;
               tmp3=0;
               break;
             }
             case 2: 
             { // AABB
               tmp1 = NCA;
               tmp2 = NMO_FULL+NCB;
               break;
             }
           };
           if( a<=tmp1 && b<=tmp1 && c<=tmp2 && d<=tmp2 ) n2_c++;
           else if( a>tmp1 && b>tmp1 && c>tmp2 && d>tmp2 ) n2++; 
           else n2_m++;
         }
       } else if(a>0 && b>0 && c>0) { // factorized 2-electron ME (i,k|n), where (i,k|j,l) = sum_n (i,k|n)*(l,j|n)*  
         if( std::abs(val) > cutoff1bar ) {
           if( !goodSpinSector(a-1,a-1,b-1,b-1,NMO_FULL) ) {
             app_error()<<" Problems in SparseGeneralHamiltonian::countElementsFromFCIDUMP. Inconsistent factorized two body term in integral file: " <<a <<" " <<b <<" " <<c <<" " <<d <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<" " <<val <<std::endl;
             return false;
           }
           int tmp1,tmp2,tmp3,tmp4;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2= NCA;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO_FULL+NCB;
               tmp3=0;
               break;
             }
           };
           if( a<=tmp1 && b<=tmp1 ) n3_c++;
           else if( a>tmp1 && b>tmp1 ) n3++;
           else n3_m++;
           if(cp > n3_vec) n3_vec = cp;
         }
       } else if(a>0 && b>0) { // 1-electron ME (a|b) = <a|b>
         if( std::abs(val) > cutoff1bar ) {

           if(spinRestricted) {
             if( a <= NCA && b <= NCA ) n1_c++;
             else if( a>NCA && b>NCA  ) n1++; 
           } else {
             if(uhf_block == 3) {
               if( a <= NCA && b <= NCA ) n1_c++;
               else if( a>NCA && b>NCA  ) n1++; 
             } else {
               if( a <= NMO_FULL+NCB && b <= NMO_FULL+NCB ) n1_c++;
               else if( a>NMO_FULL+NCB && b>NMO_FULL+NCB  ) n1++; 
             }
           }
 
         }
       } else if(a==0 && b==0 && c==0 && d==0) {
         if( val==0.0 && !spinRestricted ) {
           uhf_block++;
         }
       } else if(a!=0 && b==0 && c==0 && d==0) {
       } else {
         app_error()<<"Error with ASCII integral file, bad indexing. " <<std::endl;
         return false;
       }
    }
    return true;
  }

 // This routine assumes that states are ordered (after mapping) in the following way:
 // { core, active+virtual, ignored}, so make sure mappings obey this.   
 // An easy way to implement this is to have a list of all core states and 
 // initialize the mapping (from 1:NC) with the index of the core states.
 // Then complete the mapping with every other orbital index not in the core list.
 // This way you always end up with a mapping with the correct format. 
  bool SparseGeneralHamiltonian::readElementsFromFCIDUMP(std::ifstream& in, 
         std::vector<s2D<ValueType> >& V1,
         std::vector<s2D<ValueType> >& V1_c,
         SMDenseVector<s4D<ValueType> >& V2,
         std::vector<s4D<ValueType> >& V2_c,
         std::vector<s4D<ValueType> >& V2_m,
         ValueSMSpMat&  V3,
         ValueSMSpMat&  V3_c,
         ValueSMSpMat&  V3_m,
         std::map<IndexType,IndexType>& orbMapA, std::map<IndexType,IndexType>& orbMapB) {   

     IndexType a,b,c,d, cntS=0, cntD=0,q1;
     IndexType ap,bp,cp,dp,ab,cd;

     std::vector<s2D<ValueType> >::iterator V1_it = V1.begin();
     SMDenseVector<s4D<ValueType> >::iterator V2_it;
     if(V2.isAllocated()) V2_it = V2.begin();

     std::vector<s2D<ValueType> >::iterator V1_itc = V1_c.begin();
     std::vector<s4D<ValueType> >::iterator V2_itc = V2_c.begin();
     std::vector<s4D<ValueType> >::iterator V2_itm = V2_m.begin();

     ValueType val;
     int uhf_block=0;
     while(!in.eof()) {
       in>>val >>ap >>bp >>cp >>dp;
       if(in.fail()) break;

       // to reduce problems from numerical roundoff 
       if(std::abs(toComplex(val).imag()) < 1e-8) setImag(val,0.0);

       if(spinRestricted) { 
        a=orbMapA[ap];
        b=orbMapA[bp];
        c=orbMapA[cp];
        d=orbMapA[dp];
       } else {
         if(uhf_block == 0) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapA[cp];
           d=orbMapA[dp];
         } else if(uhf_block==1) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==2) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           c=orbMapB[cp];
           d=orbMapB[dp];
         } else if(uhf_block==3) {
           a=orbMapA[ap];
           b=orbMapA[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=3.\n"; 
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n"; 
              return false;
           }
           c=d=0;
         } else if(uhf_block==4) {
           a=orbMapB[ap];
           b=orbMapB[bp];
           if(cp != 0 && dp != 0) {
              app_error()<<"Error: In UHF FCIDUMP, c d should be zero in uhf_block=4.\n"; 
              app_error()<<"val, a, b, c, d: " <<val <<" " <<ap <<" " <<bp <<" " <<cp <<" " <<dp <<"\n"; 
              return false;
           }
           c=d=0;
         }
       }

       if(a>0 && b>0 && c > 0 && d>0) {  // 2-electron ME (ab|cd) = <ac|bd>
#if true
         if( std::abs(val) > cutoff1bar ) {

           int tmp1,tmp2,tmp4=0,tmp3=0;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2 = NCA;
               tmp3 = tmp4 = NCA;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO_FULL+NCB;
               tmp3 = tmp4 = NCA+NCB; 
               break;
             }
             case 2:
             { // AABB
               tmp1 = NCA;
               tmp3 = NCA;
               tmp2 = NMO_FULL+NCB;
               tmp4 = NCB+NCA;
               break;
             }
             default:
             {
               app_error()<<"Problems reading FCIDUMP: " 
               <<ap <<" " 
               <<bp <<" " 
               <<cp <<" " 
               <<dp <<" " 
               <<val <<std::endl;
               return false;
             }
           };
           // later on
           // if( isCore[a] && isCore[b] && isCore[c] && isCore[d] ) 
           if( a<=tmp1 && b<=tmp1 && c<=tmp2 && d<=tmp2 ) 
             *(V2_itc++) = find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-1,c-1,b-1,d-1,val) ); 
           else if( a>tmp1 && b>tmp1 && c>tmp2 && d>tmp2 ) { 
             if(head_of_nodes) *(V2_it++) = find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-tmp3-1,c-tmp4-1,b-tmp3-1,d-tmp4-1,val) ); 
           } else 
             *(V2_itm++) = find_smaller_equivalent_OneBar_for_integral_list( std::make_tuple(a-1,c-1,b-1,d-1,val) ); 
         }
#else
// reserved for non c++11 compliant compiler
#endif
       } else if(a > 0 && b > 0 && c>0) { // factorized 2-electron integral 
         if( std::abs(val) > cutoff1bar ) {

           int tmp1,tmp2,tmp4=0,tmp3=0;
           switch(uhf_block) {
             case 0:
             {  // AAAA
               tmp1 = tmp2 = NCA;
               tmp3 = tmp4 = NCA;
               break;
             }
             case 1:
             { // BBBB
               tmp1 = tmp2 = NMO_FULL+NCB;
               tmp3 = tmp4 = NCA+NCB;
               break;
             }
             default:
             {
               app_error()<<"Problems reading FCIDUMP: "
               <<ap <<" "
               <<bp <<" "
               <<cp <<" "
               <<dp <<" "
               <<val <<std::endl;
               return false;
             }
           };

           if(head_of_nodes) { 
            if( a<=tmp1 && b<=tmp1 )
              V3_c.add((a-1)*NMO+Index2Col(b-1),cp-1,val);
            else if( a>tmp1 && b>tmp1 )
              V3.add((a-tmp3-1)*NMO+Index2Col(b-tmp3-1),cp-1,val);
            else
              V3_m.add((a-1)*NMO+Index2Col(b-1),cp-1,val);
           }

         }
       } else if(a > 0 && b > 0) { // 1-electron ME (a|b) = <a|b>
#if true
         if( std::abs(val) > cutoff1bar ) { 

           if(a > b) {q1=b;b=a;a=q1;val = myconj(val);}

           if(spinRestricted) {
             if( a <= NCA && b <= NCA ) {
               std::get<0>( *V1_itc ) = a-1;
               std::get<1>( *V1_itc ) = b-1;
               std::get<2>( *(V1_itc++) ) = val;
             } else if( a>NCA && b>NCA  ) {
               std::get<0>( *V1_it ) = a-NCA-1;
               std::get<1>( *V1_it ) = b-NCA-1;
               std::get<2>( *(V1_it++) ) = val;
             } 
           } else {
             if(uhf_block == 3) {
               if( a <= NCA && b <= NCA ) { 
                 std::get<0>( *V1_itc ) = a-1;
                 std::get<1>( *V1_itc ) = b-1;
                 std::get<2>( *(V1_itc++) ) = val;
               } else if( a>NCA && b>NCA  ) { 
                 std::get<0>( *V1_it ) = a-NCA-1;
                 std::get<1>( *V1_it ) = b-NCA-1;
                 std::get<2>( *(V1_it++) ) = val;
               }
             } else {
               if( a <= NMO+NCB && b <= NMO+NCB ) { 
                 std::get<0>( *V1_itc ) = a-1;
                 std::get<1>( *V1_itc ) = b-1;
                 std::get<2>( *(V1_itc++) ) = val;
               } else if( a>NMO+NCB && b>NMO+NCB  ) { 
                 std::get<0>( *V1_it ) = a-NCB-NCA-1;
                 std::get<1>( *V1_it ) = b-NCB-NCA-1;
                 std::get<2>( *(V1_it++) ) = val;
               }
             }
           }
         }
#else
// reserved for non c++11 compliant compiler
#endif
       } else if(a==0 && b==0 && c==0 && d==0) {
         if( val==0.0 && !spinRestricted ) {
           uhf_block++;
         } else {
           NuclearCoulombEnergy = val;
         } 
       } else if(a!=0 && b==0 && c==0 && d==0) {
         // ignore, these are the eigenvalues of the Fock operator printed by VASP  
       } else {
         app_error()<<"Error with ASCII integral file, bad indexing. " <<std::endl;
         return false;
       }
     } 

/*
     return (V2_it == V2.end()) &&
            (V1_it == V1.end()) &&
            (V2_itc == V2_c.end()) &&
            (V2_itm == V2_m.end()) &&
            (V1_itc == V1_c.end());
*/
     return true;
  }

  bool SparseGeneralHamiltonian::initFromHDF5(const std::string& fileName)
  {

   Timer.reset("Generic");
   Timer.start("Generic");

    // no hamiltonian distribution yet
    min_i = 0;
    max_i = NMO;
    if(!spinRestricted) max_i *= 2;  // or 4?

    app_log()<<" Initializing Hamiltonian from file: " <<fileName <<std::endl;

    V2.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2"),TG.getNodeCommLocal());
    V2_2bar.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2_2bar"),TG.getNodeCommLocal());
    V2_full.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2_full"),TG.getNodeCommLocal());
    V2_fact.setup(head_of_nodes,"SparseGeneralHamiltonian_V2_fact",TG.getNodeCommLocal());

    if(rotation != "" && head_of_nodes ) {

      rotationMatrix.resize(NMO,NMO);
      std::ifstream in(rotation.c_str());
      if(in.fail()) {
        app_error()<<" Error opening rotation file. \n";
        return false;
      }

      for(int i=0; i<NMO; i++) 
       for(int j=0; j<NMO; j++) 
         in>>rotationMatrix(i,j);

      in.close();
    } 

    if(myComm->rank() == 0) {

      hdf_archive dump(myComm);
      if(!dump.open(fileName,H5F_ACC_RDONLY,false)) {
        app_error()<<" Error opening integral file in SparseGeneralHamiltonian. \n";
        return false;
      }

      std::string path = "/Hamiltonian/SparseGeneralHamiltonian";
      if(!dump.is_group( path )) {
        app_error()<<" ERROR: H5Group /Hamiltonian/SparseGeneralHamiltonian does not exists in restart file. \n";
        return false;
      }

      if(!dump.push("Hamiltonian",false)) {
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Group not Hamiltonian found. \n"; 
        return false;
      }
      if(!dump.push("SparseGeneralHamiltonian",false)) { 
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Group not SparseGeneralHamiltonian found. \n"; 
        return false;
      }

      std::vector<int> Idata(8);
      if(!dump.read(Idata,"dims")) return false;

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
      myComm->bcast(Idata);
      spinRestricted = (Idata[6]==0)?(true):(false);
      factorizedHamiltonian = (Idata[7]>0); 
      int nvecs = Idata[7];

      H1.resize(Idata[0]);

      occup_alpha.resize(NAEA);
      occup_beta.resize(NAEB);
      Idata.resize(NAEA+NAEB);
      if(!dump.read(Idata,"occups")) { 
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading occups dataset. \n"; 
        return false;
      }
      for(int i=0; i<NAEA; i++) occup_alpha[i] = Idata[i];
      for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) occup_beta[j] = Idata[i];

      std::vector<ValueType> Rdata(2);
      if(!dump.read(Rdata,"Energies")) { 
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading  dataset. \n"; 
        return false;
      }
      NuclearCoulombEnergy = Rdata[0];
      FrozenCoreEnergy = Rdata[0];

      myComm->bcast(occup_alpha);
      myComm->bcast(occup_beta);
      myComm->bcast(NuclearCoulombEnergy);
      myComm->bcast(FrozenCoreEnergy);

      std::vector<OrbitalType> ivec;    
      ivec.resize(2*H1.size());
      if(!dump.read(ivec,"H1_indx")) return false;
      for(int i=0, j=0; i<H1.size(); i++, j+=2)        
        H1[i] = std::forward_as_tuple(ivec[j],ivec[j+1],0);  

      std::vector<ValueType> vvec;
      vvec.resize(H1.size());
      if(!dump.read(vvec,"H1")) return false;
      for(int i=0; i<H1.size(); i++)
        std::get<2>(H1[i]) = vvec[i];

      std::sort (H1.begin(), H1.end(),mySort);
      myComm->bcast<char>(reinterpret_cast<char*>(H1.data()),H1.size()*sizeof(s2D<ValueType>));

      // setup and communicate min_i/max_i here if nnodes_per_TG > 1

      if(factorizedHamiltonian) {

        std::vector<double> residual(nvecs);
        if(!dump.read(residual,"V2fact_vec_residual")) {
          app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2fact_vec_residual dataset. \n";      
          return false;
        }

        double cut = cutoff_cholesky;
        int nvecs_after_cutoff=std::count_if(residual.begin(), residual.end(), 
           [cut] (double i) { return i>cut; } );
        myComm->bcast(nvecs_after_cutoff);

        std::vector<int> sz(nvecs);
        if(!dump.read(sz,"V2fact_vec_sizes")) {
          app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2fact_vec_sizes dataset. \n";      
          return false;
        }
        int nmax = *std::max_element(sz.begin(),sz.end());
        ivec.reserve(nmax);
        vvec.reserve(nmax);

        // eliminate std::vectors with small residuals
        std::vector<int> sz2;
        sz2.reserve(nvecs_after_cutoff);
        cholesky_residuals.reserve(nvecs_after_cutoff);
        for(int i=0; i<nvecs; i++) 
          if(residual[i] > cutoff_cholesky) {
            cholesky_residuals.push_back(residual[i]); 
            sz2.push_back(sz[i]);
          }      
        myComm->bcast(sz2);

        int NMO2 = NMO*NMO;
        if(!spinRestricted) NMO2 *= 2;
        V2_fact.setDims(NMO2,nvecs_after_cutoff);

        // this is an upper bound
        int mysz = std::accumulate(sz2.begin(),sz2.end(),0); 
        V2_fact.reserve(mysz);
        ValueMatrix Temp,Temp1;
        if(rotation!="") {
          if(!spinRestricted) 
            APP_ABORT("Error: rotation only implemented with spinRestricted. \n\n\n"); 
          Temp.resize(NMO,NMO); 
          Temp1.resize(NMO,NMO); 
          ivec.reserve(NMO2);
          vvec.reserve(NMO2);
        }
         

        Timer.reset("Generic3");
        Timer.reset("Generic4");

        for(int i=0; i<nvecs; i++) {
    
          if( residual[i] <= cutoff_cholesky ) continue;

          ivec.resize(sz[i]);
          vvec.resize(sz[i]);
          
          Timer.start("Generic3");
          if(!dump.read(ivec,std::string("V2fact_index_")+std::to_string(i))) {
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2fact_index_" <<std::to_string(i) <<" dataset. \n";      
            app_error()<<" Expected size: " <<sz[i] <<"\n";
            return false;
          } 
          if(!dump.read(vvec,std::string("V2fact_vals_")+std::to_string(i))) { 
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2fact_vals_ dataset. \n";      
            return false;
          }
          Timer.stop("Generic3");

          if(rotation!="") {
            Temp=0;
            for(int k=0; k<sz[i]; k++) 
              Temp(ivec[k]/NMO,ivec[k]%NMO) = vvec[k];
            DenseMatrixOperators::product(NMO,NMO,NMO,Temp.data(),NMO,rotationMatrix.data(),NMO,Temp1.data(),NMO); 
            DenseMatrixOperators::product_AhB(NMO,NMO,NMO,ValueType(1.0),rotationMatrix.data(),NMO,Temp1.data(),NMO,ValueType(0.0),Temp.data(),NMO); 
            ivec.clear();
            vvec.clear(); 
            for(int j=0; j<NMO; j++) 
             for(int k=0; k<NMO; k++) {
               if(std::abs(Temp(j,k)) > cutoff2bar) {
                 ivec.push_back(j*NMO+k);
                 vvec.push_back(Temp(j,k));
               } 
             } 
            app_log()<<" Change in Chol Vect. number of terms, old:" <<sz[i] <<"  new:" <<ivec.size() <<std::endl; 
 
          }

          Timer.start("Generic4");
          myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
          myComm->bcast(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);
          Timer.stop("Generic4");
 
//          if(nnodes_per_TG == 1 || (i>=min_i && i < max_i)) {
            for(int k=0; k<ivec.size(); k++) {
              if(std::abs(vvec[k]) > cutoff2bar ) 
                V2_fact.add(ivec[k],i,vvec[k],false);
            }
//          }

        }

        Timer.reset("Generic2");
        Timer.start("Generic2");
        V2_fact.compress_parallel(TG.getNodeCommLocal());
        Timer.stop("Generic2");
        app_log()<<" -- Average time to read std::vector from h5 file: " <<Timer.average("Generic3") <<"\n";
        app_log()<<" -- Average time to bcast std::vector: " <<Timer.average("Generic4") <<"\n";
        app_log()<<" -- Time to compress factorized Hamiltonian from h5 file: " <<Timer.average("Generic2") <<"\n";
        app_log()<<" Memory used by factorized 2-el integral table: " <<V2_fact.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;
      } else {
        // now V2
        Idata.resize(NMO); 
        if(!dump.read(Idata,"V2_location_of_indexes")) {
          app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_location_of_indexes dataset. \n";
          return false;
        } 
        myComm->bcast(Idata);
        int ntmax=Idata[0], nttot=0;
        if(min_i == 0) nttot += Idata[0];
        for(int i=1; i<NMO; i++) {
          if( i >= min_i && i < max_i ) nttot += Idata[i]-Idata[i-1];    
          if( Idata[i]-Idata[i-1] > ntmax) ntmax = Idata[i]-Idata[i-1];
        }
        V2.resize(nttot);
        SMDenseVector<s4D<ValueType> >::iterator V2_it = V2.begin();

        ivec.reserve(3*ntmax);
        vvec.reserve(ntmax);

        for(int k=0; k<NMO; k++) {
          int ik0 = (k==0)?0:Idata[k-1];
          if(Idata[k]==ik0) continue;
          ivec.clear();
          ivec.resize(3*(Idata[k]-ik0));
          if(!dump.read(ivec,std::string("V2_indx_")+std::to_string(k))) {
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_indx_ dataset. \n";
            return false;
          }
          myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);

          vvec.clear();
          vvec.resize(Idata[k]-ik0);
          if(!dump.read(vvec,std::string("V2_")+std::to_string(k))) {
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_ dataset. \n";
            return false;
          }
          myComm->bcast(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);

          if( k >= min_i && k < max_i) {
            std::vector<OrbitalType>::iterator iti = ivec.begin(); 
            for(std::vector<ValueType>::iterator itv = vvec.begin(); itv < vvec.end(); itv++, iti+=3)
              *(V2_it++) = std::make_tuple(k,*(iti),*(iti+1),*(iti+2),*itv); 
          }
        }
        std::sort (V2.begin(), V2.end(),mySort);
        app_log()<<" Memory used by 2-el integral table: " <<V2.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;
      }

      dump.pop();
      dump.pop();

      dump.close();


    } else {

      std::vector<int> Idata(8);
      myComm->bcast(Idata);
      spinRestricted = (Idata[6]==0)?(true):(false);
      if(NMO < 0) NMO = Idata[3];
      if(NAEA < 0) NAEA = Idata[4];
      if(NAEB < 0) NAEB = Idata[5];
      if(Idata[3] != NMO) {
        app_error()<<" ERROR: NMO differs from value in integral file. \n";
        APP_ABORT("ERROR: NMO differs from value in integral file. \n"); 
        return false;
      }
      if(Idata[4] != NAEA) {
        app_error()<<" ERROR: NMO differs from value in integral file. \n";
        APP_ABORT(" ERROR: NMO differs from value in integral file. \n");
        return false;
      }
      if(Idata[5] != NAEB) {
        app_error()<<" ERROR: NMO differs from value in integral file. \n";
        APP_ABORT(" ERROR: NMO differs from value in integral file. \n");
        return false;
      }
      factorizedHamiltonian = (Idata[7]>0); 
      int nvecs = Idata[7];

      occup_alpha.resize(NAEA);
      occup_beta.resize(NAEB);

      H1.resize(Idata[0]);

      myComm->bcast(occup_alpha);
      myComm->bcast(occup_beta);
      myComm->bcast(NuclearCoulombEnergy);
      myComm->bcast(FrozenCoreEnergy);

      myComm->bcast<char>(reinterpret_cast<char*>(H1.data()),H1.size()*sizeof(s2D<ValueType>));

      // recieve min_i/max_i here if nnodes_per_TG > 1

      if(factorizedHamiltonian) {

        int nvecs_after_cutoff;
        myComm->bcast(nvecs_after_cutoff);

        std::vector<int> sz(nvecs_after_cutoff);
        myComm->bcast(sz);

        int NMO2 = NMO*NMO;
        if(!spinRestricted) NMO2 *= 2;
        V2_fact.setDims(NMO2,nvecs_after_cutoff);

        // right now this is an upper bound, since terms below the 2bar cutoff can be ignored 
        int mysz = std::accumulate(sz.begin(),sz.end(),0); 
        V2_fact.reserve(mysz);

        if(head_of_nodes) {
          int nmax = *std::max_element(sz.begin(),sz.end());
          std::vector<IndexType> ivec;
          std::vector<ValueType> vvec;
          ivec.reserve(nmax);
          vvec.reserve(nmax);
          for(int i=0; i<nvecs_after_cutoff; i++) {
    
            ivec.resize(sz[i]);
            vvec.resize(sz[i]);
          
            myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
            myComm->bcast(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);
 
//            if(nnodes_per_TG == 1 || (i>=min_i && i < max_i)) {
              for(int k=0; k<ivec.size(); k++) {
                if(std::abs(vvec[k]) > cutoff2bar ) 
                  V2_fact.add(ivec[k],i,vvec[k],false); 
              }
//            }

          }

        }
        V2_fact.compress_parallel(TG.getNodeCommLocal());

      } else {

        // now V2
        Idata.resize(NMO); 
        myComm->bcast(Idata);
        int ntmax=Idata[0], nttot=0;
        if(min_i == 0) nttot += Idata[0];
        for(int i=1; i<NMO; i++) {
          if( i >= min_i && i < max_i ) nttot += Idata[i]-Idata[i-1];    
          if( Idata[i]-Idata[i-1] > ntmax) ntmax = Idata[i]-Idata[i-1];
        }
        V2.resize(nttot);

        if(head_of_nodes) { 
          SMDenseVector<s4D<ValueType> >::iterator V2_it = V2.begin();
          std::vector<OrbitalType> ivec;    
          std::vector<ValueType> vvec; 
          ivec.reserve(3*ntmax);
          vvec.reserve(ntmax);

          for(int k=0; k<NMO; k++) {
            int ik0 = (k==0)?0:Idata[k-1];
            if(Idata[k]==ik0) continue;
            ivec.clear();
            ivec.resize(3*(Idata[k]-ik0));
            myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
        
            vvec.clear();
            vvec.resize(Idata[k]-ik0);
            myComm->bcast(vvec.data(),vvec.size(),MPI_COMM_HEAD_OF_NODES);
        
            if( k >= min_i && k < max_i) {
              std::vector<OrbitalType>::iterator iti = ivec.begin();
              for(std::vector<ValueType>::iterator itv = vvec.begin(); itv < vvec.end(); itv++, iti+=3)
                *(V2_it++) = std::make_tuple(k,*(iti),*(iti+1),*(iti+2),*itv);           
            } 
          }
          std::sort (V2.begin(), V2.end(),mySort);
        }
      }
    }
    myComm->barrier();

    Timer.stop("Generic");
    app_log()<<" -- Time to initialize Hamiltonian from h5 file: " <<Timer.average("Generic") <<"\n";

    if(rank() == 0) {
      ascii_write();
    }

    return true;

  }

  void SparseGeneralHamiltonian::ascii_write() {

    if(ascii_write_file == std::string("")) return;
 
    std::ofstream out(ascii_write_file.c_str());
    if(out.fail()) {
      app_error()<<" Error opening restart file in SparseGeneralHamiltonian::ascii_write() . " <<std::endl; 
      return;
    }

// hack to write sparse matrix in ascii for testing
/*
   if(!factorizedHamiltonian) 
     APP_ABORT("Error: Matrix dump only for factorized hamiltonian \n\n\n");
   if(rank()!=0) return;
   int nvecs=V2_fact.rows();
   out<<"# nrows: " <<nvecs <<" ncols: " <<V2_fact.cols() <<"\n";
   for(ComplexSMSpMat::indxType i=0; i<nvecs; i++) {
     int k1 = *(V2_fact.rowIndex_begin()+i+1); 
     for(int k=*(V2_fact.rowIndex_begin()+i); k<k1; k++) 
       out<<i <<" " <<*(V2_fact.cols_begin()+k) <<" "
          <<std::setprecision(12) <<*(V2_fact.vals_begin()+k) <<"\n";
   }
   out<<std::endl;
   out.close();
   APP_ABORT(" FINISHED DUMPING ASCII FILE. \n\n\n");   
*/
 
    out<<"&FCI NORB=" <<NMO <<",NELEC=" <<NAEA+NAEB <<",MS2=" <<MS2 <<"," <<"\n"
       <<"ISYM=" <<ISYM <<"," <<"\n/\n"; 

    out.setf(std::ios::scientific, std::ios::floatfield);
    if(factorizedHamiltonian) {
      out<<" factorizedHamiltonian output not yet implemented " <<std::endl;
    } else {

      SMDenseVector<s4D<ValueType> >::iterator it=V2.begin();
      for( ; it!=V2.end(); it++) {
        if(std::get<0>(*it) >= NMO) break;
        out<<std::setprecision(16) <<std::get<4>(*it) <<std::setw(5) <<std::get<0>(*it)+1 <<std::setw(5) <<std::get<2>(*it)+1 
            <<std::setw(5) <<std::get<1>(*it)+1 <<std::setw(5) <<std::get<3>(*it)+1 <<"\n";
      }
      if(!spinRestricted) {
        out<<" UHF output not yet implemented " <<std::endl;
      }

    }

    s2Dit it=H1.begin();
    for( ; it!=H1.end(); it++) {
      if(std::get<0>(*it) >= NMO) break; 
      out<<std::setprecision(16) <<std::get<2>(*it) <<std::setw(5) <<std::get<0>(*it)+1 <<std::setw(5) <<std::get<1>(*it)+1 <<"    0    0\n"; 
    } 
    if(!spinRestricted) {
      out<<" 0.0                    0    0    0    0\n";
      for( ; it!=H1.end(); it++) {
        if(std::get<0>(*it) >= NMO) break;
        out<<std::setprecision(16) <<std::get<2>(*it) <<" " <<std::get<0>(*it)-NMO+1 <<" " <<std::get<1>(*it)-NMO+1 <<"   0   0\n";
      }
    }
    out<<std::setprecision(16) <<NuclearCoulombEnergy <<"    0    0    0    0\n";

    out.close();

  } 

  // this needs to be modified to wread/write in blocks
  // to avoid having to allocate the full array twice
  // Also, you could add support in hdf_archive for tuples 
  void SparseGeneralHamiltonian::hdf_write() {

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

    std::vector<int> Idata(8);
    Idata[0]=H1.size();
    if(factorizedHamiltonian)
      Idata[1]=V2_fact.size();
    else
      Idata[1]=V2.size();
    Idata[2]=0;
    Idata[3]=NMO;
    Idata[4]=NAEA;
    Idata[5]=NAEB;
    Idata[6]=spinRestricted?(0):(1);
    Idata[7]=0;
    if(factorizedHamiltonian)
      Idata[7] = V2_fact.cols(); 
    dump.write(Idata,"dims");

    Idata.resize(NAEA+NAEB);
    for(int i=0; i<NAEA; i++) Idata[i] = occup_alpha[i];
    for(int i=NAEA, j=0; i<NAEA+NAEB; i++, j++) Idata[i] = occup_beta[j];
    dump.write(Idata,"occups");

    std::vector<ValueType> Rdata(2);
    Rdata[0] = NuclearCoulombEnergy;
    Rdata[1] = FrozenCoreEnergy; 
    dump.write(Rdata,"Energies");

    // write H1 
    std::vector<OrbitalType> ivec; 
    ivec.resize(2*H1.size());
    for(int i=0, j=0; i<H1.size(); i++, j+=2) 
      std::tie (ivec[j],ivec[j+1],std::ignore) = H1[i];   
    dump.write(ivec,"H1_indx");

    std::vector<ValueType> vvec;
    vvec.resize(H1.size());
    for(int i=0; i<H1.size(); i++)
      std::tie (std::ignore,std::ignore,vvec[i]) = H1[i];
    dump.write(vvec,"H1");

    if(factorizedHamiltonian) {

      if(cholesky_residuals.size() != V2_fact.cols()) {
        app_error()<<"Error: Incorrect size of cholesky_residuals: " <<cholesky_residuals.size() <<" " <<V2_fact.cols() <<std::endl;  
        APP_ABORT("Error. \n\n\n");
      }
      dump.write(cholesky_residuals,std::string("V2fact_vec_residual"));

      V2_fact.transpose();

      int nvecs=V2_fact.rows();
      std::vector<int> sz(nvecs);
      for(int i=0; i<nvecs; i++) 
        sz[i]=*(V2_fact.rowIndex_begin()+i+1) - *(V2_fact.rowIndex_begin()+i);
      int nmax = *std::max_element(sz.begin(),sz.end());
      dump.write(sz,std::string("V2fact_vec_sizes"));
      std::vector<IndexType> ivec;
      std::vector<ComplexType> vvec;
      ivec.reserve(nmax);
      vvec.reserve(nmax);
      for(ComplexSMSpMat::indxType i=0; i<nvecs; i++) {

        ivec.resize(sz[i]);
        vvec.resize(sz[i]);
        std::copy(V2_fact.cols_begin()+(*(V2_fact.rowIndex_begin()+i)),V2_fact.cols_begin()+(*(V2_fact.rowIndex_begin()+i+1)),ivec.begin());
        std::copy(V2_fact.vals_begin()+(*(V2_fact.rowIndex_begin()+i)),V2_fact.vals_begin()+(*(V2_fact.rowIndex_begin()+i+1)),vvec.begin());

        dump.write(ivec,std::string("V2fact_index_")+std::to_string(i));
        dump.write(vvec,std::string("V2fact_vals_")+std::to_string(i));

      }

      V2_fact.transpose();

    } else { 
      // write V2
      // terms with i index are found between [Idata[i-1],Idata[i])  
      Idata.resize(NMO);
      int ntmax=0;
      OrbitalType i0,j0,k0,l0,iold=0;
      for(int i=0; i<V2.size(); i++) {
        std::tie (i0,j0,k0,l0,std::ignore) = V2[i];
        if(iold != i0) {
          for(int k=iold; k<i0; k++) Idata[k]=i; 
          iold=i0;
        }
      } 
      for(int k=iold; k<NMO; k++)  Idata[k]=V2.size();
      ntmax = Idata[0];
      for(int i=1; i<NMO; i++) {
        if(ntmax < (Idata[i]-Idata[i-1]))
          ntmax = Idata[i]-Idata[i-1];
      }

      dump.write(Idata,"V2_location_of_indexes");

      ivec.reserve(3*ntmax);
      vvec.reserve(ntmax);

      for(int k=0; k<NMO; k++) { 
        int ik0 = (k==0)?0:Idata[k-1]; 
        if(Idata[k]==ik0) continue; 
        ivec.clear();
        // no point in storing i
        ivec.resize(3*(Idata[k]-ik0));
        for (int i=ik0,jk=0; i<Idata[k]; i++,jk+=3) 
          std::tie (i0,ivec[jk],ivec[jk+1],ivec[jk+2],std::ignore) = V2[i];   
        dump.write(ivec,std::string("V2_indx_")+std::to_string(k));

        vvec.clear();
        vvec.resize(Idata[k]-ik0);
        for(int i=ik0; i<Idata[k]; i++) 
          std::tie (std::ignore,std::ignore,std::ignore,std::ignore,vvec[i-ik0]) = V2[i];
        dump.write(vvec,std::string("V2_")+std::to_string(k));
      }
    }

    dump.pop();
    dump.pop();

    dump.flush();
    dump.close();

  }

  void SparseGeneralHamiltonian::calculateOneBodyPropagator(RealType cut, RealType dt, ComplexMatrix& Hadd , std::vector<s2D<ComplexType> >& Pkin) 
  {

    int NMO2 = spinRestricted?(NMO):(2*NMO); 

    ComplexMatrix v(NMO2,NMO),P(NMO2,NMO);
    
    // v is the dense representation of H1+H1add+v0
    // 1. Hadd should be the contribution from mean-field substraction, otherwise zero  
    for(int i=0; i<NMO; i++) { 
     v(i,i) = Hadd(i,i); 
     for(int j=i+1; j<NMO; j++) { 
       if( std::abs( Hadd(i,j) - myconj(Hadd(j,i)) ) > 1e-8 ) {
         app_error()<<" Error during construction of 1-body propagator. Hadd is not hermitian. \n";
         app_error()<<i <<" " <<j <<" " <<Hadd(i,j) <<" " <<Hadd(j,i) <<std::endl;
         APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
       }
       v(i,j) = 0.5*(Hadd(i,j)+myconj(Hadd(j,i))); 
       v(j,i) = myconj(v(i,j)); 
     }
    }
    if(!spinRestricted) {
      for(int i=0; i<NMO; i++) {
       v(i+NMO,i) = Hadd(i+NMO,i);
       for(int j=i+1; j<NMO; j++) {
         if( std::abs( Hadd(i+NMO,j) - myconj(Hadd(j+NMO,i)) ) > 1e-8 ) {
           app_error()<<" Error during construction of 1-body propagator. Hadd is not hermitian. \n";
           app_error()<<i+NMO <<" " <<j <<" " <<Hadd(i+NMO,j) <<" " <<Hadd(j+NMO,i) <<std::endl;
           APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
         }
         v(i+NMO,j) = 0.5*(Hadd(i+NMO,j)+myconj(Hadd(j+NMO,i)));
         v(j+NMO,i) = myconj(v(i+NMO,j));
       }
      }
    }

    if(!DenseMatrixOperators::isHermitian(NMO,v.data(),NMO)) {
      app_log()<<"Hadd: \n"  <<Hadd <<std::endl;
      app_error()<<" Error during construction of 1-body propagator. Hadd is not hermitian. \n";
      APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
    }  
    if(!spinRestricted) { 
      if(!DenseMatrixOperators::isHermitian(NMO,v.data()+NMO*NMO,NMO)) {
        app_error()<<" Error during construction of 1-body propagator. Hadd(beta) is not hermitian. \n";
        APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
      }  
    }  

    // 2. Add H1
    for(s2Dit it = H1.begin(); it != H1.end(); it++) { 
      v( std::get<0>(*it), Index2Col(std::get<1>(*it)) ) += std::get<2>(*it);
      if(std::get<0>(*it) != std::get<1>(*it)) 
        v( std::get<1>(*it), Index2Col(std::get<0>(*it)) ) += myconj(std::get<2>(*it));
    }


// FIX FIX FIX: Add this contribution to Hadd on the propagator. Make the hamiltonian return this 
// during the construction of the HS potential and keep a copy with the prop. Call is v0 or something like that.
// That way you do not need to read the 2-electron intehrals if reading wfn and prop from files. 
      // adding contribution from transformation of H2 into quadratic form
      // 3.  -0.5* sum_{i,l,sigma} (sum_j <i_sigma,j_sigma|j_sigma,l_sigma> ) c+i_sigma cl_sigma 
      // NOTE: DO NOT SEE A REASON TO ADD PERMUTED CONTRIBUTION LIKE SHIWEI DOES, THEY ARE THE SAME # BY SYMMETRY 
    if(factorizedHamiltonian) {
      if(distribute_Ham) 
        APP_ABORT("Error: distribute_Ham not yet implemented in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n\n\n "); 
      for(int i=0; i<NMO; i++)
       for(int l=i; l<NMO; l++) {
         ValueType vl = ValueType(0);
         for(int j=0; j<NMO; j++) 
           vl += H(i,j,j,l);  
         v(i,l) -= 0.5*vl; 
         if(i!=l) v(l,i) -= 0.5*myconj(vl); 
         if(!spinRestricted) {
           vl=ValueType(0);
           for(int j=NMO; j<2*NMO; j++) 
             vl += H(i+NMO,j,j,l+NMO);
           v(i+NMO,Index2Col(l+NMO)) -= 0.5*vl; 
           if(i!=l) v(l+NMO,Index2Col(l+NMO)) -= 0.5*myconj(vl); 
         } 
       }
    } else { 
      std::vector<s4D<ValueType> > v2sym; 
      for(s4Dit it = V2.begin(); it != V2.end(); it++) {
        // dumb and slow, but easy for now
        find_equivalent_OneBar_for_integral_list(*it,v2sym); 
        for(int n=0; n<v2sym.size(); n++) {   
          IndexType i,j,k,l;
          ValueType V;
          std::tie (i,j,k,l,V) = v2sym[n];   
          int sector = getSpinSector(i,j,k,l);
          // <i,j | j,l> -> v(i,l)  
          if( (sector==0 || sector==3) && j==k) 
              v(i,Index2Col(l)) -= 0.5*V;
        }
      }  
    }         

      // 4. scale by -0.5*dt
      for(int i=0; i<NMO2; i++) 
       for(int j=0; j<NMO; j++) 
         v(i,j) *= -0.5*dt; 

      if(!DenseMatrixOperators::isHermitian(NMO,v.data(),NMO)) {
        app_error()<<" Error during construction of 1-body propagator. v is not hermitian. \n";
        APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
      }  

      if(!spinRestricted && !DenseMatrixOperators::isHermitian(NMO,v.data()+NMO*NMO,NMO)) {
        app_error()<<" Error during construction of 1-body propagator. v(beta) is not hermitian. \n";
        APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
      }  

      // 5. exponentiate v
      if(!DenseMatrixOperators::exponentiateHermitianMatrix(NMO,v.data(),NMO,P.data(),NMO)) {
        app_error()<<" Error during exponentiation of one body hamiltonian. \n";
        APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
      }  

      if(!spinRestricted && !DenseMatrixOperators::exponentiateHermitianMatrix(NMO,v.data()+NMO*NMO,NMO,P.data()+NMO*NMO,NMO)) {
        app_error()<<" Error during exponentiation of one body hamiltonian. \n";
        APP_ABORT("Error in SparseGeneralHamiltonian::calculateOneBodyPropagator. \n");
      }  

      // 6. get non-zero components, generate sparse representation 
      int cnt=0;
      for(int i=0; i<NMO; i++)
       for(int j=i; j<NMO; j++)  
         if(std::abs(P(i,j)) > cut) cnt = ( (i==j)?(cnt+1):(cnt+2) ); 
      if(!spinRestricted) {
        for(int i=NMO; i<2*NMO; i++)
         for(int j=i; j<2*NMO; j++)
           if(std::abs(P(i,j-NMO)) > cut) cnt = ( (i==j)?(cnt+1):(cnt+2) );
      }

      Pkin.clear();
      Pkin.reserve(cnt);
      for(IndexType i=0; i<NMO; i++)
       for(IndexType j=i; j<NMO; j++)  
         if(std::abs(P(i,j)) > cut) {
           if(i==j) {
              Pkin.push_back(std::forward_as_tuple(i,j,P(i,j)));
            } else {
              Pkin.push_back(std::forward_as_tuple(i,j,P(i,j)));
              Pkin.push_back(std::forward_as_tuple(j,i,myconj(P(i,j))));
            }
         } 
     std::sort(Pkin.begin(),Pkin.end(),mySort);
     int ntermsalpha = Pkin.size();
     if(!spinRestricted) {
        for(IndexType i=NMO; i<2*NMO; i++)
         for(IndexType j=i; j<2*NMO; j++)
           if(std::abs(P(i,j-NMO)) > cut) {
             if(i==j) {
                Pkin.push_back(std::forward_as_tuple(i-NMO,j-NMO,P(i,j-NMO)));
              } else {
                Pkin.push_back(std::forward_as_tuple(i-NMO,j-NMO,P(i,j-NMO)));
                Pkin.push_back(std::forward_as_tuple(j-NMO,i-NMO,myconj(P(i,j-NMO))));
              }
           }
     }     
     std::sort(Pkin.begin()+ntermsalpha,Pkin.end(),mySort);

     // 7. test and time if you want
     int nterms = 20;
     Timer.reset("Generic2");
     Timer.start("Generic2");

     ComplexMatrix S(NMO,NAEA),Snew(NMO,NAEA);
     for(int i=0; i<NAEA; i++) S(i,i)=1.0;

     for(int i=0; i<nterms; i++)
       DenseMatrixOperators::product(NMO,NAEA,NMO,P.data(),NMO,S.data(),NAEA,Snew.data(),NAEA); 

     Timer.stop("Generic2");
     app_log()<<" -- Average time for dense x dense MM product: " <<Timer.average("Generic2")/nterms <<"\n";

     Timer.reset("Generic2");
     Timer.start("Generic2");

     for(int i=0; i<nterms; i++)
       SparseMatrixOperators::product_SD(NAEA,Pkin.data(),ntermsalpha,S.data(),NAEA,Snew.data(),NAEA);

     Timer.stop("Generic2");
     app_log()<<" -- Average time for sparse x dense MM product: " <<Timer.average("Generic2")/nterms <<"\n";


  }

  void SparseGeneralHamiltonian::calculateHSPotentials_Diagonalization(RealType cut, RealType dt, ComplexSMSpMat& Spvn, TaskGroup& TGprop, std::vector<int>& nvec_per_node,  bool parallel)
  {

    if(factorizedHamiltonian) {
      calculateHSPotentials_FactorizedHam(cut,dt,Spvn,TGprop,nvec_per_node,parallel);
      return;
    }

    if(parallel) {
      APP_ABORT("Error: calculateHSPotentials_Diagonalization not implemented in parallel. \n");
    }

    int rnk=0;
#if defined(USE_MPI)
    rnk = rank();
#endif

      int NMO2 = NMO*NMO; 
      if(!spinRestricted) NMO2*=2;
      ValueMatrix Vuv(NMO2);
      for(int i=0; i<NMO2; i++)
       for(int j=0; j<NMO2; j++)
         Vuv(i,j)=ValueType(0.0);

      s4Dit it = V2.begin();
      std::vector<s4D<ValueType> > vs4D;
      vs4D.reserve(24);
      while(it != V2.end()) {
        find_equivalent_OneBar_for_integral_list(*it++,vs4D);
        for(int n=0; n<vs4D.size(); n++) {

          IndexType i,j,k,l,ik,lj;
          ValueType V; 
          std::tie (i,j,k,l,V) = vs4D[n];

          // do I need to check that (i,k) and (j,l) belong to the 
          // same spin sector? This would be a problem with the ME from the beginning  
          ik = i*NMO+Index2Col(k);
          lj = l*NMO+Index2Col(j);
          Vuv(ik,lj) = V;
        }
      } 

      Timer.reset("Generic");       
      Timer.start("Generic");       

      ValueMatrix eigVec(NMO2);
      RealVector eigVal(NMO2);
      if(!DenseMatrixOperators::symEigenSysAll(NMO2,Vuv.data(),NMO2,eigVal.data(),eigVec.data(),NMO2) ) {
        app_error()<<"Problems with eigenvalue/eigenstd::vector calculation in calculateHSPotentials_Diagonalization.\n";
        APP_ABORT("Problems with eigenvalue/eigenstd::vector calculation in calculateHSPotentials_Diagonalization.\n");
      }

      Timer.stop("Generic");       
      if(rnk==0) app_log()<<" -- Time to solve eigenvalue problem: " <<Timer.average("Generic") <<"\n";

      if(printEig) { 
       app_log()<<"Eigenvalues of Hamiltonian factorization: \n"; 
       for(int i=0; i<NMO2; i++) 
        if(std::abs(eigVal[i]) > std::abs(cutoff_cholesky)) 
         app_log()<<i <<" " <<eigVal[i] <<"\n";
       app_log()<<std::endl; 
      }
        

/*
      for(int i=0; i<NMO2; i++) { 
        if(eigVal[i] < 0.0) {
          app_log()<<"Found negative eigenvalue in calculateHSPotentials_Diagonalization: " <<eigVal[i] <<std::endl;
        }
      }        
*/

      int cnt1=0;
      int cnt2=0;
      for(int i=0; i<NMO2; i++) { 
       if(std::abs(eigVal[i]) > std::abs(cutoff_cholesky)) { 
         ComplexType eig = std::sqrt( ComplexType( -0.25*dt*eigVal[i],0.0) );
         int cnt3=0;
         for(int j=0; j<NMO; j++) 
          for(int k=0; k<NMO; k++) { 
            IndexType jk = j*NMO+k; 
            IndexType kj = k*NMO+j; 
            ComplexType V = eig*(eigVec(jk,i) + myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) cnt3++;
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) + myconj(eigVec(NMO+NMO+kj,i)));
              if(std::abs(V) > cut) cnt3++;
            }
          }
         if(cnt3 > 0) {
           cnt1++;
           cnt2 += cnt3;
         }

         eig = std::sqrt( ComplexType( 0.25*dt*eigVal[i],0.0) );
         cnt3=0;
         for(int j=0; j<NMO; j++) 
          for(int k=0; k<NMO; k++) { 
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ComplexType V = eig*(eigVec(jk,i) - myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) cnt3++;
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) - myconj(eigVec(NMO+NMO+kj,i)));
              if(std::abs(V) > cut) cnt3++;
            }
          }
         if(cnt3 > 0) {
           cnt1++;
           cnt2 += cnt3;
         }
       }
      }

     // later on, instead of doing all ~M^4 terms, choose a few thousand randomly
     if(test_breakup && false) {

      if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<2*NMO; i++)
       for(IndexType j=0; j<2*NMO; j++) 
        for(IndexType k=0; k<2*NMO; k++)
         for(IndexType l=0; l<2*NMO; l++) {     
           if(spinRestricted) { 
            if((i>=NMO || j>=NMO || k>=NMO || l>=NMO ))
             continue;
           } else {
            if(!goodSpinSector(i,j,k,l,NMO))
             continue;   
           }   
           ValueType v2 = H(i,j,k,l);
           ValueType v2c = 0.0;
           int kp = Index2Col(k);
           int jp = Index2Col(j);
           for(int n=0; n<NMO2; n++) if(std::abs(eigVal[n]) > std::abs(cutoff_cholesky)) v2c += eigVal[n]*eigVec(i*NMO+kp,n)*myconj(eigVec(l*NMO+jp,n));
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
      double scl = spinRestricted?(1.0):(4.0);
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/scl <<" " <<max <<" \n ********************************************\n"<<std::endl; 

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }

      Spvn.setDims(NMO2,cnt1);
      Spvn.allocate_serial(cnt2);
      
      ComplexType isqrtdt = ComplexType(0.0,std::sqrt(dt)*0.5);
      ComplexType sqrtdt = ComplexType(std::sqrt(dt)*0.5,0.0);


      cnt1=0;
      int cntn=0;
      for(int i=0; i<NMO2; i++) {
       if(std::abs(eigVal[i]) > std::abs(cutoff_cholesky)) {
         ComplexType eig = std::sqrt( ComplexType( -0.25*dt*eigVal[i],0.0) );
         int cnt3=0;
         for(int j=0; j<NMO; j++)
          for(int k=0; k<NMO; k++) {
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ComplexType V = eig*(eigVec(jk,i) + myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) {
              cnt3++;
              //V*=isqrtdt;
              Spvn.add(jk,cnt1,V);
            }
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) + myconj(eigVec(NMO*NMO+kj,i)));
              if(std::abs(V) > cut) {
                cnt3++;
                //V*=isqrtdt;
                Spvn.add(NMO*NMO+jk,cnt1,V);
              }
            }
          }
         if(cnt3 > 0) {
           cnt1++; 
         }
         eig = std::sqrt( ComplexType( 0.25*dt*eigVal[i],0.0) );
         cnt3=0;
         for(int j=0; j<NMO; j++)
          for(int k=0; k<NMO; k++) {
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ComplexType V = eig*(eigVec(jk,i) - myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) {
              cnt3++;
              //V*=sqrtdt;
              Spvn.add(jk,cnt1,V);
            }
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) - myconj(eigVec(NMO*NMO+kj,i)));
              if(std::abs(V) > cut) {
                cnt3++;
                //V*=sqrtdt;
                Spvn.add(NMO*NMO+jk,cnt1,V);
              }
            }
          }
         if(cnt3 > 0) {
           cnt1++; 
           cntn++;
         }
       }
      }

      app_log()<<"Number of positive and negative Cholesky std::vectors: " <<cnt1-cntn <<" " <<cntn <<std::endl; 

      app_log()<<"Number of HS potentials: " <<Spvn.cols() <<std::endl;
      app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;
      app_log()<<"Compressing Spvn. \n";
      Spvn.compress();
      app_log()<<"Done Compressing Spvn. \n";

     if(test_breakup) {

      if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      int* cols = Spvn.column_data();
      int* rows = Spvn.row_data();
      int* indx = Spvn.row_index();
      ComplexType* vals = Spvn.values(); 

      ComplexMatrix v2prod(NMO2);
      for(int i=0; i<NMO2; i++) 
       for(int j=0; j<NMO2; j++) 
        v2prod(i,j) = SparseMatrixOperators::product_SpVSpV<ComplexType>(indx[i+1]-indx[i],cols+indx[i],vals+indx[i],indx[j+1]-indx[j],cols+indx[j],vals+indx[j]);   

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<2*NMO; i++)
       for(IndexType j=0; j<2*NMO; j++)
        for(IndexType k=0; k<2*NMO; k++)
         for(IndexType l=0; l<2*NMO; l++) {
           if(spinRestricted) {
            if((i>=NMO || j>=NMO || k>=NMO || l>=NMO ))
             continue;
           } else {
            if(!goodSpinSector(i,j,k,l,NMO))
             continue;
           }
           ComplexType v2 = -dt*H(i,j,k,l);
           ComplexType v2c = v2prod( i*NMO+Index2Col(k) , l*NMO+Index2Col(j)   );
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
      double scl = spinRestricted?(1.0):(4.0);
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/scl <<" " <<max <<" \n ********************************************\n"<<std::endl;

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }



  }

  void SparseGeneralHamiltonian::calculateHSPotentials_FactorizedHam(RealType cut, RealType dt, ComplexSMSpMat& Spvn, TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool parallel)
  {

    if(nnodes_per_TG>1) {
      APP_ABORT("Error: calculateHSPotentials_FactorizedHam not implemented with distributed Hamiltonian (nnodes_per_TG(hamiltonian)>1) . \n");
    }

    /********************************************************************
    *  You get 2 potentials per Cholesky std::vector   
    *
    *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
    ********************************************************************/

    ComplexType isqrtdt = ComplexType(0.0,std::sqrt(dt)*0.5);
    ComplexType sqrtdt = ComplexType(std::sqrt(dt)*0.5,0.0);
    int NMO2 = NMO*NMO;
    if(!spinRestricted) NMO2*=2; 

    int n3Vect = V2_fact.cols();
    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

    // transpose V2_fact here and again at the end
    Timer.reset("Generic");
    Timer.start("Generic");
    if(head_of_nodes) 
      V2_fact.transpose(); 
    Timer.stop("Generic");
    app_log()<<" Time to transpose V2_fact inside calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

    int* cols = V2_fact.column_data();
    int* rows = V2_fact.row_data();
    int* indx = V2_fact.row_index();
    ValueType* vals = V2_fact.values();
    ValueMatrix Ln;

    if(head_of_nodes) Ln.resize((spinRestricted?(NMO):(2*NMO)) ,NMO);  

    std::vector<int> cnt_per_vec(2*n3Vect);
    int nnodes = TGprop.getNNodesPerTG();
    int cv0=0,cvN=2*n3Vect;
    std::vector<int> sets;
    std::vector<int> sz_per_node(nnodes);
    int node_number = TGprop.getLocalNodeNumber();

    int cnt=0;
    // count terms and distribute std::vectors 
    if(rank()==0) {
      // generate sparse version
      for(int n=0; n<n3Vect; n++) { 
       int np=0, nm=0;
       Ln = ValueType(0);
       for(int p=indx[n]; p<indx[n+1]; p++) {
         int i = cols[p]/NMO;
         int k = cols[p]%NMO;
         Ln(i,k) = vals[p];
       }
       // if hamiltonian is distributed, there needs to be communication here
       for(IndexType i=0; i<NMO; i++) 
        for(IndexType k=0; k<NMO; k++) { 
          // v+
          if(std::abs( (Ln(i,k) + myconj(Ln(k,i))) ) > cut) 
            np++;
          // v-
          if(std::abs( (Ln(i,k) - myconj(Ln(k,i))) ) > cut)  
            nm++;
          if(!spinRestricted) {
            if(std::abs( (Ln(i+NMO,k) + myconj(Ln(k+NMO,i))) ) > cut) 
              np++;
            if(std::abs( (Ln(i+NMO,k) - myconj(Ln(k+NMO,i))) ) > cut) 
              nm++;
          }
        }
       cnt_per_vec[2*n] = np;
       cnt_per_vec[2*n+1] = nm;
      }
    }

    if(parallel) {
      myComm->bcast(cnt_per_vec);
      int nvec = std::count_if(cnt_per_vec.begin(),cnt_per_vec.end(),
               [] (int i) { return i>0; } );
      Spvn.setDims(NMO2,nvec);


      nvec_per_node.resize(nnodes);
      if(nnodes==1) {
        cv0=0;
        cvN=2*n3Vect;
        cnt = std::accumulate(cnt_per_vec.begin(),cnt_per_vec.end(),0);
        Spvn.reserve(cnt);
        nvec_per_node[0] = nvec;
        sz_per_node[0] = cnt;
      } else {
        sets.resize(nnodes+1);
        if(rank()==0) {
          // partition std::vectors over nodes in TG
          std::vector<int> blocks(2*n3Vect+1);
          blocks[0]=0;
          cnt=0;
          for(int i=0; i<2*n3Vect; i++) {
            cnt+=cnt_per_vec[i];
            blocks[i+1] = cnt;
          }
          balance_partition_ordered_set(2*n3Vect,blocks.data(),sets);
          myComm->bcast(sets.data(),sets.size(),MPI_COMM_HEAD_OF_NODES);
          cv0 = sets[node_number];
          cvN = sets[node_number+1];
          app_log()<<" Partitioning of Cholesky Vectors: ";
          for(int i=0; i<=nnodes; i++)
            app_log()<<sets[i] <<" ";
          app_log()<<std::endl;
          app_log()<<" Number of terms in each partitioning: ";
          for(int i=0; i<nnodes; i++)
            app_log()<<accumulate(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],0) <<" ";
          app_log()<<std::endl;
          // since many std::vectors might have zero size and will be discarded below,
          // count only non-zero            
          for(int i=0; i<nnodes; i++)
            nvec_per_node[i] = std::count_if(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],
               [] (int i) { return i>0; } );
          for(int i=0; i<nnodes; i++)
            sz_per_node[i] = std::accumulate(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],0);
        } else if(head_of_nodes) {
          myComm->bcast(sets.data(),sets.size(),MPI_COMM_HEAD_OF_NODES);
          cv0 = sets[node_number];
          cvN = sets[node_number+1];
        }
        myComm->bcast(nvec_per_node);
        myComm->bcast(sz_per_node);
        Spvn.reserve(sz_per_node[node_number]);
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
      Spvn.setDims(NMO2,nvec);
      Spvn.allocate_serial(cnt);
    }            

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 
     
    if(head_of_nodes) {
      cnt=std::accumulate(nvec_per_node.begin(),nvec_per_node.begin()+node_number,0);
      for(int n=0; n<n3Vect; n++) { 
       if( cnt_per_vec[2*n]==0 && cnt_per_vec[2*n+1]==0 ) continue;
       if( !(2*n>=cv0 && 2*n<cvN) && !((2*n+1)>=cv0 && (2*n+1)<cvN)  ) continue;
       int np=0;
       Ln = ValueType(0);
       for(int p=indx[n]; p<indx[n+1]; p++) {
         int i = cols[p]/NMO;
         int k = cols[p]%NMO;
         Ln(i,k) = vals[p];
       }

       if(2*n>=cv0 && 2*n<cvN && cnt_per_vec[2*n]>0) {
         // v+
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) { 
           ComplexType V = (Ln(i,k) + myconj(Ln(k,i))); 
           if(std::abs(V) > cut) { 
             V*=isqrtdt;
             Spvn.add(i*NMO+k,cnt,V);
             ++np;
           } 
           if(!spinRestricted) {
             V = (Ln(i+NMO,k) + myconj(Ln(k+NMO,i)));
             if(std::abs(V) > cut) {
               V*=isqrtdt;
               Spvn.add(NMO*NMO+i*NMO+k,cnt,V);
               ++np;
             }
           }
         }
         ++cnt;
         if(np==0)
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n");
       }
       if((2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) {
         np=0;
         // v-
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) {
           ComplexType V = (Ln(i,k) - myconj(Ln(k,i))); 
           if(std::abs(V) > cut) {
             V*=sqrtdt;
             Spvn.add(i*NMO+k,cnt,V);
             ++np;
           }
           if(!spinRestricted) {
             V = (Ln(i+NMO,k) - myconj(Ln(k+NMO,i)));
             if(std::abs(V) > cut) {
               V*=sqrtdt;
               Spvn.add(NMO*NMO+i*NMO+k,cnt,V);
               ++np;
             }
           }
         }
         ++cnt;
         if(np==0)
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n");
       }
      }
    }

    if(rank()==0 && nnodes>1 && parallel) {
      app_log()<<" Partition of Cholesky Vectors: 0 ";
      cnt=0;
      for(int i=0; i<nnodes; i++) {
        cnt+=nvec_per_node[i];
        app_log()<<cnt <<" ";
      }
      app_log()<<std::endl;
      app_log()<<" Number of terms in Spvn per node in TG: ";
      for(int i : sz_per_node ) app_log()<<i <<" ";
      app_log()<<std::endl;
    }

    app_log()<<"Number of HS potentials: " <<Spvn.cols() <<std::endl;
    app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;

    app_log()<<"Compressing Spvn. \n";

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

    Timer.reset("Generic");
    Timer.start("Generic");
    Spvn.compress();
    Timer.stop("Generic");
    app_log()<<" Time to compress Spvn in calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;

    Timer.reset("Generic");
    Timer.start("Generic");
    if(head_of_nodes) V2_fact.transpose(); 
    Timer.stop("Generic");
    app_log()<<" Time to transpose V2_fact inside calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

  }


  void SparseGeneralHamiltonian::calculateHSPotentials(RealType cut, RealType dt, ComplexSMSpMat& Spvn, TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool parallel)
  {

    if(factorizedHamiltonian) {
      calculateHSPotentials_FactorizedHam(cut,dt,Spvn,TGprop,nvec_per_node,parallel);
      return;
    }

    int rnk=0;
    rnk = rank();
    cholesky_residuals.clear();
    cholesky_residuals.reserve(2*NMO*NMO);

    if(spinRestricted) {

      /********************************************************************
      *               Calculate Cholesky decomposition 
      *
      *   1. The mapping of the 2-el repulsion integrals to a 2-D matrix
      *      is done as follows:
      *         V(i,j,k,l) -->  V( i*NMO+k,  l*NMO+j )
      ********************************************************************/

     Timer.reset("Generic");
     Timer.start("Generic");

     int npr = myComm->size(), rk = myComm->rank(); 
     int ik0=0,ik1=NMO*NMO-1;
     if(min_i != 0 || max_i != NMO ) {
       APP_ABORT("Error: calculateHSPotentials doesn't work with distributed hamiltonian yet. \n\n\n");
     }
     if(parallel) {
       // requirements during parallelization
       // 1. ik0 >= min_i*NMO 
       // 2. ik1 <= max_i*NMO 
 
       // only works with full sparse hamiltonian, make more generic later when distribution of 
       // hamiltonian is well defined 
       int nt = (NMO*NMO)/npr, next = (NMO*NMO)%npr;
       if(rk < next) {
         ik0 = rk*(nt+1); 
         ik1 = ik0 + nt; 
       } else {
         ik0 = next*(nt+1) + (rk-next)*nt;
         ik1 = ik0 + nt - 1;
       } 
     }
     if(ik0 < min_i*NMO || ik1 > max_i*NMO) {
       APP_ABORT("Error: Problems in calculateHSPotentials. Distribution of hamiltonian is not consistent \nwith necessary distribution for parallel evaluation of Cholesky factorization. \n"); 
     }
     int nterms = ik1-ik0+1; 

     if(nterms < 1) {
       APP_ABORT("Error: Too many processors in parallel calculation of HS potential. Try reducing the number of cores, calculating in serial or reading from a file. \n\n\n ");
     }

     // will store full Cholesky std::vectors now and keep sparse versions only at the end.
     // want to avoid possible numerical issues from truncation
     std::vector< std::vector<ValueType> > L;
     L.reserve(NMO*NMO);

     std::vector<ValueType> Lnmax(NMO*NMO);
     std::vector<s2D<RealType> > IKLmax(npr); 
     s2D<RealType> mymax;
     int maxloc=0;

     // to store diagonal elements to avoid search, since they are used often 
     std::vector<ValueType> Duv(nterms);
     for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
      for(IndexType k=0; k<NMO; k++,nt++) {
        if(nt<ik0 || nt>ik1) continue;
        Duv[ik] = H(i,k,k,i);  
        if(toComplex(Duv[ik]).real() < RealType(0.0)) {
          app_log()<<" WARNING: Found negative Duv: " <<i <<" " <<k <<" " <<Duv[ik] <<std::endl;  
          if(zero_bad_diag_2eints) 
            Duv[ik] = ValueType(0);
        }
        ik++;
      }
      if(nt>ik1) break;
     }

     // D(ik,lj) = H(i,j,k,l) - sum_p Lp(ik) Lp*(lj)
     // Diagonal:  D(ik,ik) = H(i,k,k,i) - sum_p Lp(ik) Lp*(ik) 
     RealType max=0;
     IndexType ii=-1,kk=-1;
     mymax = std::make_tuple(-1,-1,0);
     for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
      for(IndexType k=0; k<NMO; k++,nt++) {
        if(nt<ik0 || nt>ik1) continue;
        if( std::abs(Duv[ik]) > max) {
          max = std::get<2>(mymax) =std::abs(Duv[ik]);  
          ii=std::get<0>(mymax)=i;
          kk=std::get<1>(mymax)=k;
        } 
        ik++;
      }
      if(nt>ik1) break;
     }
     if(ii<0 || kk<0) {
      app_error()<<"Problems with Cholesky decomposition. \n";
      APP_ABORT("Problems with Cholesky decomposition. \n");   
     }

     if(parallel) {
       myComm->allgather(reinterpret_cast<char*>(&mymax),reinterpret_cast<char*>(IKLmax.data()),sizeof(s2D<RealType>)); 
       ii=-1;kk=-1;
       max=0;
       for(int i=0; i<npr; i++) {
         if(std::get<2>(IKLmax[i])>max) {
           ii=std::get<0>(IKLmax[i]);
           kk=std::get<1>(IKLmax[i]);
           max=std::get<2>(IKLmax[i]);
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

     if(test_2eint && !parallel) {

       for(s4Dit it = V2.begin(); it != V2.end(); it++) {
         IndexType i,j,k,l;
         ValueType w1,w2,w3;
         std::tie (i,j,k,l,w1) = *it;  
 
         w2 = H(i,i,k,k);
         w3 = H(j,j,l,l);
         if( std::abs(w1) > std::sqrt(std::abs(w2*w3)) ) {
           app_log()<<" Problems with positive-definiteness: " 
                    <<i <<" " <<j <<" " <<k <<" " <<l <<" " 
                    <<w1 <<" " <<w2 <<" " <<w3 <<std::endl; 
         }

       } 

     }

     int cnt_energy_increases=0;
     RealType max_old;
     while(max > cutoff_cholesky) {

       RealType oneOverMax = 1/std::sqrt(std::abs(max));

       cholesky_residuals.push_back(std::abs(max));

       // calculate new cholesky std::vector based on (ii,kk)
       L.push_back(std::vector<ValueType>(nterms));  
       std::vector<ValueType>& Ln = L.back();
       std::vector<ValueType>::iterator it = Ln.begin();

       if(rk==maxloc) {
         for(int n=0; n<L.size()-1; n++)
           Lnmax[n] = L[n][ii*NMO+kk-ik0]; 
       }
       if(parallel && L.size()>1) {
         myComm->bcast(Lnmax.data(),L.size()-1,maxloc,myComm->getMPI()); 
       }

       Timer.start("Generic1");
       for(IndexType i=0, ik=0, nt=0; i<NMO; i++) {
         for(IndexType k=0; k<NMO; k++, nt++) {
           if(nt<ik0 || nt>ik1) continue;
           *(it++) = H(i,kk,k,ii);
         }
         if(nt>ik1) break;
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
       it = Ln.begin();
       Timer.stop("Generic1");
       
       Timer.start("Generic2");
       max_old = max;
       IndexType ii0=ii,kk0=kk;
       max=0;
       ii=-1;
       kk=-1;
       mymax = std::make_tuple(-1,-1,0);
       for(IndexType i=0,ik=0,nt=0; i<NMO; i++) {
        for(IndexType k=0; k<NMO; k++,nt++) {
         if(nt<ik0 || nt>ik1) continue;
         Duv[ik] -= Ln[ik]*myconj(Ln[ik]);  
         if(zero_bad_diag_2eints) {
           if( std::abs(Duv[ik]) > max && toComplex(Duv[ik]).real() > 0) {
             max = std::get<2>(mymax) =std::abs(Duv[ik]);  
             ii=std::get<0>(mymax)=i;
             kk=std::get<1>(mymax)=k;
           }
         } else {
           if( std::abs(Duv[ik]) > max) {
             max = std::get<2>(mymax) =std::abs(Duv[ik]);  
             ii=std::get<0>(mymax)=i;
             kk=std::get<1>(mymax)=k;
           }
         }
         ik++;
        }
        if(nt>ik1) break;
       }
       if(parallel) {
         myComm->allgather(reinterpret_cast<char*>(&mymax),reinterpret_cast<char*>(IKLmax.data()),sizeof(s2D<RealType>)); 
         ii=-1;kk=-1;
         max=0;
         for(int i=0; i<npr; i++) {
           if(std::get<2>(IKLmax[i])>max) {
             ii=std::get<0>(IKLmax[i]);
             kk=std::get<1>(IKLmax[i]);
             max=std::get<2>(IKLmax[i]);
             maxloc=i;
           }
         }
       }
       if(ii<0 || kk<0) {
        app_error()<<"Problems with Cholesky decomposition. \n";
        APP_ABORT("Problems with Cholesky decomposition. \n");
       }
       Timer.stop("Generic2");
       if(myComm->rank()==0 && printEig)
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
     app_log()<<" Found: " <<L.size() <<" Cholesky std::vectors with a cutoff of: " <<cutoff_cholesky <<std::endl;   

     Timer.stop("Generic");
     if(rnk==0) app_log()<<" -- Time to generate Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     if(test_breakup && !parallel ) {

     if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0,nt=0,ik=0; i<NMO; i++)
       for(IndexType j=0; j<NMO; j++) 
        for(IndexType k=0; k<NMO; k++)
         for(IndexType l=0; l<NMO; l++,nt++) {     
           if(nt<ik0||nt>ik1) continue;
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
       if(rnk==0) app_log()<<" -- Time to test Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     }

      /********************************************************************
      *  You get 2 potentials per Cholesky std::vector   
      *
      *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
      ********************************************************************/

      ComplexType isqrtdt = ComplexType(0.0,std::sqrt(dt)*0.5);
      ComplexType sqrtdt = ComplexType(std::sqrt(dt)*0.5,0.0);

      std::vector<ValueType> Lcomm;
      std::vector<int> cnts,displ;
      std::vector<int> cnt_per_vec(2*L.size());
      if(parallel) {
        Lcomm.resize(NMO*NMO); 
        if(rank()==0) {
          cnts.resize(npr);
          displ.resize(npr);
          int nt = (NMO*NMO)/npr, next = (NMO*NMO)%npr;
          for(int i=0; i<npr; i++) { 
            if(i < next) {
              ik0 = i*(nt+1);
              ik1 = ik0 + nt;
            } else {
              ik0 = next*(nt+1) + (i-next)*nt;
              ik1 = ik0 + nt - 1;
            }
            cnts[i]  = ik1-ik0+1;
            displ[i] = ik0;
          }
        } 
      }


      int cnt=0;
      // generate sparse version
      for(int n=0; n<L.size(); n++) { 
       ValueType* Ls; 
       if(parallel) {
         myComm->gatherv(L[n],Lcomm,cnts,displ,0);
         if(rank()==0) Ls = Lcomm.data();
       } else {
         Ls = L[n].data();
       } 
       if(rank()==0) {
         int np=0, nm=0;
         for(IndexType i=0; i<NMO; i++) 
          for(IndexType k=0; k<NMO; k++) { 
            // v+
            //if(std::abs( (L[n][i*NMO+k] + myconj(L[n][k*NMO+i])) ) > cut) {
            if(std::abs( (Ls[i*NMO+k] + myconj(Ls[k*NMO+i])) ) > cut) 
              np++;
            // v-
            //if(std::abs( (L[n][i*NMO+k] - myconj(L[n][k*NMO+i])) ) > cut) { 
            if(std::abs( (Ls[i*NMO+k] - myconj(Ls[k*NMO+i])) ) > cut)  
              nm++;
          }
         cnt_per_vec[2*n] = np;
         cnt_per_vec[2*n+1] = nm;
       } 
      }

      int nnodes = TGprop.getNNodesPerTG();
      int cv0=0,cvN=2*L.size();
      std::vector<int> sets;
      std::vector<int> sz_per_node(nnodes);
      int node_number = TGprop.getLocalNodeNumber();
      if(parallel) {
        myComm->bcast(cnt_per_vec);
        int nvec = std::count_if(cnt_per_vec.begin(),cnt_per_vec.end(),
                 [] (int i) { return i>0; } ); 
        Spvn.setDims(NMO*NMO,nvec);

        nvec_per_node.resize(nnodes);
        if(nnodes==1) {
          cv0=0;
          cvN=2*L.size();
          cnt = std::accumulate(cnt_per_vec.begin(),cnt_per_vec.end(),0);
          Spvn.reserve(cnt);
          nvec_per_node[0] = nvec; 
          sz_per_node[0] = cnt;
        } else {
          sets.resize(nnodes+1);
          if(rank()==0) {
            // partition std::vectors over nodes in TG
            std::vector<int> blocks(2*L.size()+1); 
            blocks[0]=0;
            cnt=0;
            for(int i=0; i<2*L.size(); i++) {
              cnt+=cnt_per_vec[i];
              blocks[i+1] = cnt; 
            }
            balance_partition_ordered_set(2*L.size(),blocks.data(),sets);
            myComm->bcast(sets.data(),sets.size(),MPI_COMM_HEAD_OF_NODES);
            cv0 = sets[node_number];
            cvN = sets[node_number+1];
          
            // since many std::vectors might have zero size and will be discarded below,
            // count only non-zero
            for(int i=0; i<nnodes; i++) 
              nvec_per_node[i] = std::count_if(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],
                 [] (int i) { return i>0; } );      
            for(int i=0; i<nnodes; i++) 
              sz_per_node[i] = std::accumulate(cnt_per_vec.begin()+sets[i],cnt_per_vec.begin()+sets[i+1],0);
          } else if(head_of_nodes) {
            myComm->bcast(sets.data(),sets.size(),MPI_COMM_HEAD_OF_NODES);
            cv0 = sets[node_number];
            cvN = sets[node_number+1];
          }
          myComm->bcast(nvec_per_node);  
          myComm->bcast(sz_per_node);  
          Spvn.reserve(sz_per_node[node_number]); 
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
        Spvn.setDims(NMO*NMO,nvec);
        Spvn.allocate_serial(cnt);
      } 

      if(parallel) myComm->barrier();

      cnt=std::accumulate(nvec_per_node.begin(),nvec_per_node.begin()+node_number,0);
      for(int n=0; n<L.size(); n++) { 
       if( cnt_per_vec[2*n]==0 && cnt_per_vec[2*n+1]==0 ) continue;
       ValueType* Ls;
       if(parallel) {
         myComm->gatherv(L[n],Lcomm,cnts,displ,0);
         if(head_of_nodes) {
           myComm->bcast(Lcomm.data(),Lcomm.size(),MPI_COMM_HEAD_OF_NODES);
           Ls = Lcomm.data();
         }
       } else {
         Ls = L[n].data();
       }
       if(head_of_nodes && 2*n>=cv0 && 2*n<cvN && cnt_per_vec[2*n]>0) {
         int np=0;
         // v+
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) { 
           ComplexType V = (Ls[i*NMO+k] + myconj(Ls[k*NMO+i])); 
           if(std::abs(V) > cut) { 
             V*=isqrtdt;
             Spvn.add(i*NMO+k,cnt,V);
             ++np;
           } 
         }
         ++cnt;
         if(np==0) 
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n"); 
       }
       if(head_of_nodes && (2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) {
         int np=0;
         // v-
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) { 
           ComplexType V = (Ls[i*NMO+k] - myconj(Ls[k*NMO+i]));
           if(std::abs(V) > cut) {
             V*=sqrtdt;
             Spvn.add(i*NMO+k,cnt,V);
             ++np;
           }
          }
         ++cnt;
         if(np==0) 
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n"); 
       }
      }

      if(rank()==0 && nnodes>1 && parallel) { 
        app_log()<<" Partition of Cholesky Vectors: 0 ";
        cnt=0;
        for(int i=0; i<nnodes; i++) { 
          cnt+=nvec_per_node[i]; 
          app_log()<<cnt <<" ";  
        }
        app_log()<<std::endl;
        app_log()<<" Number of terms in Spvn per node in TG: ";
        for(int i : sz_per_node ) app_log()<<i <<" "; 
        app_log()<<std::endl;
      }

      app_log()<<"Number of HS potentials: " <<Spvn.cols() <<std::endl;
      app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;
      app_log()<<"Compressing Spvn. \n";
      if(head_of_nodes) Spvn.compress();
      app_log()<<"Done Compressing Spvn. \n";
      // setup communication if parallel  

      if(parallel) myComm->barrier();

    } else {

      /********************************************************************
      *               Calculate Cholesky decomposition 
      *
      *   1. The mapping of the 2-el repulsion integrals to a 2-D matrix
      *      is done as follows:
      *         V(i,j,k,l) -->  V( i*NMO+k,  l*NMO+j )
      ********************************************************************/

     Timer.reset("Generic");
     Timer.start("Generic");

     // will store full Cholesky std::vectors now and keep sparse versions only at the end.
     // want to avoid possible numerical issues from truncation
     std::vector< std::vector<ValueType> > L;
     L.reserve(2*NMO*NMO);

     // to store diagonal elements to avoid search, since they are used often 
     ValueMatrix Duv(2*NMO,NMO);
     for(IndexType i=0; i<NMO; i++)
      for(IndexType k=i; k<NMO; k++) { 
        Duv(i,k) = H(i,k,k,i);
        if(i!=k) Duv(k,i) = Duv(i,k);
        if(i!=k) if(toComplex(Duv(i,k)).imag() > 1e-12)
         app_error()<<" Found std::complex Duv(i,k) term: " <<i <<" " <<k <<" " <<Duv(i,k) <<std::endl; 
     }

     for(IndexType i=NMO; i<2*NMO; i++)
      for(IndexType k=i; k<2*NMO; k++) { 
        Duv(i,k-NMO) = H(i,k,k,i);
        if(i!=k) Duv(k,i-NMO) = Duv(i,k-NMO);
        if(i!=k) if(toComplex(Duv(i,k-NMO)).imag() > 1e-12)
         app_error()<<" Found std::complex Duv(i,k) term: " <<i <<" " <<k <<" " <<Duv(i,k-NMO) <<std::endl; 
     }

     // D(ik,lj) = H(i,j,k,l) - sum_p Lp(ik) Lp*(lj)
     // Diagonal:  D(ik,ik) = H(i,k,k,i) - sum_p Lp(ik) Lp*(ik) 
     RealType max=0;
     IndexType ii=-1,kk=-1;
     for(IndexType i=0; i<NMO; i++)
      for(IndexType k=0; k<NMO; k++) {
       if( std::abs(Duv(i,k)) > max) {
         max = std::abs(Duv(i,k));  
         ii=i;
         kk=k;
       } 
      }
     for(IndexType i=NMO; i<2*NMO; i++)
      for(IndexType k=0; k<NMO; k++) {
       if( std::abs(Duv(i,k)) > max) {
         max = std::abs(Duv(i,k));
         ii=i;
         kk=k+NMO;
       }
      }
     if(ii<0 || kk<0) {
      app_error()<<"Problems with Cholesky decomposition. \n";
      APP_ABORT("Problems with Cholesky decomposition. \n");   
     }

     RealType max_old;
     while(max > cutoff_cholesky) {

       RealType oneOverMax = 1/std::sqrt(std::abs(max));

       // calculate new cholesky std::vector based on (ii,kk)
       L.push_back(std::vector<ValueType>(2*NMO*NMO));  
       std::vector<ValueType>& Ln = L.back();

       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) {
          ValueType& Lnik = Ln[i*NMO+k]; 
          Lnik = H(i,kk,k,ii);
          for(int n=0; n<L.size()-1; n++)
            Lnik -= L[n][i*NMO+k]*myconj(L[n][ii*NMO+Index2Col(kk)]);  
          Lnik *= oneOverMax;
        } 
       for(IndexType i=NMO; i<2*NMO; i++)
        for(IndexType k=NMO; k<2*NMO; k++) {
          ValueType& Lnik = Ln[i*NMO+Index2Col(k)];
          Lnik = H(i,kk,k,ii);
          for(int n=0; n<L.size()-1; n++)
            Lnik -= L[n][i*NMO+Index2Col(k)]*myconj(L[n][ii*NMO+Index2Col(kk)]);
          Lnik *= oneOverMax;
        }
       
       max_old = max;
       IndexType ii0=ii,kk0=kk;
       max=0;
       ii=-1;
       kk=-1;
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) {
         Duv(i,k) -= Ln[i*NMO+k]*myconj(Ln[i*NMO+k]);  
         if( std::abs(Duv(i,k)) > max) {
           max = std::abs(Duv(i,k)); 
           ii=i;
           kk=k;
         }
        }
       for(IndexType i=NMO; i<2*NMO; i++)
        for(IndexType k=0; k<NMO; k++) {
         Duv(i,k) -= Ln[i*NMO+k]*myconj(Ln[i*NMO+k]);  
         if( std::abs(Duv(i,k)) > max) {
           max = std::abs(Duv(i,k));
           ii=i;
           kk=k+NMO;
         }
        }
       if(max > max_old) {
         app_error()<<"ERROR: Problems with convergence of Cholesky decomposition. \n" 
           <<"Number of std::vectors found so far: " <<L.size() <<"\n"
           <<"Current value of truncation error: " <<max_old <<" " <<max <<std::endl;  
           APP_ABORT("Problems with convergence of Cholesky decomposition.\n"); 
       }
     }
     app_log()<<" Found: " <<L.size() <<" Cholesky std::vectors with a cutoff of: " <<cutoff_cholesky <<std::endl;   

     Timer.stop("Generic");
     if(rnk==0) app_log()<<" -- Time to generate Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     if(test_breakup && false) {

     if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<2*NMO; i++)
       for(IndexType j=0; j<2*NMO; j++) 
        for(IndexType k=0; k<2*NMO; k++)
         for(IndexType l=0; l<2*NMO; l++) {     
           if(!goodSpinSector(i,j,k,l,NMO))
            continue;
           ValueType v2 = H(i,j,k,l);
           ValueType v2c = 0.0;
           int kp = Index2Col(k);
           int jp = Index2Col(j);
           for(int n=0; n<L.size(); n++) v2c += L[n][i*NMO+kp]*myconj(L[n][l*NMO+jp]);
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
         }
      app_log()<<"\n ********************************************\n Average error due to truncated Cholesky factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/4.0 <<" " <<max <<" \n********************************************\n"<<std::endl; 

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     }

      /********************************************************************
      *  You get 2 potentials per Cholesky std::vector   
      *
      *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
      ********************************************************************/

      ComplexType isqrtdt = ComplexType(0.0,std::sqrt(dt)*0.5);
      ComplexType sqrtdt = ComplexType(std::sqrt(dt)*0.5,0.0);

      int cnt=0;
      int cntp=0, cntm=0;
      // generate sparse version
      for(int n=0; n<L.size(); n++) { 
       int np=0, nm=0;
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) { 
          // v+
          if(std::abs( (L[n][i*NMO+k] + myconj(L[n][k*NMO+i])) ) > cut) {
            cnt++; 
            np++;
          }
          if(std::abs( (L[n][NMO*NMO+i*NMO+k] + myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) { 
            cnt++; 
            np++;
          }
          // v-
          if(std::abs( (L[n][i*NMO+k] - myconj(L[n][k*NMO+i])) ) > cut) { 
            cnt++; 
            nm++;
          }
          if(std::abs( (L[n][NMO*NMO+i*NMO+k] - myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) { 
            cnt++; 
            nm++;
          }
        }
        if(np>0) ++cntp;
        if(nm>0) ++cntm;
      }

      Spvn.setDims(2*NMO*NMO,cntp+cntm);
      Spvn.allocate_serial(cnt);

      cnt=0;
      for(int n=0; n<L.size(); n++) { 
       int np=0;
       // v+
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) { 
         if(std::abs( (L[n][i*NMO+k] + myconj(L[n][k*NMO+i])) ) > cut) { 
           Spvn.add(i*NMO+k,cnt,isqrtdt*(L[n][i*NMO+k] + myconj(L[n][k*NMO+i])));
           ++np;
         }
         if(std::abs( (L[n][NMO*NMO+i*NMO+k] + myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) {
           Spvn.add(NMO*NMO+i*NMO+k,cnt,isqrtdt*(L[n][NMO*NMO+i*NMO+k] + myconj(L[n][NMO*NMO+k*NMO+i])));
           ++np;
         }
        }
       if(np>0) ++cnt;
       np=0;
       // v-
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) { 
         if(std::abs( (L[n][i*NMO+k] - myconj(L[n][k*NMO+i])) ) > cut) { 
           Spvn.add(i*NMO+k,cnt,sqrtdt*(L[n][i*NMO+k] - myconj(L[n][k*NMO+i])));
           ++np;
         }
         if(std::abs( (L[n][NMO*NMO+i*NMO+k] - myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) {
           Spvn.add(NMO*NMO+i*NMO+k,cnt,sqrtdt*(L[n][NMO*NMO+i*NMO+k] - myconj(L[n][NMO*NMO+k*NMO+i])));
           ++np;
         }
        }
       if(np>0) ++cnt;
      }

      app_log()<<"Number of HS potentials: " <<Spvn.cols() <<std::endl;
      app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;

      app_log()<<"Compressing Spvn. \n";
      Spvn.compress();
      app_log()<<"Done Compressing Spvn. \n";

    } 

     if(test_breakup) {

      if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      int* cols = Spvn.column_data();
      int* rows = Spvn.row_data();
      int* indx = Spvn.row_index();
      ComplexType* vals = Spvn.values();

      int NMO2 = spinRestricted?(NMO*NMO):(2*NMO*NMO);

      ComplexMatrix v2prod(NMO2);
      for(int i=0; i<NMO2; i++)
       for(int j=0; j<NMO2; j++)
        v2prod(i,j) = SparseMatrixOperators::product_SpVSpV<ComplexType>(indx[i+1]-indx[i],cols+indx[i],vals+indx[i],indx[j+1]-indx[j],cols+indx[j],vals+indx[j]);

      RealType s=0.0;
      RealType max=0.0;
      for(IndexType i=0; i<2*NMO; i++)
       for(IndexType j=0; j<2*NMO; j++)
        for(IndexType k=0; k<2*NMO; k++)
         for(IndexType l=0; l<2*NMO; l++) {
           if(spinRestricted) {
            if((i>=NMO || j>=NMO || k>=NMO || l>=NMO ))
             continue;
           } else {
            if(!goodSpinSector(i,j,k,l,NMO))
             continue;
           }
           ComplexType v2 = -dt*H(i,j,k,l);
           ComplexType v2c = v2prod( i*NMO+Index2Col(k) , l*NMO+Index2Col(j)   );
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
      double scl = spinRestricted?(1.0):(4.0);
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/scl <<" " <<max <<" \n ********************************************\n"<<std::endl;

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }

  }

  bool SparseGeneralHamiltonian::parse(xmlNodePtr cur)
  {

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur; 
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    std::string order("no"); 
    std::string bkp("no"); 
    std::string str1("no"); 
    std::string str2("no"); 
    std::string str3("no"); 
    std::string str4("yes"); 
    ParameterSet m_param;
    m_param.add(order,"orderStates","std::string");    
    m_param.add(cutoff1bar,"cutoff_1bar","double");
    m_param.add(cutoff2bar,"cutoff_2bar","double");
    m_param.add(cutoff_cholesky,"cutoff_decomp","double");
    m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
    m_param.add(cutoff_cholesky,"cutoff_factorization","double");
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(filename,"filename","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    m_param.add(ascii_write_file,"ascii_write_file","std::string");
    m_param.add(rotation,"rotation","std::string");
    m_param.add(nnodes_per_TG,"nnodes_per_TG","int");
    m_param.add(nnodes_per_TG,"nnodes","int");
    m_param.add(nnodes_per_TG,"nodes","int");
    m_param.add(bkp,"test_breakup","std::string");
    m_param.add(str1,"printEig","std::string");
    m_param.add(str2,"test_2eint","std::string");
    m_param.add(str3,"fix_2eint","std::string");
    m_param.add(str4,"test_algo","std::string");
    m_param.put(cur);

    orderStates=false;
    std::transform(order.begin(),order.end(),order.begin(),(int (*)(int))tolower);
    std::transform(filetype.begin(),filetype.end(),filetype.begin(),(int (*)(int))tolower);
    std::transform(bkp.begin(),bkp.end(),bkp.begin(),(int (*)(int))tolower);
    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int))tolower);
    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int))tolower);
    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int))tolower);
    std::transform(str4.begin(),str4.end(),str4.begin(),(int (*)(int))tolower);
    if(order == "yes" || order == "true") orderStates = true;  
    if(bkp == "yes" || bkp == "true") test_breakup = true;  
    if(str1 == "yes" || str1 == "true") printEig = true;  
    if(str2 == "yes" || str2 == "true") test_2eint = true;  
    if(str3 == "yes" || str3 == "true") zero_bad_diag_2eints = true;  
    if(str4 == "no" || str4 == "false") test_algo = false;  
   
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
  bool SparseGeneralHamiltonian::checkObject() 
  {
    return true;
  }

  // For a given quartet ijkl and the 3 (for real) different permutations 
  // in the reduced list, add all possible contributions to the Hamiltonian
  // taking into account cutoff and (ik)<->(jl) symmetry  
  // J1 = <ij|kl>
  // J2 = <ij|lk>
  // J3 = <ik|jl> 
  //  routine currently assumes that ijkl is the smallest of the 3 non-symmetric terms 
  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_closed_shell(OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, double cut, std::vector<s4D<ValueType> >& v)
  {
    // simple algorithm for now
    // 1. add all contributions blindly
    // 2. check for repeated and remove
    // 3. apply symmetry elimination
   
    v.reserve(24);  
    if( isComplex( J1 ) ) {
      APP_ABORT("Not yet working for std::complex \n\n\n");
    }
 
    ValueType J1J2 = 4*J1-2*J2;
    if(std::abs(J1J2) > cut) {
      push_ijkl(i,j,k,l,J1J2,v);
      push_ijkl(k,l,i,j,J1J2,v);  
    }
    ValueType J1J3 = 4*J1-2*J3;
    if(std::abs(J1J3) > cut) {
      push_ijkl(i,l,k,j,J1J3,v);  
      push_ijkl(j,k,l,i,J1J3,v);  
    }

    ValueType J2J1 = 4*J2-2*J1;
    if(std::abs(J2J1) > cut) {
      push_ijkl(i,j,l,k,J2J1,v);
      push_ijkl(k,l,j,i,J2J1,v);
    }
    ValueType J2J3 = 4*J2-2*J3;
    if(std::abs(J2J3) > cut) {
      push_ijkl(i,k,l,j,J2J3,v);
      push_ijkl(j,l,k,i,J2J3,v);
    }
       
    ValueType J3J1 = 4*J3-2*J1;
    if(std::abs(J3J1) > cut) {
      push_ijkl(i,l,j,k,J3J1,v);
      push_ijkl(k,j,l,i,J3J1,v);
    }    
    ValueType J3J2 = 4*J3-2*J2;
    if(std::abs(J3J2) > cut) {
      push_ijkl(i,k,j,l,J3J2,v);
      push_ijkl(l,j,k,i,J3J2,v);
    }    

    // order to remove consecutive repreated 
    std::sort (v.begin(), v.end(), mySort);
    // remove consecutive repeated
    std::vector<s4D<ValueType> >::iterator it;
    it = std::unique(v.begin(),v.end());
    // resize array, since unique does not remove elements
    // just reorders them  
    v.resize( std::distance(v.begin(),it) );

  }

  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_spinRestricted(OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, double cut, std::vector<s4D<ValueType> >& v) 
  { 
    // simple algorithm for now
    // 1. add all contributions blindly
    // 2. check for repeated and remove
    // 3. apply symmetry elimination
   
    v.reserve(48);  
    v.clear();
    if( isComplex( J1 ) ) {
      APP_ABORT("Not yet working for std::complex \n\n\n");
    }
  
    int i2 = i+NMO; 
    int j2 = j+NMO; 
    int k2 = k+NMO; 
    int l2 = l+NMO; 
    // 2bar terms
    ValueType J1J2 = J1-J2;
    if(std::abs(J1J2) > cut) {
      push_ijkl(i,j,k,l,J1J2,v);
      push_ijkl(k,l,i,j,J1J2,v);  
      push_ijkl(i2,j2,k2,l2,J1J2,v);
      push_ijkl(k2,l2,i2,j2,J1J2,v);  
    }
    ValueType J1J3 = J1-J3;
    if(std::abs(J1J3) > cut) {
      push_ijkl(i,l,k,j,J1J3,v);  
      push_ijkl(j,k,l,i,J1J3,v);  
      push_ijkl(i2,l2,k2,j2,J1J3,v);  
      push_ijkl(j2,k2,l2,i2,J1J3,v);  
    }

    ValueType J2J1 = J2-J1;
    if(std::abs(J2J1) > cut) {
      push_ijkl(i,j,l,k,J2J1,v);
      push_ijkl(k,l,j,i,J2J1,v);
      push_ijkl(i2,j2,l2,k2,J2J1,v);
      push_ijkl(k2,l2,j2,i2,J2J1,v);
    }
    ValueType J2J3 = J2-J3;
    if(std::abs(J2J3) > cut) {
      push_ijkl(i,k,l,j,J2J3,v);
      push_ijkl(j,l,k,i,J2J3,v);
      push_ijkl(i2,k2,l2,j2,J2J3,v);
      push_ijkl(j2,l2,k2,i2,J2J3,v);
    }
       
    ValueType J3J1 = J3-J1;
    if(std::abs(J3J1) > cut) {
      push_ijkl(i,l,j,k,J3J1,v);
      push_ijkl(k,j,l,i,J3J1,v);
      push_ijkl(i2,l2,j2,k2,J3J1,v);
      push_ijkl(k2,j2,l2,i2,J3J1,v);
    }    
    ValueType J3J2 = J3-J2;
    if(std::abs(J3J2) > cut) {
      push_ijkl(i,k,j,l,J3J2,v);
      push_ijkl(l,j,k,i,J3J2,v);
      push_ijkl(i2,k2,j2,l2,J3J2,v);
      push_ijkl(l2,j2,k2,i2,J3J2,v);
    }    

    // 1bar terms
    if(std::abs(J1) > cut) {
      push_ijkl(i,j2,k,l2,J1,v);
      push_ijkl(k,l2,i,j2,J1,v);
      push_ijkl(i,l2,k,j2,J1,v);
      push_ijkl(j,k2,l,i2,J1,v);
      push_ijkl(j,i2,l,k2,J1,v);
      push_ijkl(l,k2,j,i2,J1,v);
      push_ijkl(l,i2,j,k2,J1,v);
      push_ijkl(k,j2,i,l2,J1,v);
    }       
    if(std::abs(J2) > cut) {
      push_ijkl(i,j2,l,k2,J2,v);
      push_ijkl(k,l2,j,i2,J2,v);
      push_ijkl(i,k2,l,j2,J2,v);
      push_ijkl(j,l2,k,i2,J2,v);
      push_ijkl(j,i2,k,l2,J2,v);
      push_ijkl(l,k2,i,j2,J2,v);
      push_ijkl(k,i2,j,l2,J2,v);
      push_ijkl(l,j2,i,k2,J2,v);
    }    
    if(std::abs(J3) > cut) {
      push_ijkl(i,l2,j,k2,J3,v);
      push_ijkl(k,j2,l,i2,J3,v);
      push_ijkl(i,k2,j,l2,J3,v);
      push_ijkl(l,j2,k,i2,J3,v);
      push_ijkl(l,i2,k,j2,J3,v);
      push_ijkl(j,k2,i,l2,J3,v);
      push_ijkl(k,i2,l,j2,J3,v);
      push_ijkl(j,l2,i,k2,J3,v);
    }

    // order to remove consecutive repreated 
    std::sort (v.begin(), v.end(), mySort);
    // remove consecutive repeated
    std::vector<s4D<ValueType> >::iterator it;
    it = std::unique(v.begin(),v.end());
    // resize array, since unique does not remove elements
    // just reorders them  
    v.resize( std::distance(v.begin(),it) );

  }

  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_general(OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, double cut, std::vector<s4D<ValueType> >& v) 
  {
    APP_ABORT("Finsigh implementation. \n\n\n");
  }

  // looks for all equivalent terms associated with Vijkl 
  // Applies std::complex conjugate when needed
  // Quite inefficient, so don't use outside initialization  
  void SparseGeneralHamiltonian::find_equivalent_OneBar_for_integral_list(s4D<ValueType> Vijkl, std::vector<s4D<ValueType> >& v)
  {
    v.reserve(24);
    v.clear();
    IndexType i,j,k,l;
    ValueType V;
    std::tie (i,j,k,l,V) = Vijkl;
    if( isComplex( std::get<4>(Vijkl) ) ) {
      // only alpha/alpha sector is stored
      // so there should bever be any beta index here, should I check???
      // <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
      s3D<IndexType> ijkl = std::make_tuple(i,j,k,l); 
      s3D<IndexType> jilk = std::make_tuple(j,i,l,k); 
      s3D<IndexType> klij = std::make_tuple(k,l,i,j); 
      s3D<IndexType> lkji = std::make_tuple(l,k,j,i); 

      v.push_back(Vijkl);
      if( jilk != ijkl ) v.push_back(std::make_tuple(j,i,l,k,V));  
      if( klij != ijkl && klij != jilk ) v.push_back(std::make_tuple(k,l,i,j,std::conj(V)));   
      if( lkji != ijkl && lkji != jilk && lkji != klij ) v.push_back(std::make_tuple(l,k,j,i,std::conj(V)));   
      // order, not needed but good measure 
      std::sort (v.begin(), v.end(),mySort);
    } else {
      // only alpha/alpha sector is stored
      // so there should bever be any beta index here, should I check???
      //  For one-bar integrals: <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
      //                                 = <ji|lk> = <li|jk> = <jk|li> = <lk|ji>

      v.push_back(std::make_tuple(i,j,k,l,V));  
      v.push_back(std::make_tuple(k,j,i,l,V));  
      v.push_back(std::make_tuple(i,l,k,j,V));  
      v.push_back(std::make_tuple(k,l,i,j,V));  
      v.push_back(std::make_tuple(j,i,l,k,V));  
      v.push_back(std::make_tuple(l,i,j,k,V));  
      v.push_back(std::make_tuple(j,k,l,i,V));  
      v.push_back(std::make_tuple(l,k,j,i,V));  
      // order to remove consecutive repreated 
       

      std::sort (v.begin(), v.end(), mySort);
      //std::sort (v.begin(), v.end());
      // remove consecutive repeated
      std::vector<s4D<ValueType> >::iterator it;
      it = std::unique(v.begin(),v.end());
      // resize array, since unique does not remove elements
      // just reorders them  
      v.resize( std::distance(v.begin(),it) );
    } 

  }

  // for symmetric terms (e.g. ik/jl - jl/ik, we keep 1 and multiply V by 2.
  void SparseGeneralHamiltonian::find_equivalent_OneBar_for_hamiltonian_generation(s4D<ValueType> Vijkl, std::vector<s4D<ValueType> >& v)
  {
    v.reserve(24);
    v.clear();
    IndexType i,j,k,l;
    ValueType V;
    std::tie (i,j,k,l,V) = Vijkl;
    if( isComplex( std::get<4>(Vijkl) ) ) {
      // <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
      s3D<IndexType> ijkl = std::make_tuple(i,j,k,l); 
      s3D<IndexType> jilk = std::make_tuple(j,i,l,k); 
      s3D<IndexType> klij = std::make_tuple(k,l,i,j); 
      s3D<IndexType> lkji = std::make_tuple(l,k,j,i); 

      
      if( jilk != ijkl ) {
        v.push_back(std::make_tuple(i,j,k,l,static_cast<RealType>(2.0)*V));  
      } else 
        v.push_back(Vijkl);

      if( klij != ijkl && klij != jilk ) {

        // need to add klij as conj(V). How do I add it? 
        if(klij != lkji) {  // add once with factor of 2
          v.push_back(std::make_tuple(k,l,i,j,std::conj(static_cast<RealType>(2.0)*V)));   
          if(ijkl == lkji)  {
            std::cerr<<" Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl == lkji): " <<i <<" " <<j <<" " <<k <<" " <<l <<"  " <<V  <<std::endl;
            APP_ABORT("Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl == lkji) \n"); 
          }
        } else 
          v.push_back(std::make_tuple(k,l,i,j,std::conj(V)));

      } else {
        // just checking
        if( klij == jilk && toComplex(V).imag() > 1e-8  ) {  
          std::cerr<<" Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (klij == jilk): " <<i <<" " <<j <<" " <<k <<" " <<l <<"  " <<V  <<std::endl;
          APP_ABORT("Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (klij == jilk) \n"); 
        } 
      } 

    } else {

      v.push_back(std::make_tuple(i,j,k,l,V));  
      v.push_back(std::make_tuple(k,j,i,l,V));  
      v.push_back(std::make_tuple(i,l,k,j,V));  
      v.push_back(std::make_tuple(k,l,i,j,V));  
      v.push_back(std::make_tuple(j,i,l,k,V));  
      v.push_back(std::make_tuple(l,i,j,k,V));  
      v.push_back(std::make_tuple(j,k,l,i,V));  
      v.push_back(std::make_tuple(l,k,j,i,V));  
      // order to remove consecutive repreated 
      std::sort (v.begin(), v.end(), mySort);
      // remove consecutive repeated
      std::vector<s4D<ValueType>>::iterator it;
      it = std::unique(v.begin(),v.end());
      // resize array, since unique does not remove elements
      // just reorders them  
      v.resize( std::distance(v.begin(),it) );

      //return;

      // look for symmetric terms
      it = v.begin();
      std::vector<s4D<ValueType>>::iterator it2;
      do {
        it2 = it+1; 
        while(it2!=v.end()) {
          // <i1j1|k1l1> -> i1k1/j1l1, so look for 
          //   i1==j2, k1==l2, j1==i2, l1==k2 
          if( std::get<0>(*it)==std::get<1>(*it2) && std::get<2>(*it)==std::get<3>(*it2) && std::get<1>(*it)==std::get<0>(*it2) && std::get<3>(*it)==std::get<2>(*it2) ) { 
            it2=v.erase(it2); // since it2 > it, we can safely erase 
            std::get<4>(*it) *= static_cast<RealType>(2.0); 
          } else {
            ++it2; 
          }
        }
        it++;
      } while(it!=v.end());

    } 
  }

  void SparseGeneralHamiltonian::find_equivalent_TwoBar_for_integral_list(s4D<ValueType> Vijkl, std::vector<s4D<ValueType> >& v)
  {

    v.reserve(24);
    v.clear();
    IndexType i,j,k,l;
    ValueType V;
    std::tie (i,j,k,l,V) = Vijkl;

    // <ij||kl> -> (ijkl)
    // Symmetries:
    // (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)* 
    // This method doesn't work here, because I modify V. 
    //v.push_back(std::make_tuple(i,j,k,l,V));  
    //v.push_back(std::make_tuple(j,i,l,k,V));  
    //v.push_back(std::make_tuple(k,l,i,j,myconj(V)));  
    //v.push_back(std::make_tuple(l,k,j,i,myconj(V)));  
    //v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));  
    //v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));  
    //v.push_back(std::make_tuple(l,k,i,j,myconj(static_cast<RealType>(-1.0)*V)));  
    //v.push_back(std::make_tuple(k,l,j,i,myconj(static_cast<RealType>(-1.0)*V)));  
    
    auto ijkl = std::make_tuple(i,j,k,l); 
    auto jilk = std::make_tuple(j,i,l,k); 
    auto klij = std::make_tuple(k,l,i,j); 
    auto lkji = std::make_tuple(l,k,j,i); 
    auto ijlk = std::make_tuple(i,j,l,k); 
    auto jikl = std::make_tuple(j,i,k,l); 
    auto lkij = std::make_tuple(l,k,i,j); 
    auto klji = std::make_tuple(k,l,j,i); 

    // doing it by hand 
    v.push_back(std::make_tuple(i,j,k,l,V));  
    // slow and inefficient, but EASY to write!!!
    if( ijkl != jilk ) 
      v.push_back(std::make_tuple(j,i,l,k,V));
    if( klij != ijkl && klij != jilk ) 
      v.push_back(std::make_tuple(k,l,i,j,myconj(V)));  
    if( lkji != ijkl && lkji != jilk && lkji != klij ) 
      v.push_back(std::make_tuple(l,k,j,i,myconj(V)));  
    if( ijlk != lkji && ijlk != ijkl && ijlk != jilk && ijlk != klij ) 
      v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));  
    if( jikl != ijlk && jikl != lkji && jikl != ijkl && jikl != jilk && jikl != klij ) 
      v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));  
    if( lkij != jikl && lkij != ijlk && lkij != lkji && lkij != ijkl && lkij != jilk && lkij != klij ) 
      v.push_back(std::make_tuple(l,k,i,j,myconj(static_cast<RealType>(-1.0)*V)));  
    if( klji != lkij && klji != jikl && klji != ijlk && klji != lkji && klji != ijkl && klji != jilk && klji != klij ) 
      v.push_back(std::make_tuple(k,l,j,i,myconj(static_cast<RealType>(-1.0)*V)));  


    // just in case 
    std::sort (v.begin(), v.end(), mySort);

  }

  // to do:
  //  3. create list of 2-bar integrals based on sparse storage of smallest element per symmetry set
  void SparseGeneralHamiltonian::find_equivalent_TwoBar_for_hamiltonian_generation(s4D<ValueType> Vijkl, std::vector<s4D<ValueType> >& v)
  {

    v.reserve(24);
    v.clear();
    IndexType i,j,k,l;
    ValueType V;
    std::tie (i,j,k,l,V) = Vijkl;

    // <ij||kl> -> (ijkl)
    // Symmetries:
    // (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)* 
    // This method doesn't work here, because I modify V. 
    //v.push_back(std::make_tuple(i,j,k,l,V));  
    //v.push_back(std::make_tuple(j,i,l,k,V));  
    //v.push_back(std::make_tuple(k,l,i,j,myconj(V)));  
    //v.push_back(std::make_tuple(l,k,j,i,myconj(V)));  
    //v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));  
    //v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));  
    //v.push_back(std::make_tuple(l,k,i,j,myconj(static_cast<RealType>(-1.0)*V)));  
    //v.push_back(std::make_tuple(k,l,j,i,myconj(static_cast<RealType>(-1.0)*V)));  
    
    auto ijkl = std::make_tuple(i,j,k,l); 
    auto jilk = std::make_tuple(j,i,l,k); 
    auto klij = std::make_tuple(k,l,i,j); 
    auto lkji = std::make_tuple(l,k,j,i); 
    auto ijlk = std::make_tuple(i,j,l,k); 
    auto jikl = std::make_tuple(j,i,k,l); 
    auto lkij = std::make_tuple(l,k,i,j); 
    auto klji = std::make_tuple(k,l,j,i); 

    // slow and inefficient, but EASY to write!!!
    if( ijkl != jilk ) 
      v.push_back(std::make_tuple(i,j,k,l,static_cast<RealType>(2.0)*V));  
    else
      v.push_back(std::make_tuple(i,j,k,l,V));  

    bool t=false;
    if( klij != ijkl && klij != jilk ) { 
      t=true;
      v.push_back(std::make_tuple(k,l,i,j,myconj(V)));  
    }
    if( lkji != ijkl && lkji != jilk && lkji != klij ) { 
      if(t) 
        std::get<4>(v.back()) *= static_cast<RealType>(2.0);
      else
        v.push_back(std::make_tuple(l,k,j,i,myconj(V)));  
    }

    t=false;
    if( ijlk != lkji && ijlk != ijkl && ijlk != jilk && ijlk != klij ) {
      t=true; 
      v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));  
    }
    if( jikl != ijlk && jikl != lkji && jikl != ijkl && jikl != jilk && jikl != klij ) {
      if(t)
        std::get<4>(v.back()) *= static_cast<RealType>(2.0);
      else
        v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));  
    }

    t=false;
    if( lkij != jikl && lkij != ijlk && lkij != lkji && lkij != ijkl && lkij != jilk && lkij != klij ) { 
      t=true;
      v.push_back(std::make_tuple(l,k,i,j,myconj(static_cast<RealType>(-1.0)*V)));  
    }
    if( klji != lkij && klji != jikl && klji != ijlk && klji != lkji && klji != ijkl && klji != jilk && klji != klij ) {
      if(t) 
        std::get<4>(v.back()) *= static_cast<RealType>(2.0);
      else 
        v.push_back(std::make_tuple(k,l,j,i,myconj(static_cast<RealType>(-1.0)*V)));  
    }

    std::sort (v.begin(), v.end(), mySort);

  }

  

  // looks for the indexes of all equivalent terms in V2 and returns
  // the smaller one (based on s4D ordering). 
  s4D<ValueType> SparseGeneralHamiltonian::find_smaller_equivalent_OneBar_for_integral_list(s4D<ValueType> ijkl)
  {
    std::vector<s4D<ValueType> > v;
    find_equivalent_OneBar_for_integral_list(ijkl,v);
    std::sort (v.begin(), v.end(), mySort);
    return v[0]; 
  }

  // looks for the indexes of all equivalent terms in V2_2bar and returns
  // the smaller one (based on s4D ordering). 
  s4D<ValueType> SparseGeneralHamiltonian::find_smaller_equivalent_TwoBar_for_integral_list(s4D<ValueType> ijkl)
  {
    std::vector<s4D<ValueType> > v;
    find_equivalent_TwoBar_for_integral_list(ijkl,v);
    std::sort (v.begin(), v.end(), mySort);
    return v[0]; 
  }

  bool SparseGeneralHamiltonian::createHamiltonianForPureDeterminant( std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b, std::vector<s1D<ValueType> >& hij, std::vector<s2D<ValueType> >& Vijkl, const RealType cut, bool closed_shell)
  {

    app_error()<<" Error:  SparseGeneralHamiltonian::createHamiltonianForPureDeterminant with s2D has been disabled. \n"; 
    return false;

  }

  bool SparseGeneralHamiltonian::createHamiltonianForPureDeterminant( std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b, std::vector<s1D<ValueType> >& hij, ComplexSMSpMat& Vijkl, const RealType cut, bool closed_shell)
  {

    //  For alpha-alpha and beta-beta store two-bar integrals directly
    //  For alpha-beta and beta-alpha, store one bar integrals 
    //
    // Symmetries for real orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij> = <lk||ji>  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji> = -<lk||ij> 
    //  For one-bar integrals: <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
    //                                 = <ji|lk> = <li|jk> = <jk|li> = <lk|ji> 
    //
    // Symmetries for Complex orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij>* = <lk||ji>*  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji>* = -<lk||ij>* 
    //  For one-bar integrals: <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
    //  Notice that in this case:   <ij|kl> != <kj|il> and other permutations            
    // 

#ifdef AFQMC_DEBUG
    app_log()<<" In SparseGeneralHamiltonian::createHamiltonianForPureDeterminant." <<std::endl; 
#endif

    hij.clear();

    ValueType V;
    std::vector<s4D<ValueType> > vs4D;  
    s4D<ValueType> s;
 
    // First count how many elements we need
    int cnt1=0; 
    s2Dit it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;
      if( std::abs(V) <= cut ) continue; 
      // I can assume that i<=j
      if(spinRestricted) {
        if( i == j ) {
          if(occ_a[i]) cnt1++; 
          if(!closed_shell) if(occ_b[i+NMO]) cnt1++; 
        } else {
          if(occ_a[i]) cnt1++; 
          if(occ_a[j]) cnt1++; 
          if(!closed_shell) if(occ_b[i+NMO]) cnt1++; 
          if(!closed_shell) if(occ_b[j+NMO]) cnt1++; 
        }
      } else { 
        if( i == j ) { 
          if(occ_a[i] || occ_b[i]) cnt1++; 
        } else { 
          if(occ_a[i] || occ_b[i]) cnt1++; 
          if(occ_a[j] || occ_b[j]) cnt1++; 
    }
      }
    }  

    hij.resize(cnt1);
    s1Dit ith = hij.begin();
    it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;

      if( std::abs(V) <= cut ) continue;

      if(spinRestricted) {
        if(closed_shell) {
          if( i == j ) {  // ii alpha / ii beta 
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+i , 2*V );
          } else  {  // ij/ji both (alpha/alpha) and (beta/beta)
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+j , 2*V );
            if(occ_a[j]) *ith++ = std::make_tuple( j*NMO+i , 2*myconj(V) );
          }
        } else {
          if( i == j ) {  // ii alpha / ii beta 
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+i , V );
            if(occ_b[i+NMO]) *ith++ = std::make_tuple( NMO*NMO+i*NMO+i , V );
          } else  {  // ij/ji both (alpha/alpha) and (beta/beta)
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+j , V );
            if(occ_a[j]) *ith++ = std::make_tuple( j*NMO+i , myconj(V) );
            if(occ_b[i+NMO]) *ith++ = std::make_tuple( NMO*NMO+i*NMO+j , V );
            if(occ_b[j+NMO]) *ith++ = std::make_tuple( NMO*NMO+j*NMO+i , myconj(V) );
          }
        }
      } else {
        if( i == j ) {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j), V );
        } else {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j), V );
          if(occ_a[j] || occ_b[j]) *ith++ = std::make_tuple( Index2Mat(j,i), myconj(V) );
        }
      }
    }

    std::sort (hij.begin(), hij.end(),mySort);
    ith = std::unique(hij.begin(),hij.end(),myEqv);
    if(ith != hij.end()) {
      app_log()<<" \n\n\n*************************************************************\n"
             <<"  Error! Found repeated terms in construction of hij for Pure hamiltonian. \n"
             <<"  This should only happen if the integral file contains symmetry equivalent terms. \n"
             <<" *************************************************************\n\n\n";
      return false;
    }
    // Done with hij
    
    long cnt2=0, number_of_terms=0; 
    // add one-bar terms (mixed spin) 
    if(factorizedHamiltonian) {

      if(test_algo) { 

        app_log()<<"Testing new algorithm in createHamiltonian with factorizedHamiltonian=true. \n";

        Timer.reset("Generic");
        Timer.start("Generic");

        ValueType zero = ValueType(0);

        std::vector<s4D<ValueType> > vs4D;
        vs4D.reserve(24);
        if( isComplex(zero) ) {
          APP_ABORT("New algorithm doesn't yet work for std::complex matrix elements. \n");
        }

        cnt2=0; 
        long cnter=0, npr = myComm->size(), rk = myComm->rank();
        OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
        ValueType J1,J2,J3,fct;

        if(!spinRestricted)
          APP_ABORT("Error: Need to finish implementation for !spinRestricted. \n\n\n");

        ValueMatrix Diag(NMO,NMO); 
        for(IndexType i=0; i<NMO; i++) 
        for(IndexType k=i; k<NMO; k++, cnter++) { 
          if( cnter%npr != rk ) continue;
          Diag(i,k) =  H(i,k,i,k);
          Diag(k,i) = Diag(i,k); 
        }
        myComm->allreduce(Diag);
/*
        if( isComplex(Diag(0,0)) )
          MPI_Allreduce(MPI_IN_PLACE, Diag.data(),  2*Diag.size(), MPI_DOUBLE, MPI_SUM,
                myComm->getMPI());
        else
          MPI_Allreduce(MPI_IN_PLACE, Diag.data(),  Diag.size(), MPI_DOUBLE, MPI_SUM,
                myComm->getMPI());
*/

        bool occi, occj, occk, occl; 
        cnter=0;
        // Approximation: (similar to non factorized case)
        //   - if <ij|kl>  <  cutoff2bar, set to zero 
        for(IndexType i=0; i<NMO; i++)  {
         occi = occ_a[i]||occ_b[i];
        for(IndexType j=i; j<NMO; j++,cnter++)  {
         if( cnter%npr != rk ) continue;
         occj = occ_a[j]||occ_b[j];
        for(IndexType k=j; k<NMO; k++)  {
         occk = occ_a[k]||occ_b[k];
        for(IndexType l=k; l<NMO; l++)  {

            occl = occ_a[l]||occ_b[l];
            if( ! (occi || occj || occk || occl) ) continue; 

            J1=J2=J3=zero; 
            // J1 = <ij|kl>   
            // J1 < sqrt( <ik|ik> * <jl|jl> )  
            if( std::sqrt( std::abs(Diag(i,k)*Diag(j,l)) ) > cutoff2bar ) 
              J1 = H(i,j,k,l);  

            // J2 = <ij|lk> or <ik|lj>   
            if(i==j || l==k) {
              J2=J1; 
            } else { 
              if( std::sqrt( std::abs(Diag(i,l)*Diag(j,k)) ) > cutoff2bar ) 
                J2 = H(i,j,l,k);  
            }

            // J3 = <ik|jl> or <il|jk>  
            if(j==k) {
              J3=J1; 
            } else { 
              if( std::sqrt( std::abs(Diag(i,j)*Diag(k,l)) ) > cutoff2bar ) 
                J3 = H(i,k,j,l);  
            }

            if( std::abs(J1)<cutoff2bar && std::abs(J2)<cutoff2bar && std::abs(J3)<cutoff2bar) continue; 

            vs4D.clear();
            if(closed_shell) {
              find_all_contributions_to_hamiltonian_closed_shell(i,j,k,l,J1,J2,J3,cut,vs4D); 
            } else if(spinRestricted) {
              find_all_contributions_to_hamiltonian_spinRestricted(i,j,k,l,J1,J2,J3,cut,vs4D); 
            } else {
              find_all_contributions_to_hamiltonian_general(i,j,k,l,J1,J2,J3,cut,vs4D); 
              APP_ABORT(" Error: Finish implementation \n\n\n");
            }
            cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);                    

        }
        }
        }
        }

        myComm->allreduce(cnt2);
        number_of_terms = cnt2;
        app_log()<<" Number of terms in hamiltonian: " <<cnt2 <<std::endl;
        Vijkl.reserve(cnt2);
        cnt2=0; 

        cnter=0;
        // Approximation: (similar to non factorized case)
        //   - if <ij|kl>  <  cutoff2bar, set to zero 
        for(IndexType i=0; i<NMO; i++)  {
         occi = occ_a[i]||occ_b[i];
        for(IndexType j=i; j<NMO; j++,cnter++)  {
         if( cnter%npr != rk ) continue;
         occj = occ_a[j]||occ_b[j];
        for(IndexType k=j; k<NMO; k++)  {
         occk = occ_a[k]||occ_b[k];
        for(IndexType l=k; l<NMO; l++)  {

            occl = occ_a[l]||occ_b[l];
            if( ! (occi || occj || occk || occl) ) continue; 

            J1=J2=J3=zero; 
            // J1 = <ij|kl>   
            // J1 < sqrt( <ik|ik> * <jl|jl> )  
            if( std::sqrt( std::abs(Diag(i,k)*Diag(j,l)) ) > cutoff2bar ) 
              J1 = H(i,j,k,l);  

            // J2 = <ij|lk> or <ik|lj>   
            if(i==j || l==k) {
              J2=J1; 
            } else { 
              if( std::sqrt( std::abs(Diag(i,l)*Diag(j,k)) ) > cutoff2bar ) 
                J2 = H(i,j,l,k);  
            }

            // J3 = <ik|jl> or <il|jk>  
            if(j==k) {
              J3=J1; 
            } else { 
              if( std::sqrt( std::abs(Diag(i,j)*Diag(k,l)) ) > cutoff2bar ) 
                J3 = H(i,k,j,l);  
            }

            if( std::abs(J1)<cutoff2bar && std::abs(J2)<cutoff2bar && std::abs(J3)<cutoff2bar) continue; 

            vs4D.clear();
            if(closed_shell) {
              find_all_contributions_to_hamiltonian_closed_shell(i,j,k,l,J1,J2,J3,cut,vs4D); 
            } else if(spinRestricted) {
              find_all_contributions_to_hamiltonian_spinRestricted(i,j,k,l,J1,J2,J3,cut,vs4D); 
            } else {
              find_all_contributions_to_hamiltonian_general(i,j,k,l,J1,J2,J3,cut,vs4D); 
              APP_ABORT(" Error: Finish implementation \n\n\n");
            }
            cnt2+=add_allowed_terms(vs4D,occ_a,occ_b, Vijkl, true);                    
        }
        }
        }
        }
        myComm->barrier();
        Timer.stop("Generic");
        app_log()<<"Time to generate full Hamiltonian from factorized form: " <<Timer.total("Generic") <<std::endl; 
 
      } else {

        Timer.reset("Generic");
        Timer.start("Generic");
        // count number of terms
        cnt2 = 0;
        std::vector<s4D<ValueType> > vs4D;
        vs4D.reserve(24);
        std::vector<s4D<ValueType> > v2s4D;
        v2s4D.reserve(24);
        std::vector<s4D<ValueType> > v3s4D;
        v3s4D.reserve(24);
        // parallelize these loops later and reduce final count  
        int cnter=0, npr = myComm->size(), rk = myComm->rank(); 
        for(IndexType i=0; i<NMO; i++)  {
        for(IndexType j=0; j<NMO; j++,cnter++)  {
         if( cnter%npr != rk ) continue;  
        for(IndexType k=0; k<NMO; k++)  {
        for(IndexType l=0; l<NMO; l++)  {

          if(closed_shell) {
            // add 2*<ij|kl>-<ij|lk>          
// THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!         
// closed_shell integrals have a different symmetry 
            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl == vs4D[0]) {
              ValueType v, v1 = H(i,j,k,l);
              if(i!=j && k!=l) v = 2*v1 - H(i,j,l,k);
              else v = v1; 
              if(std::abs(v) > cut) { 
                v2s4D.clear();
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl,v2s4D);
                cnt2 += count_allowed_terms(v2s4D,occ_a,occ_a);
              }
            } 
          } else if(spinRestricted) {
            ValueType v1, v2, v3, v4;
  
            bool done=false;
            // < i j | k l> , alpha/beta and beta/alpha sectors
            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl == vs4D[0]) {
              v1 = H(i,j,k,l);
              done=true;
              if(std::abs(v1) > cut)  {  
                // transform into alpha/beta and beta/alpha sectors  
                v2s4D.clear();
                v3s4D.clear();
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl,v2s4D); 
                for(int i=0; i<v2s4D.size(); i++) {
                  v3s4D.push_back( std::forward_as_tuple(std::get<0>(v2s4D[i]),std::get<1>(v2s4D[i])+NMO,std::get<2>(v2s4D[i]),std::get<3>(v2s4D[i])+NMO,std::get<4>(v2s4D[i]))  );
                  v3s4D.push_back( std::forward_as_tuple(std::get<0>(v2s4D[i])+NMO,std::get<1>(v2s4D[i]),std::get<2>(v2s4D[i])+NMO,std::get<3>(v2s4D[i]),std::get<4>(v2s4D[i]))  );
                }
                // eliminate unnecessary terms 
                cnt2 += count_allowed_terms(v3s4D,occ_a,occ_b);
              }
            }
          
            // < i j || k l> , alpha and beta sectors  
            if(i!=j && k!=l) {

              vs4D.clear();
              find_equivalent_TwoBar_for_integral_list(ijkl,vs4D);
              std::sort (vs4D.begin(), vs4D.end(), mySort);
              if(ijkl == vs4D[0]) {
                if(done) v2 = v1 - H(i,j,l,k);
                else v2 = H_2bar(i,j,k,l);
                if(std::abs(v2) > cut)  {  
                  // transform into alpha/beta and beta/alpha sectors  
                  v2s4D.clear();
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
                  for(int i=0; i<v2s4D.size(); i++) {
                    std::get<0>(v2s4D[i]) += NMO;
                    std::get<1>(v2s4D[i]) += NMO;
                    std::get<2>(v2s4D[i]) += NMO;
                    std::get<3>(v2s4D[i]) += NMO;
                  }
                  cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
                }
              }
            }
          } else {

            // < i j | k l> , alpha/beta and beta/alpha sectors
            s4D<ValueType> ijkl_ab = std::make_tuple(i,j+NMO,k,l+NMO,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl_ab,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl_ab == vs4D[0]) {
              ValueType v1 = H(i,j+NMO,k,l+NMO);
              if(std::abs(v1) > cut)  {  
                // transform into alpha/beta and beta/alpha sectors  
                v2s4D.clear();
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl_ab,v2s4D); 
                // eliminate unnecessary terms 
                cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
              }
            }

            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            // < i j || k l> , alpha and beta sectors  
            if(i!=j && k!=l) {

              vs4D.clear();
              find_equivalent_TwoBar_for_integral_list(ijkl,vs4D);
              std::sort (vs4D.begin(), vs4D.end(), mySort);
              if(ijkl == vs4D[0]) {
                // alpha
                ValueType v1 = H_2bar(i,j,k,l);
                if(std::abs(v1) > cut)  {  
                  v2s4D.clear();
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
                }
                v1 = H_2bar(i+NMO,j+NMO,k+NMO,l+NMO);
                if(std::abs(v1) > cut)  {  
                  v2s4D.clear();
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
                }
              }
            }
          } // closed_shell/spinRestricted/etc
        }
        }
        }
        }

        std::cout.flush();
        myComm->allreduce(cnt2);
        number_of_terms = cnt2;
        Vijkl.reserve(cnt2);

        cnter = 0;
        cnt2=0;
        for(IndexType i=0; i<NMO; i++)  {
        for(IndexType j=0; j<NMO; j++, cnter++)  {
         if( cnter%npr != rk ) continue;  
        for(IndexType k=0; k<NMO; k++)  {
        for(IndexType l=0; l<NMO; l++)  {

          if(closed_shell) {
            // add 2*<ij|kl>-<ij|lk>          
      APP_ABORT("ERROR: CLosed_shell hamiltonian not working. \n\n\n");      
            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl == vs4D[0]) {
              ValueType v, v1 = H(i,j,k,l);
              if(i!=j && k!=l) v = 2*v1 - H(i,j,l,k);
              else v = v1; 
              if(std::abs(v) > cut) { 
                v2s4D.clear();
                std::get<4>(ijkl) = ValueType(v); 
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl,v2s4D);
                cnt2+=add_allowed_terms(v2s4D,occ_a,occ_a,Vijkl,true);
              }
            } 
          } else if(spinRestricted) {
            ValueType v1, v2, v3, v4;

            bool done=false;
            // < i j | k l> , alpha/beta and beta/alpha sectors
            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl == vs4D[0]) {
              v1 = H(i,j,k,l);
              done=true;
              if(std::abs(v1) > cut)  {  
                // transform into alpha/beta and beta/alpha sectors  
                v2s4D.clear();
                v3s4D.clear();
                std::get<4>(ijkl) = ValueType(v1); 
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl,v2s4D); 
                for(int i=0; i<v2s4D.size(); i++) {
                  v3s4D.push_back( std::forward_as_tuple(std::get<0>(v2s4D[i]),std::get<1>(v2s4D[i])+NMO,std::get<2>(v2s4D[i]),std::get<3>(v2s4D[i])+NMO,std::get<4>(v2s4D[i]))  );
                  v3s4D.push_back( std::forward_as_tuple(std::get<0>(v2s4D[i])+NMO,std::get<1>(v2s4D[i]),std::get<2>(v2s4D[i])+NMO,std::get<3>(v2s4D[i]),std::get<4>(v2s4D[i]))  );
                }
                // eliminate unnecessary terms 
                cnt2+=add_allowed_terms(v3s4D,occ_a,occ_b,Vijkl,true);
              }
            }
          
            // < i j || k l> , alpha and beta sectors  
            if(i!=j && k!=l) {
              std::get<4>(ijkl) = ValueType(0);
              vs4D.clear();
              find_equivalent_TwoBar_for_integral_list(ijkl,vs4D);
              std::sort (vs4D.begin(), vs4D.end(), mySort);
              if(ijkl == vs4D[0]) {
                if(done) v2 = v1 - H(i,j,l,k);
                else v2 = H_2bar(i,j,k,l);
                if(std::abs(v2) > cut)  {  
                  // transform into alpha/beta and beta/alpha sectors  
                  v2s4D.clear();
                  std::get<4>(ijkl) = ValueType(v2);
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2+=add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl,true);
                  for(int i=0; i<v2s4D.size(); i++) {
                    std::get<0>(v2s4D[i]) += NMO;
                    std::get<1>(v2s4D[i]) += NMO;
                    std::get<2>(v2s4D[i]) += NMO;
                    std::get<3>(v2s4D[i]) += NMO;
                  }
                  cnt2+=add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl,true);
                }
              }
            }
          } else {

            // < i j | k l> , alpha/beta and beta/alpha sectors
            s4D<ValueType> ijkl_ab = std::make_tuple(i,j+NMO,k,l+NMO,static_cast<ValueType>(0.0));
            vs4D.clear();
            find_equivalent_OneBar_for_integral_list(ijkl_ab,vs4D);
            std::sort (vs4D.begin(), vs4D.end(), mySort);
            if(ijkl_ab == vs4D[0]) {
              ValueType v1 = H(i,j+NMO,k,l+NMO);
              if(std::abs(v1) > cut)  {  
                // transform into alpha/beta and beta/alpha sectors  
                v2s4D.clear();
                std::get<4>(ijkl_ab) = ValueType(v1);
                find_equivalent_OneBar_for_hamiltonian_generation(ijkl_ab,v2s4D); 
                // eliminate unnecessary terms 
                cnt2+=add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl,true);
              }
            }

            s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,static_cast<ValueType>(0.0));
            // < i j || k l> , alpha and beta sectors  
            if(i!=j && k!=l) {

              vs4D.clear();
              find_equivalent_TwoBar_for_integral_list(ijkl,vs4D);
              std::sort (vs4D.begin(), vs4D.end(), mySort);
              if(ijkl == vs4D[0]) {
                // alpha
                ValueType v1 = H_2bar(i,j,k,l);
                if(std::abs(v1) > cut)  {  
                  v2s4D.clear();
                  std::get<4>(ijkl) = ValueType(v1);
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2+=add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl,true);
                }
                v1 = H_2bar(i+NMO,j+NMO,k+NMO,l+NMO);
                if(std::abs(v1) > cut)  {  
                  v2s4D.clear();
                  std::get<4>(ijkl) = ValueType(v1);
                  find_equivalent_TwoBar_for_hamiltonian_generation(ijkl,v2s4D);
                  cnt2+=add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl,true);
                }
              }
            }
          } // closed_shell/spinRestricted/etc
        }
        }
        }
        }
        myComm->barrier();
        Timer.stop("Generic");
        app_log()<<"Time to generate full Hamiltonian from factorized form: " <<Timer.total("Generic") <<std::endl; 

        myComm->allreduce(cnt2);
        app_log()<<"cnt2: " <<cnt2 <<std::endl;

      } // test_algo

      myComm->barrier();

      Timer.reset("Generic");
      Timer.start("Generic");
      if(head_of_nodes) {

        int rk,npr,ptr,n0;
        MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk); 
        MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr); 
        n0 = Vijkl.size(); // my number of terms, always from zero to n0
        ptr = n0; // position to copy elements to
        std::vector<int> size(npr);
        size[rk] = n0;
        myComm->gsum(size,MPI_COMM_HEAD_OF_NODES);
        int ntot = 0; 
        for(int i=0; i<npr; i++) ntot+=size[i];
        if(ntot != number_of_terms) { 
          app_error()<<" Problems setting up hamiltonian for wavefunction from factorized form. Inconsistent number of terms. " <<std::endl;
          return false; 
        }
        if( Vijkl.capacity() < ntot) {  
          app_error()<<" Problems setting up hamiltonian for wavefunction from factorized form. Inconsistent std::vector capacity. " <<Vijkl.capacity() <<" " <<ntot <<std::endl;
          return false; 
        }
        Vijkl.resize_serial(number_of_terms);
        for(int i=0; i<npr; i++) {
          if(i==rk) { // I send
            myComm->bcast<int>(Vijkl.row_data(),n0,i,MPI_COMM_HEAD_OF_NODES);
            myComm->bcast<int>(Vijkl.column_data(),n0,i,MPI_COMM_HEAD_OF_NODES);
            myComm->bcast(Vijkl.values(),n0,i,MPI_COMM_HEAD_OF_NODES);
          } else { // I reveive
            myComm->bcast<int>(Vijkl.row_data()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
            myComm->bcast<int>(Vijkl.column_data()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
            myComm->bcast(Vijkl.values()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
            ptr+=size[i]; 
          }
        }
        std::cout<<"rank, size: " <<rk <<" " <<Vijkl.size() <<std::endl;
      }
      myComm->barrier();
      Timer.stop("Generic");
      app_log()<<"Time to communicate Hamiltonian: " <<Timer.total("Generic") <<std::endl; 

    } else {  //factorizedHamiltonian 

// You do not need V2_2bar at all.
// Notice that for every ijkl on the list, ijlk should also be on the list 
// It should be the smallest member of its symmetry group (at least for real. 
// Think it through for std::complex. !!!
// So all you need to do is to generate the IJ array you used in the generalSD
// case and for every ijkl, search the remainig ij** sector (which you will have in
// the parallel case) for the exchange part. If it is non-zero it will be there.
// Then add all the necessary permutations (including the term ijlk-ijkl )
// This can be done in parallel with locks like the rest of the code does
// This also works on the generalSD routine, avoids V2_full 


// for parallel algorithm with distributed matrices, count the number of terms
// in the hamiltonian and divide evenly among nodes on the TG.
// Then the cores on the local node work on the local segment 

     if(test_algo) {

      app_log()<<"Testing new algorithm in createHamiltonian. \n";

      Timer.reset("Generic");
      Timer.start("Generic");

      std::vector<s4D<ValueType> > vs4D;  
      vs4D.reserve(24);
      long N = spinRestricted?NMO:2*NMO;
      if(IJ.size() == 0) {
        if(spinRestricted)
          IJ.resize(NMO*(NMO+1)/2+1);
        else
          IJ.resize(2*NMO*(2*NMO+1)/2+1);
        long nt=0;
        OrbitalType i,j,k,l;
        OrbitalType i0=std::get<0>(V2[0]);
        OrbitalType j0=std::get<1>(V2[0]);
        long p2, p1 = mapUT(i0,j0,N);  
        for(long ii=0; ii<=p1; ii++) IJ[ii]=0;
        ValueType V;
        for(s4Dit it = V2.begin(); it != V2.end(); it++, nt++) {
          std::tie (i,j,k,l,V) = *it;
          if(i!=i0 || j!=j0) {
            p2 = mapUT(i,j,N);
            for(long ii=p1+1; ii<=p2; ii++) IJ[ii]=nt;
            i0=i;
            j0=j;
            p1=p2;
          }
        }
        for(long ii=p1+1; ii<IJ.size(); ii++) IJ[ii]=nt;
      }

// right now, the algorithm will add all the terms associated with a given quartet (ijkl)
// at once. For the given ijkl, 3 (6) possible combinations can appear in the list for real (complex)  
//  Only do the smallest of the 3 (6) possible combinations to avoid duplicated 

      if( isComplex( std::get<4>(*(V2.begin())) ) ) {
        APP_ABORT("New algorithm doesn't yet work for std::complex matrix elements. \n");
      }

      ValueType zero = ValueType(0);
      cnt2=0; 
      long npr = myComm->size(), rk = myComm->rank();
      OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
      ValueType J1,J2,J3,fct;
      // to avoid race conditions, only the smallest from the symmetry inequivalent terms 
      // is treated directly  
      for(long p=0, nt=0; p<IJ.size()-1; p++) {
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

          s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,J1); 
          s4D<ValueType> ijlk = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(i,j,l,k,J1)); 
          s4D<ValueType> ikjl = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(i,k,j,l,J1)); 
          if( mySort(ijlk,ijkl) || mySort(ikjl,ijkl) ) continue; 
          J2=J3=zero; 
          // look for <ij|lk>
          if(i==j || l==k) {
            J2=J1; 
          } else { 
            long p0 = mapUT(i,std::get<1>(ijlk),N);
            s4Dit it1 = V2.begin()+IJ[p0];
            s4Dit it2 = V2.begin()+IJ[p0+1];
            s4Dit first = std::lower_bound(it1,it2,ijlk,
               [] (const s4D<ValueType>& a, const s4D<ValueType>& b) 
               {return (std::get<2>(a)<std::get<2>(b)) || 
                       (!(std::get<2>(b)<std::get<2>(a))&&(std::get<3>(a)<std::get<3>(b)));} );
            if (first!=it2 && (std::get<2>(ijlk)==std::get<2>(*first)) && (std::get<3>(ijlk)==std::get<3>(*first))) {
              J2 = std::get<4>(*first);           
            }
          }
          // look for <ik|jl>
          if(j==k) {
            J3=J1; 
          } else { 
            long p0 = mapUT(i,std::get<1>(ikjl),N);
            s4Dit it1 = V2.begin()+IJ[p0];
            s4Dit it2 = V2.begin()+IJ[p0+1];
            s4Dit first = std::lower_bound(it1,it2,ikjl,
             [] (const s4D<ValueType>& a, const s4D<ValueType>& b) 
               {return (std::get<2>(a)<std::get<2>(b)) || 
                       (!(std::get<2>(b)<std::get<2>(a))&&(std::get<3>(a)<std::get<3>(b)));} );
            if (first!=it2 && (std::get<2>(ikjl)==std::get<2>(*first)) && (std::get<3>(ikjl)==std::get<3>(*first))) {
              J3 = std::get<4>(*first); 
            }
          }

          vs4D.clear();
          if(closed_shell) {
            find_all_contributions_to_hamiltonian_closed_shell(i,j,k,l,J1,J2,J3,cut,vs4D); 
          } else if(spinRestricted) {
            find_all_contributions_to_hamiltonian_spinRestricted(i,j,k,l,J1,J2,J3,cut,vs4D); 
          } else {
            find_all_contributions_to_hamiltonian_general(i,j,k,l,J1,J2,J3,cut,vs4D); 
            APP_ABORT(" Error: Finish implementation \n\n\n");
          }
          cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);                    
        }
      }

      myComm->allreduce(cnt2);
      Vijkl.allocate(cnt2+1000);
      cnt2=0; 

      for(long p=0, nt=0; p<IJ.size()-1; p++) {
        // from n->m, I have all non-zero (k,l) for a given (i,j), with i<=j  
        long n = IJ[p];
        long m = IJ[p+1];
        if(n==m) continue; 
        nt++;  // found good term, increase counter
        if( nt%npr != rk ) continue;
        s4Dit end = V2.begin()+m; 
        for(s4Dit it = V2.begin()+n; it != end; it++) {
          // J1 = <ij|kl>   
          // J2 = <ij|lk>   
          // J3 = <ik|jl>   
          std::tie (i,j,k,l,J1) = *it;        
          s4D<ValueType> ijkl = std::make_tuple(i,j,k,l,J1); 
          // this can be made faster
          s4D<ValueType> ijlk = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(i,j,l,k,J1)); 
          // this can be made faster
          s4D<ValueType> ikjl = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(i,k,j,l,J1)); 
          if( mySort(ijlk,ijkl) || mySort(ikjl,ijkl) ) continue; 
          J2=J3=zero; 
          // look for <ij|lk>
          if(i==j || l==k) {
            J2=J1; 
          } else { 
            long p0 = mapUT(i,std::get<1>(ijlk),N);
            s4Dit it1 = V2.begin()+IJ[p0];
            s4Dit it2 = V2.begin()+IJ[p0+1];
            s4Dit first = std::lower_bound(it1,it2,ijlk,
               [] (const s4D<ValueType>& a, const s4D<ValueType>& b) 
               {return (std::get<2>(a)<std::get<2>(b)) || 
                       (!(std::get<2>(b)<std::get<2>(a))&&(std::get<3>(a)<std::get<3>(b)));} );
            if (first!=it2 && (std::get<2>(ijlk)==std::get<2>(*first)) && (std::get<3>(ijlk)==std::get<3>(*first))) {
              J2 = std::get<4>(*first);           
            }
          }
          // look for <ik|jl>
          if(j==k) {
            J3=J1; 
          } else { 
            long p0 = mapUT(i,std::get<1>(ikjl),N);
            s4Dit it1 = V2.begin()+IJ[p0];
            s4Dit it2 = V2.begin()+IJ[p0+1];
            s4Dit first = std::lower_bound(it1,it2,ikjl,
             [] (const s4D<ValueType>& a, const s4D<ValueType>& b) 
               {return (std::get<2>(a)<std::get<2>(b)) || 
                       (!(std::get<2>(b)<std::get<2>(a))&&(std::get<3>(a)<std::get<3>(b)));} );
            if (first!=it2 && (std::get<2>(ikjl)==std::get<2>(*first)) && (std::get<3>(ikjl)==std::get<3>(*first))) {
              J3 = std::get<4>(*first); 
            }
          }

          vs4D.clear();
          if(closed_shell) {
            find_all_contributions_to_hamiltonian_closed_shell(i,j,k,l,J1,J2,J3,cut,vs4D); 
          } else if(spinRestricted) {
            find_all_contributions_to_hamiltonian_spinRestricted(i,j,k,l,J1,J2,J3,cut,vs4D); 
          } else {
            find_all_contributions_to_hamiltonian_general(i,j,k,l,J1,J2,J3,cut,vs4D); 
            APP_ABORT(" Error: Finish implementation \n\n\n");
          }
          cnt2+=add_allowed_terms(vs4D,occ_a,occ_b, Vijkl, true);                    
        }
      }

      myComm->barrier();

      if(!communicate_Vijkl(Vijkl)) return false;


      Timer.stop("Generic");
      app_log()<<"Time to generate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;

  
      } else {  // test_algo

//*
      s4Dit it2 = V2.begin();
      std::vector<s4D<ValueType> > v2s4D;  
      v2s4D.reserve(24);
      if(!head_of_nodes) it2 = V2.end();
      while(it2 != V2.end()) {
       if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }
       if(spinRestricted) {

        vs4D.clear();
        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);
        if(false && closed_shell) {
//this is a possible problem, since both 1bar and 2 bar terms get combined into 1, so they will appear as repeated
//and one of tjhem will be erased later on  
          cnt2 += count_allowed_terms(vs4D,occ_a,occ_a);
        } else {
          // transform into alpha/beta and beta/alpha sectors  
          v2s4D.clear();
          for(int i=0; i<vs4D.size(); i++) {
// in principle, this is adding extra terms, it is possible to add just i,j+N,k,l+N with a factor of 2
            v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i]),std::get<1>(vs4D[i])+NMO,std::get<2>(vs4D[i]),std::get<3>(vs4D[i])+NMO,std::get<4>(vs4D[i]))  ); 
            v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i])+NMO,std::get<1>(vs4D[i]),std::get<2>(vs4D[i])+NMO,std::get<3>(vs4D[i]),std::get<4>(vs4D[i]))  ); 
          }
          // eliminate unnecessary terms 
          cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);
        }  

       } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            break;
          } 
          case 1: 
          { 
            vs4D.clear();
            find_equivalent_OneBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). I should never get here. \n"; 
            return false;  
            break;
          } 
          case 3: 
          { 
            break;
          } 
        } 
        it2++;
       }
      }

      // for now, generating V2_2bar here.
      // modify to avoid generation of V2_2bar, which is a waste of space 
      if( !generate_V2_2bar() ) {
        app_error()<<" Error generating 2bar integrals. \n";
        return false;
      }

      it2 = V2_2bar.begin();
      if(!head_of_nodes) it2 = V2_2bar.end();
      // add two-bar terms (same spin) 
      while(it2 != V2_2bar.end()) {
       if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }
       if(spinRestricted) {
        vs4D.clear();
        find_equivalent_TwoBar_for_hamiltonian_generation(*it2++,vs4D);
        if(false && closed_shell) {
          cnt2 += count_allowed_terms(vs4D,occ_a,occ_a);
        } else {
          cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
          for(int i=0; i<vs4D.size(); i++) {
            std::get<0>(vs4D[i]) += NMO;
            std::get<1>(vs4D[i]) += NMO;
            std::get<2>(vs4D[i]) += NMO;
            std::get<3>(vs4D[i]) += NMO;
          }
          cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
        }
       } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
          case 1: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). Found sector 1 term in V2_2bar. Wy???? Bug Bug Bug ???!!! \n" <<std::endl; 
            return false;  
            APP_ABORT("");
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). Found sector 1 term in V2_2bar. Wy???? Bug Bug Bug ???!!! \n" <<std::endl; 
            return false;  
            APP_ABORT("");
            break;
          } 
          case 3: 
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
        } 
        it2++;
       }
      }

      Vijkl.allocate(cnt2);
      it2 = V2.begin();
      if(!head_of_nodes) it2 = V2.end();
      while(it2 != V2.end()) {

       if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }

       if(false && closed_shell) {

        vs4D.clear();
        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);
        for(int i=0; i<vs4D.size(); i++) 
          std::get<4>(vs4D[i]) *= 2.0;
        // eliminate unnecessary terms 
        add_allowed_terms(vs4D,occ_a,occ_a,Vijkl);

       } else if(spinRestricted) {

        vs4D.clear();
        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);
        // transform into alpha/beta and beta/alpha sectors  
        v2s4D.clear();
        for(int i=0; i<vs4D.size(); i++) {
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i]),std::get<1>(vs4D[i])+NMO,std::get<2>(vs4D[i]),std::get<3>(vs4D[i])+NMO,std::get<4>(vs4D[i]))  ); 
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i])+NMO,std::get<1>(vs4D[i]),std::get<2>(vs4D[i])+NMO,std::get<3>(vs4D[i]),std::get<4>(vs4D[i]))  ); 
        }
        // eliminate unnecessary terms 
        add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl);


       } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            break;
          } 
          case 1: 
          { 
            vs4D.clear();
            find_equivalent_OneBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). I should never get here. \n"; 
            return false;  
            break;
          } 
          case 3: 
          { 
            break;
          } 
        } 
        it2++;
       }
      }

      it2 = V2_2bar.begin();
      if(!head_of_nodes) it2 = V2_2bar.end();
      while(it2 != V2_2bar.end()) {

       if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }

       if(false && closed_shell) {

        // gets alpha/alpha, add by hand beta/beta
        vs4D.clear();
        find_equivalent_TwoBar_for_hamiltonian_generation(*it2++,vs4D);
        for(int i=0; i<vs4D.size(); i++) 
          std::get<4>(vs4D[i]) *= 2.0;
        add_allowed_terms(vs4D,occ_a,occ_a,Vijkl);
      
       } else if(spinRestricted) {

        // gets alpha/alpha, add by hand beta/beta
        vs4D.clear();
        find_equivalent_TwoBar_for_hamiltonian_generation(*it2++,vs4D);
        add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
        for(int i=0; i<vs4D.size(); i++) { 
          std::get<0>(vs4D[i]) += NMO;
          std::get<1>(vs4D[i]) += NMO;
          std::get<2>(vs4D[i]) += NMO;
          std::get<3>(vs4D[i]) += NMO;
        }
        add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);

      } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
          case 1: 
          { 
            break;
          } 
          case 2: 
          { 
            break;
          } 
          case 3: 
          { 
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
        } 
        it2++;
       }
      }
     } // test_algo
// */

    } // factrorizedHamiltonnian


#ifdef AFQMC_DEBUG
    app_log()<<" Done generating sparse hamiltonians. " <<std::endl; 
    app_log()<<" Compressing sparse hamiltonians. " <<std::endl; 
#endif

    if(head_of_nodes) {

      Timer.reset("Generic");
      Timer.start("Generic");

      if(!Vijkl.remove_repeated_and_compress()) {
        APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
      }

      Timer.stop("Generic");
      app_log()<<"Time to remove_repeated_and_compress Hamiltonian: " <<Timer.total("Generic") <<std::endl;


#ifdef AFQMC_DEBUG
      app_log()<<" Done compressing sparse hamiltonians. " <<std::endl; 
#endif

      myComm->barrier();
      myComm->barrier();

    } else {
      myComm->barrier();
      if(!Vijkl.initializeChildren()) return false;
      myComm->barrier();
    }

    return true;  

  }

  bool SparseGeneralHamiltonian::communicate_Vijkl(ComplexSMSpMat& Vijkl) 
  {

    myComm->barrier();

    if(head_of_nodes) {

      int rk,npr,ptr,n0;
      MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk);
      MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr);
      n0 = Vijkl.size(); // my number of terms, always from zero to n0
      ptr = n0; // position to copy elements to 
      std::vector<int> size(npr);
      size[rk] = n0;
      myComm->gsum(size,MPI_COMM_HEAD_OF_NODES);
      int ntot = 0;
      for(int i=0; i<npr; i++) ntot+=size[i];
      if(ntot > Vijkl.capacity()) {
        app_error()<<" Problems gathering hamiltonian. Capacity of std::vector is not sufficient: " <<ntot <<" " <<Vijkl.capacity() <<" \n";
        return false;
      }
      Vijkl.resize_serial(ntot);
      for(int i=0; i<npr; i++) {
        if(i==rk) { // I send
          myComm->bcast<int>(Vijkl.row_data(),n0,i,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<int>(Vijkl.column_data(),n0,i,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast(Vijkl.values(),n0,i,MPI_COMM_HEAD_OF_NODES);
        } else { // I reveive
          myComm->bcast<int>(Vijkl.row_data()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast<int>(Vijkl.column_data()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
          myComm->bcast(Vijkl.values()+ptr,size[i],i,MPI_COMM_HEAD_OF_NODES);
          ptr+=size[i];
        }
      }
    }
    myComm->barrier();
    return true; 
  }

  bool SparseGeneralHamiltonian::createHamiltonianForPureDeterminant( std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b, std::vector<s1D<ValueType> >& hij, ComplexSpMat& Vijkl, const RealType cut, bool closed_shell)
  {

    if(factorizedHamiltonian) {
      std::cerr<<" Error: Should not be here, SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(ComplexSpMat Vijkl) has been disabled. " <<std::endl; 
      return false;
    }


// the compression is taking a very long time, possibly buggy. Check!!!! 
    std::vector<s2D<ValueType> > Vijkl_;
    createHamiltonianForPureDeterminant(occ_a,occ_b,hij,Vijkl_,cut, closed_shell);

    Vijkl.initFroms2D(Vijkl_,true);
     
    return true;

    //  For alpha-alpha and beta-beta store two-bar integrals directly
    //  For alpha-beta and beta-alpha, store one bar integrals 
    //
    // Symmetries for real orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij> = <lk||ji>  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji> = -<lk||ij> 
    //  For one-bar integrals: <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
    //                                 = <ji|lk> = <li|jk> = <jk|li> = <lk|ji> 
    //
    // Symmetries for Complex orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij>* = <lk||ji>*  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji>* = -<lk||ij>* 
    //  For one-bar integrals: <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
    //  Notice that in this case:   <ij|kl> != <kj|il> and other permutations            
    // 

#ifdef AFQMC_DEBUG
    app_log()<<" In SparseGeneralHamiltonian::createHamiltonianForPureDeterminant." <<std::endl; 
#endif

    hij.clear();
    Vijkl.clear();

    ValueType V;
    std::vector<s4D<ValueType> > vs4D;  
    s4D<ValueType> s;
 
    // First count how many elements we need
    int cnt1=0, cnt2=0; 
    s2Dit it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;
      if( std::abs(V) <= cut ) continue; 
      // I can assume that i<=j
      if(spinRestricted) {
        if( i == j ) {
          if(occ_a[i]) cnt1++; 
          if(occ_b[i+NMO]) cnt1++; 
        } else {
          if(occ_a[i]) cnt1++; 
          if(occ_a[j]) cnt1++; 
          if(occ_b[i+NMO]) cnt1++; 
          if(occ_b[j+NMO]) cnt1++; 
        }
      } else { 
        if( i == j ) { 
          if(occ_a[i] || occ_b[i]) cnt1++; 
        } else { 
          if(occ_a[i] || occ_b[i]) cnt1++; 
          if(occ_a[j] || occ_b[j]) cnt1++; 
        }
      }
    }  

    // add one-bar terms (mixed spin) 
    s4Dit it2 = V2.begin();
    std::vector<s4D<ValueType> > v2s4D;  
    v2s4D.reserve(24);
    while(it2 != V2.end()) {
      if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }
      if(spinRestricted) {

        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);
        // transform into alpha/beta and beta/alpha sectors  
        v2s4D.clear();
        for(int i=0; i<vs4D.size(); i++) {
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i]),std::get<1>(vs4D[i])+NMO,std::get<2>(vs4D[i]),std::get<3>(vs4D[i])+NMO,std::get<4>(vs4D[i]))  ); 
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i])+NMO,std::get<1>(vs4D[i]),std::get<2>(vs4D[i])+NMO,std::get<3>(vs4D[i]),std::get<4>(vs4D[i]))  ); 
        }
        // eliminate unnecessary terms 
        cnt2 += count_allowed_terms(v2s4D,occ_a,occ_b);

      } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            break;
          } 
          case 1: 
          { 
            vs4D.clear();
            find_equivalent_OneBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). I should never get here. \n"; 
            return false;  
            break;
          } 
          case 3: 
          { 
            break;
          } 
        } 
        it2++;
      }
    }

    // for now, generating V2_2bar here.
    //       // modify to avoid generation of V2_2bar, which is a waste of space 
    if( !generate_V2_2bar() ) {
      app_error()<<" Error generating 2bar integrals. \n";
      return false;
    }

    it2 = V2_2bar.begin();
    // add two-bar terms (same spin) 
    while(it2 != V2_2bar.end()) {
      if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }
      if(spinRestricted) {
        vs4D.clear();
        find_equivalent_TwoBar_for_hamiltonian_generation(*it2++,vs4D);
        cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
        for(int i=0; i<vs4D.size(); i++) {
          std::get<0>(vs4D[i]) += NMO;
          std::get<1>(vs4D[i]) += NMO;
          std::get<2>(vs4D[i]) += NMO;
          std::get<3>(vs4D[i]) += NMO;
        }
        cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
        //for(int i=0; i<vs4D.size(); i++) { 
        //  if( occ_a[get<0>(vs4D[i])] && occ_a[get<1>(vs4D[i])] ) cnt2++;
        //  if( occ_b[get<0>(vs4D[i])+NMO] && occ_b[get<1>(vs4D[i])+NMO] ) cnt2++;
        //}
      } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
          case 1: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). Found sector 1 term in V2_2bar. Wy???? Bug Bug Bug ???!!! \n" <<std::endl; 
            return false;  
            APP_ABORT("");
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). Found sector 1 term in V2_2bar. Wy???? Bug Bug Bug ???!!! \n" <<std::endl; 
            return false;  
            APP_ABORT("");
            break;
          } 
          case 3: 
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            //cnt2 += vs4D.size();
            cnt2 += count_allowed_terms(vs4D,occ_a,occ_b);
            break;
          } 
        } 
        it2++;
      }
    }

    // avoids constant resizing
    hij.reserve(cnt1);
    Vijkl.reserve(cnt2);

    s1Dit ith = hij.begin();

    it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;

      if( std::abs(V) <= cut ) continue;

      if(spinRestricted) {
        if( i == j ) {  // ii alpha / ii beta 
          if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+i , V );
          if(occ_b[i+NMO]) *ith++ = std::make_tuple( NMO*NMO+i*NMO+i , V );
        } else  {  // ij/ji both (alpha/alpha) and (beta/beta)
          if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+j , V );
          if(occ_a[j]) *ith++ = std::make_tuple( j*NMO+i , myconj(V) );
          if(occ_b[i+NMO]) *ith++ = std::make_tuple( NMO*NMO+i*NMO+j , V );
          if(occ_b[j+NMO]) *ith++ = std::make_tuple( NMO*NMO+j*NMO+i , myconj(V) );
        }
      } else {
        if( i == j ) {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j), V );
        } else {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j), V );
          if(occ_a[j] || occ_b[j]) *ith++ = std::make_tuple( Index2Mat(j,i), myconj(V) );
        }
      }
    }


    it2 = V2.begin();
    while(it2 != V2.end()) {

      if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }

      if(spinRestricted) {

        vs4D.clear();
        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);
        // transform into alpha/beta and beta/alpha sectors  
        v2s4D.clear();
        for(int i=0; i<vs4D.size(); i++) {
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i]),std::get<1>(vs4D[i])+NMO,std::get<2>(vs4D[i]),std::get<3>(vs4D[i])+NMO,std::get<4>(vs4D[i]))  ); 
          v2s4D.push_back( std::forward_as_tuple(std::get<0>(vs4D[i])+NMO,std::get<1>(vs4D[i]),std::get<2>(vs4D[i])+NMO,std::get<3>(vs4D[i]),std::get<4>(vs4D[i]))  ); 
        }
        // eliminate unnecessary terms 
        add_allowed_terms(v2s4D,occ_a,occ_b,Vijkl);


      } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            break;
          } 
          case 1: 
          { 
            find_equivalent_OneBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). I should never get here. \n"; 
            return false;  
            break;
          } 
          case 3: 
          { 
            break;
          } 
        } 

      }
    }

    it2 = V2_2bar.begin();
    while(it2 != V2_2bar.end()) {

      if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }

      if(spinRestricted) {

        // gets alpha/alpha, add by hand beta/beta
        vs4D.clear();
        vs4D.clear();
        find_equivalent_TwoBar_for_hamiltonian_generation(*it2++,vs4D);
        add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
        for(int i=0; i<vs4D.size(); i++) { 
          std::get<0>(vs4D[i]) += NMO;
          std::get<1>(vs4D[i]) += NMO;
          std::get<2>(vs4D[i]) += NMO;
          std::get<3>(vs4D[i]) += NMO;
        }
        add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);

      } else {
        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            vs4D.clear();
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
          case 1: 
          { 
            break;
          } 
          case 2: 
          { 
            break;
          } 
          case 3: 
          { 
            find_equivalent_TwoBar_for_hamiltonian_generation(*it2,vs4D);
            add_allowed_terms(vs4D,occ_a,occ_b,Vijkl);
            break;
          } 
        } 
        it2++;
      }
    }

#ifdef AFQMC_DEBUG
    app_log()<<" Done generating sparse hamiltonians. " <<std::endl; 
    app_log()<<" Compressing sparse hamiltonians. " <<std::endl; 
#endif

    std::sort (hij.begin(), hij.end(),mySort);
    ith = std::unique(hij.begin(),hij.end(),myEqv);
    if(ith != hij.end()) {
      app_log()<<" \n\n\n*************************************************************\n"
               <<"  Error! Found repeated terms in construction of hij for Pure hamiltonian. \n"
               <<"  This should only happen if the integral file contains symmetry equivalent terms. \n"
               <<" *************************************************************\n\n\n";
      return false;
    }

// fix bug on remove_repeated
    if(!Vijkl.remove_repeated()) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }
    Vijkl.compress();

#ifdef AFQMC_DEBUG
    app_log()<<" Done compressing sparse hamiltonians. " <<std::endl; 
#endif
    return true;  

  }

  bool SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant (int walker_type, const ComplexMatrix& A,std::vector<s1D<ComplexType> >& hij, ComplexSMSpMat& Vijkl, const RealType cut) 
  {

    ComplexMatrix M,N;
    const ComplexType one = ComplexType(1.0);
    const ComplexType zero = ComplexType(0.0);

    if(!spinRestricted) {
      app_error()<<" Error: SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant not implemented for UHF matrix elements. \n"; 
      return false;
    }

    if(walker_type == 0) {
      app_log()<<" Generating rotated hamiltonian matrices for RHF walker type. \n"; 
      if(!spinRestricted) {
        app_error()<<" Error: Can not generate RHF rotated hamiltonian from spin unrestricted hamiltonian. \n";
        return false;
      }
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

    // 1-body part 
    if(walker_type == 0) {

      ValueType V;
      M.resize(NMO,NMO);
      N.resize(NAEA,NMO);
      M = 0;
      N = 0;
      s2Dit it1 = H1.begin();
      while(it1 != H1.end()) {
        IndexType i,j;
        std::tie (i,j,V) = *it1++;
        M(i,j) = 2*V; //
        if( i!=j ) M(j,i) = 2*myconj(V);
      }
      DenseMatrixOperators::product_AhB(NAEA,NMO,NMO,one,A.data(),NAEA,M.data(),NMO,zero,N.data(),NMO);
      int cnt=0;
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(N(i,j)) > cut ) cnt++;
      hij.resize(cnt);
      std::vector<s1D<ComplexType> >::iterator ith = hij.begin(); 
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(N(i,j)) > cut ) *ith++ = std::make_tuple( i*NMO+j , N(i,j) ); 

    } else if(walker_type == 1) {

      M.resize(NMO,NMO);
      N.resize(NAEA,NMO);
      M = 0;
      N = 0;
      ValueType V;

      // alpha-alpha block
      s2Dit it1 = H1.begin();
      while(it1 != H1.end()) {
        IndexType i,j;
        std::tie (i,j,V) = *it1++;
        if(i<NMO && j<NMO) {
          M(i,j) = V; //
          if( i!=j ) M(j,i) = myconj(V);
        }
      }
      DenseMatrixOperators::product_AhB(NAEA,NMO,NMO,one,A.data(),A.cols(),M.data(),NMO,zero,N.data(),NMO);
      int cnt=0;
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(N(i,j)) > cut ) cnt++;

      // beta-beta block
      ComplexMatrix Q(NAEB,NMO);
      // for spinRestricted, M matrix is the same 
      if(!spinRestricted) {
        M=0;
        it1 = H1.begin();
        while(it1 != H1.end()) {
          IndexType i,j;
          std::tie (i,j,V) = *it1++;
          if(i>=NMO && j>=NMO) {
            M(i,j) = V; //
            if( i!=j ) M(j,i) = myconj(V);
          }
        }      
      }
      if( A.cols() == NAEA )
        DenseMatrixOperators::product_AhB(NAEB,NMO,NMO,one,A.data()+A.cols()*NMO,A.cols(),M.data(),NMO,zero,Q.data(),NMO);
      else if( A.cols() == NAEA+NAEB )
        DenseMatrixOperators::product_AhB(NAEB,NMO,NMO,one,A.data()+A.cols()*NMO+NAEA,A.cols(),M.data(),NMO,zero,Q.data(),NMO);
      else {
        app_error()<<" Error in rotation of hamiltonian. \n";  
        return false;
      }

      for(int i=0; i<NAEB; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(Q(i,j)) > cut ) cnt++; 

      hij.resize(cnt);
      std::vector<s1D<ComplexType> >::iterator ith = hij.begin();
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(N(i,j)) > cut ) *ith++ = std::make_tuple( i*NMO+j , N(i,j) );  
      for(int i=0; i<NAEB; i++)
       for(int j=0; j<NMO; j++)
        if( std::abs(Q(i,j)) > cut ) *ith++ = std::make_tuple( NMO*NMO+i*NMO+j , Q(i,j) );  

    } else if(walker_type == 2) {

    }

    std::sort (hij.begin(), hij.end(),mySort);
    std::vector<s1D<ComplexType> >::iterator ith = std::unique(hij.begin(),hij.end(),myEqv);
    if(ith != hij.end()) {
      app_log()<<" \n\n\n*************************************************************\n"
             <<"  Error! Found repeated terms in construction of hij for General hamiltonian. \n"
             <<"  This should never happen. \n"
             <<" *************************************************************\n\n\n";
      return false;
    }


    if(factorizedHamiltonian) {
  
      app_error()<<" Error: SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant not implemented with factorizedHamiltonian. \n";
      return false;  

    } else {


        // Alg:
        // 0. Erase V2. Create full sparse hamiltonian (with symmetric terms)  
        // 1. swap indexes, (i,j) <-> (k,l)
        // 2. sort
        // 3. use A*M*A algorithm
        

        Timer.reset("Generic");
        Timer.reset("Generic1");
        Timer.reset("Generic2");
        //Timer.reset("Generic3");
        Timer.start("Generic");

        std::vector<s4D<ValueType> > v2sym;
        v2sym.reserve(24);
        
        // to improve memory usage, split work based on index
        // and only generate the (i,j) sectors that you need locally 
        if(V2_full.size()==0) {
        
          app_log()<<" Generating V2_full in SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant() \n";
          Timer.start("Generic1"); 
          long cnter=0;
          if(head_of_nodes) {
            for(s4Dit it = V2.begin(); it != V2.end(); it++) {
              find_equivalent_OneBar_for_integral_list(*it,v2sym);
              for(int n=0; n<v2sym.size(); n++) {
                IndexType i0,j0,k0,l0;
                ValueType V;
                std::tie (i0,j0,k0,l0,V) = v2sym[n];
                if(k0<l0) continue;
                cnter++;
              }
            }
          }

          V2_full.reserve(cnter);

          if(head_of_nodes) {
            for(s4Dit it = V2.begin(); it != V2.end(); it++) {
              find_equivalent_OneBar_for_integral_list(*it,v2sym);
              for(int n=0; n<v2sym.size(); n++) {
                IndexType i0,j0,k0,l0;
                ValueType V;
                std::tie (i0,j0,k0,l0,V) = v2sym[n];
                if(k0<l0) continue;
                V2_full.push_back(v2sym[n]);
                std::swap( std::get<0>( V2_full.back() ), std::get<2>( V2_full.back() ) );
                std::swap( std::get<1>( V2_full.back() ), std::get<3>( V2_full.back() ) );
              }
            }
            std::sort (V2_full.begin(), V2_full.end(),mySort);
          }
          v2full_transposed=true;
          myComm->barrier();
          IndexType k0=std::get<0>(V2_full[0]); 
          IndexType l0=std::get<1>(V2_full[0]); 
          IndexType i,j,k,l;
          ValueType V;
          if(spinRestricted)
            KL.reserve(NMO*(NMO+1)/2+1);
          else
            KL.reserve(2*NMO*(2*NMO+1)/2+1);
          KL.push_back(0);
          long nt=0;
          for(s4Dit it = V2_full.begin(); it != V2_full.end(); it++, nt++) {
            std::tie (k,l,i,j,V) = *it;
            if(k!=k0 || l!=l0) {
              KL.push_back(nt);
              k0=k;
              l0=l;
            }
          }
          KL.push_back(nt);
          Timer.stop("Generic1"); 
          app_log()<<" Time to generate V2_full: " <<Timer.total("Generic1") <<std::endl; 
          app_log()<<" Size of V2_full: " <<V2_full.size()*sizeof(s4D<ValueType>)/1024.0/1024.0 <<" MB. " <<std::endl;

        }
        if(!v2full_transposed) {
          Timer.reset("Generic1");
          Timer.start("Generic1");
          if(head_of_nodes) { 
            //swap         
            for(s4Dit it = V2_full.begin(); it != V2_full.end(); it++) { 
              std::swap( std::get<0>(*it), std::get<2>(*it));  
              std::swap( std::get<1>(*it), std::get<3>(*it));  
            }
            std::sort (V2_full.begin(), V2_full.end(),mySort); 
          }
          myComm->barrier();
          IndexType k0=std::get<0>(V2_full[0]);
          IndexType l0=std::get<1>(V2_full[0]);
          IndexType i,j,k,l;
          ValueType V;
          KL.clear();
          KL.reserve(2*NMO*(2*NMO+1)/2+1);
          KL.push_back(0);
          long nt=0;
          for(s4Dit it = V2_full.begin(); it != V2_full.end(); it++, nt++) {
            std::tie (k,l,i,j,V) = *it;
            if(k!=k0 || l!=l0) {
              KL.push_back(nt);
              k0=k;
              l0=l;
            }
          }
          KL.push_back(nt);
          v2full_transposed=true;
          Timer.stop("Generic1");
          app_log()<<" Time to transpose and sort V2_full: " <<Timer.total("Generic1") <<std::endl; 
          app_log()<<" Size of V2_full: " <<V2_full.size()*sizeof(s4D<ValueType>)/1024.0/1024.0 <<" MB. " <<std::endl;
        }

      int nrows = (walker_type==0)?NMO:2*NMO;
      int ncols = (walker_type==2)?(NAEA+NAEB):NAEA;
      ComplexMatrix Ac(nrows,ncols);
      // Ac = myconj(A) 
      for(int i=0; i<nrows; i++)    
       for(int j=0; j<ncols; j++)    
        Ac(i,j) = myconj(A(i,j));
      if(NAEA>NAEB && walker_type==1) {
        // fill empty spaces with 0
        for(int i=NMO; i<2*NMO; i++)
         for(int j=NAEB; j<NAEA; j++)
          Ac(i,j) = ComplexType(0); 
      }
      int npr = myComm->size(), rk = myComm->rank();
      double fct=1.0, fct2=1.0;        
      ComplexMatrix::iterator itQa,itQb,itQab, itQab2, itAia, itAib, itAja, itAjb; 
      ComplexMatrix::iterator itBia, itBib, itBja, itBjb; 
      ComplexType s1, s2, s3, s4, s5, s6;
      long cnter = 0;
      ComplexMatrix Qa,Qb,Qab,Qab2;

      if(walker_type == 0) {
  
        Qa.resize(NAEA,NAEA);
        for(int p=0, nt=0; p<KL.size()-1; p++, nt++) {
          if( nt%npr != rk ) continue;
          long n = KL[p];
          long m = KL[p+1];
          IndexType k0 = std::get<0>(V2_full[n]);
          IndexType l0 = std::get<1>(V2_full[n]);
          IndexType i,j,k,l;
          ValueType V,fct; 
          Qa=0;
          // from n->m, I have all non-zero (i,j) for a given (k,l), with k<=l
          for(s4Dit it = V2_full.begin()+n; it != V2_full.begin()+m; it++) {
            std::tie (k,l,i,j,V) = *it;
            // V (a,b,k,l) = sum_i,j A*(i,a) * A*(j,b) * V(i,j,k,l)
            if(k0==l0) fct=1.0;
            else fct=2.0;
            itAia = Ac.begin(i);
            itAja = Ac.begin(j);
            for(int a=0; a<NAEA; a++, itAia++, itAja++) {
             s1 = fct*4*V*(*itAia);
             s2 = fct*2*V*(*itAja);
             itQa = Qa.begin(a);
             itAib = Ac.begin(i);
             itAjb = Ac.begin(j);          
             for(int b=0; b<NAEA; b++, itQa++, itAib++, itAjb++) 
               *itQa += ( s1*(*itAjb) - s2*(*itAib) );
            } 
          } 
          // Now I have Q(a,b,k0,l0);
          for(int a=0; a<NAEA; a++)
           for(int b=0; b<NAEA; b++) 
             if(std::abs(Qa(a,b)) > cut) cnter++; 
        }

      } else if(walker_type == 1) {

        Qa.resize(NAEA,NAEA);  // Vaa(a,b,k,l)   k<=l
        Qb.resize(NAEA,NAEA);  // Vbb(a,b,k,l)   k<=l
        Qab.resize(NAEA,NAEA); // Vab(a,b,k,l)   k<=l
        Qab2.resize(NAEA,NAEA); // Vab(a,b,l,k)  k<=l
        int del = (walker_type==2)?NAEA:0;
        for(int p=0, nt=0; p<KL.size()-1; p++, nt++) {
          if( nt%npr != rk ) continue;
          long n = KL[p];
          long m = KL[p+1];
          IndexType k0 = std::get<0>(V2_full[n]);
          IndexType l0 = std::get<1>(V2_full[n]);
          IndexType i,j,k,l;
          ValueType V,fct; 
          if(spinRestricted) {      
            Qa=0;
            Qb=0;
            Qab=0;
            Qab2=0;
            // from n->m, I have all non-zero (i,j) for a given (k,l), with k<=l
            for(s4Dit it = V2_full.begin()+n; it != V2_full.begin()+m; it++) {
              std::tie (k,l,i,j,V) = *it;
              // V (a,b,k,l) = sum_i,j A*(i,a) * A*(j,b) * V(i,j,k,l)
              if(k0==l0) fct=0.0;
              else fct=2.0;
              itAia = Ac.begin(i);
              itBia = Ac.begin(i+NMO)+del;
              itAja = Ac.begin(j);
              itBja = Ac.begin(j+NMO)+del;
              for(int a=0; a<NAEA; a++, itAia++, itAja++, itBia++, itBja++) {
                s1 = fct*V*(*itAia);  // Coul aa
                s2 = fct*V*(*itAja);  // Exch aa
                s3 = fct*V*(*itBia);  // Coul bb
                s4 = fct*V*(*itBja);  // Exch bb
                s5 = 2.0*V*(*itAia);  // Coul ab (kl) 
                s6 = 2.0*V*(*itAja);  // Coul ab (lk) 
                itQa = Qa.begin(a);
                itQab = Qab.begin(a);
                itQab2 = Qab2.begin(a);
                itQb = Qb.begin(a);
                itAib = Ac.begin(i);
                itAjb = Ac.begin(j);          
                itBib = Ac.begin(i+NMO)+del;
                itBjb = Ac.begin(j+NMO)+del;          
                for(int b=0; b<NAEA; b++, itQa++, itQb++, itQab++, itQab2++, itAib++, itAjb++, itBib++, itBjb++) { 
                  *itQa += ( s1*(*itAjb) - s2*(*itAib) );
                  *itQb += ( s3*(*itBjb) - s4*(*itBib) );
                  *itQab += s5*(*itBjb);
                  *itQab2 += s6*(*itBib);
                }
              }
            }
            // Now I have Q(a,b,k0,l0);
            for(int a=0; a<NAEA; a++)
             for(int b=0; b<NAEA; b++) { 
               if(std::abs(Qa(a,b)) > cut) cnter++; 
               if(std::abs(Qb(a,b)) > cut) cnter++; 
               if(std::abs(Qab(a,b)) > cut) cnter++; 
               if(std::abs(Qab2(a,b)) > cut && k0!=l0) cnter++; 
             }

          } else {
            app_error()<<" Error: Finish implementation. \n\n\n";
            return false;
            if(k0 < NMO && l0 < NMO) { // aa 
              for(int a=0; a<NAEA; a++, itAia++, itAja++) {
               s1 = fct*4*V*(*itAia);
               s2 = fct*2*V*(*itAja);
               itQa = Qa.begin(a);
               itAib = Ac.begin(i);
               itAjb = Ac.begin(j);
               for(int b=0; b<NAEA; b++, itQa++, itAib++, itAjb++)
                 *itQa += ( s1*(*itAjb) - s2*(*itAib) );
              }
            } else if(k0 >= NMO && l0 >= NMO) { // bb 
               app_error()<<" Error: Finish implementation. \n\n\n";
               return false;
            } else if(k0 < NMO && l0 >= NMO) { // ab
               app_error()<<" Error: Finish implementation. \n\n\n";
               return false;
            } else {
               app_error()<<" Error: Problems in createHamiltonianForGeneralDeterminant. \n";
               return false;
            }
          } // spinRestricted
        } // KL loop

      } else if(walker_type == 2) {
      }

      myComm->allreduce(cnter);
      Vijkl.allocate(cnter+1000);

      if(walker_type == 0) {

        for(int p=0, nt=0; p<KL.size()-1; p++, nt++) {
          if( nt%npr != rk ) continue;
          long n = KL[p];
          long m = KL[p+1];
          IndexType k0 = std::get<0>(V2_full[n]);
          IndexType l0 = std::get<1>(V2_full[n]);
          IndexType i,j,k,l;
          ValueType V,fct; 
          Qa=0;
          // from n->m, I have all non-zero (i,j) for a given (k,l), with k<=l
          for(s4Dit it = V2_full.begin()+n; it != V2_full.begin()+m; it++) {
            std::tie (k,l,i,j,V) = *it;
            // V (a,b,k,l) = sum_i,j A*(i,a) * A*(j,b) * V(i,j,k,l)
            if(k0==l0) fct=1.0;
            else fct=2.0;
            itAia = Ac.begin(i);
            itAja = Ac.begin(j);
            for(int a=0; a<NAEA; a++, itAia++, itAja++) {
             s1 = fct*4*V*(*itAia);
             s2 = fct*2*V*(*itAja);
             itQa = Qa.begin(a);
             itAib = Ac.begin(i);
             itAjb = Ac.begin(j);
             for(int b=0; b<NAEA; b++, itQa++, itAib++, itAjb++) 
               *itQa += ( s1*(*itAjb) - s2*(*itAib) );
            }
          } 
          // Now I have Q(a,b,k0,l0);
          for(int a=0; a<NAEA; a++)
           for(int b=0; b<NAEA; b++) 
            if(std::abs(Qa(a,b)) > cut) 
             Vijkl.add( a*NMO+k0, b*NMO+l0, Qa(a,b), true);
        }

      } else if(walker_type == 1) {

        int del = (walker_type==2)?NAEA:0;
        for(int p=0, nt=0; p<KL.size()-1; p++, nt++) {
          if( nt%npr != rk ) continue;
          long n = KL[p];
          long m = KL[p+1];
          IndexType k0 = std::get<0>(V2_full[n]);
          IndexType l0 = std::get<1>(V2_full[n]);
          IndexType i,j,k,l;
          ValueType V,fct; 
          if(spinRestricted) {      
            Qa=0;
            Qb=0;
            Qab=0;
            Qab2=0;
            // from n->m, I have all non-zero (i,j) for a given (k,l), with k<=l
            for(s4Dit it = V2_full.begin()+n; it != V2_full.begin()+m; it++) {
              std::tie (k,l,i,j,V) = *it;
              // V (a,b,k,l) = sum_i,j A*(i,a) * A*(j,b) * V(i,j,k,l)
              if(k0==l0) fct=0.0;
              else fct=2.0;
              itAia = Ac.begin(i);
              itBia = Ac.begin(i+NMO)+del;
              itAja = Ac.begin(j);
              itBja = Ac.begin(j+NMO)+del;
              for(int a=0; a<NAEA; a++, itAia++, itAja++, itBia++, itBja++) {
                s1 = fct*V*(*itAia);  // Coul aa
                s2 = fct*V*(*itAja);  // Exch aa
                s3 = fct*V*(*itBia);  // Coul bb
                s4 = fct*V*(*itBja);  // Exch bb
                s5 = 2.0*V*(*itAia);  // Coul ab (kl) 
                s6 = 2.0*V*(*itAja);  // Coul ab (lk) 
                itQa = Qa.begin(a);
                itQab = Qab.begin(a);
                itQab2 = Qab2.begin(a);
                itQb = Qb.begin(a);
                itAib = Ac.begin(i);
                itAjb = Ac.begin(j);          
                itBib = Ac.begin(i+NMO)+del;
                itBjb = Ac.begin(j+NMO)+del;          
                for(int b=0; b<NAEA; b++, itQa++, itQb++, itQab++, itQab2++, itAib++, itAjb++, itBib++, itBjb++) { 
                  *itQa += ( s1*(*itAjb) - s2*(*itAib) );
                  *itQb += ( s3*(*itBjb) - s4*(*itBib) );
                  *itQab += s5*(*itBjb);
                  *itQab2 += s6*(*itBib);
                }
              }
            }
            // Now I have Q(a,b,k0,l0);
            for(int a=0; a<NAEA; a++)
             for(int b=0; b<NAEA; b++) { 
               if(std::abs(Qa(a,b)) > cut) 
                 Vijkl.add( a*NMO+k0, b*NMO+l0, Qa(a,b), true);
               if(std::abs(Qb(a,b)) > cut) 
                 Vijkl.add( NMO*NMO+a*NMO+k0, NMO*NMO+b*NMO+l0, Qb(a,b), true);
               if(std::abs(Qab(a,b)) > cut)
                 Vijkl.add( a*NMO+k0, NMO*NMO+b*NMO+l0, Qab(a,b), true);
               if(std::abs(Qab2(a,b)) > cut && k0!=l0)
                 Vijkl.add( a*NMO+l0, NMO*NMO+b*NMO+k0, Qab2(a,b), true);
             }

          } else {
            app_error()<<" Error: Finish implementation. \n\n\n";
            return false;
            if(k0 < NMO && l0 < NMO) { // aa 
              for(int a=0; a<NAEA; a++, itAia++, itAja++) {
               s1 = fct*4*V*(*itAia);
               s2 = fct*2*V*(*itAja);
               itQa = Qa.begin(a);
               itAib = Ac.begin(i);
               itAjb = Ac.begin(j);
               for(int b=0; b<NAEA; b++, itQa++, itAib++, itAjb++)
                 *itQa += ( s1*(*itAjb) - s2*(*itAib) );
              }
            } else if(k0 >= NMO && l0 >= NMO) { // bb 
               app_error()<<" Error: Finish implementation. \n\n\n";
               return false;
            } else if(k0 < NMO && l0 >= NMO) { // ab
               app_error()<<" Error: Finish implementation. \n\n\n";
               return false;
            } else {
               app_error()<<" Error: Problems in createHamiltonianForGeneralDeterminant. \n";
               return false;
            }
          } // spinRestricted
        } // KL loop
      } else if(walker_type == 2) {
      }

      if(!communicate_Vijkl(Vijkl)) return false;
      myComm->barrier();

      Timer.stop("Generic");
      app_log()<<"Time to generate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    }  // factorizedHamiltonian

#ifdef AFQMC_DEBUG
    app_log()<<" Done generating sparse hamiltonians. " <<std::endl;
    app_log()<<" Compressing sparse hamiltonians. " <<std::endl;
#endif

    if(head_of_nodes) {

      if(!Vijkl.remove_repeated_and_compress()) {
        APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
      }

#ifdef AFQMC_DEBUG
      app_log()<<" Done compressing sparse hamiltonians. " <<std::endl;
#endif
    }
    myComm->barrier();
    return true;
  }


  bool SparseGeneralHamiltonian::initializeCCProjector(ComplexMatrix& Pmat, RealType cut) {

    ValueType V;
    std::vector<s4D<ValueType> > vs4D;  
    s4D<ValueType> s;

    int nvira = (NMO-NAEA);
    int nvirb = (NMO-NAEB);
    int na = nvira*NAEA;
    int nb = nvirb*NAEB;

    if(NAEA != NAEB) {
      app_error()<<"Error: NAEA != NAEB not supported in CCProjector yet. \n"; 
      return false;
    }
 
    s4Dit it2 = V2.begin();
    while(it2 != V2.end()) {

      if( std::abs( std::get<4>(*it2) ) <= cut ) { it2++; continue; }

      if(spinRestricted) {

        vs4D.clear();
        s = *it2++;
        // generate alpha/alpha sector 
        find_equivalent_OneBar_for_hamiltonian_generation(s,vs4D);

        for(int n=0; n<vs4D.size(); n++) {
          int a=std::get<0>(vs4D[n]);
          int b=std::get<1>(vs4D[n]);
          int i=std::get<2>(vs4D[n]);
          int j=std::get<3>(vs4D[n]);
          if( a >= NAEA && i<NAEA && b >= NAEA && j < NAEA ) { 
            Pmat(i*nvira+(a-NAEA),j*nvirb+(b-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(na+i*nvira+(a-NAEA),j*nvirb+(b-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(i*nvira+(a-NAEA),nb+j*nvirb+(b-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(na+i*nvira+(a-NAEA),nb+j*nvirb+(b-NAEA)) = std::get<4>(vs4D[n]); 

            Pmat(j*nvirb+(b-NAEA),i*nvira+(a-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(j*nvirb+(b-NAEA),na+i*nvira+(a-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(nb+j*nvirb+(b-NAEA),i*nvira+(a-NAEA)) = std::get<4>(vs4D[n]); 
            Pmat(nb+j*nvirb+(b-NAEA),na+i*nvira+(a-NAEA)) = std::get<4>(vs4D[n]); 
          }
        }

      } else {

        app_error()<<"Error: initializeCCProjector not implemented for UHF. DO IT!!! \n";
        return false;

        int sector = getSpinSector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2),std::get<3>(*it2));
        switch(sector) { 
          case 0:
          { 
            break;
          } 
          case 1: 
          { 
            find_equivalent_OneBar_for_hamiltonian_generation(*it2,vs4D);
            for(int n=0; n<vs4D.size(); n++) {
              int a=std::get<0>(vs4D[n]);
              int b=std::get<1>(vs4D[n]);
              int i=std::get<2>(vs4D[n]);
              int j=std::get<3>(vs4D[n]);
              if( a >= NAEA && i<NAEA && b >= NAEA && j < NAEA ) {
                Pmat((a-NAEA)*i,(b-NAEA)*j) = std::get<4>(vs4D[n]);
                Pmat(na+(a-NAEA)*i,(b-NAEA)*j) = std::get<4>(vs4D[n]);
                Pmat((a-NAEA)*i,nb+(b-NAEA)*j) = std::get<4>(vs4D[n]);
                Pmat(na+(a-NAEA)*i,nb+(b-NAEA)*j) = std::get<4>(vs4D[n]);
              }
            }
            break;
          } 
          case 2: 
          { 
            app_error()<<"Error in SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(). I should never get here. \n"; 
            return false;  
            break;
          } 
          case 3: 
          { 
            break;
          } 
        } 

      }
    }

  }

/*
  template< >
  void SQCHamMatElms<double>::calculateEigenvalueDecomposition(int N, std::vector<double>& A, std::vector<double>& eigR, std::vector<double>& eigI, std::vector<double>& V) {   

    int rnk=0;
#if defined(USE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
#endif
    if(rnk!=0) return;

#ifdef _TIMER_
    TIMER.reset(TTAG11);
    TIMER.start(TTAG11);
#endif

    char JOBVL('N');
    char JOBVR('V');
    double dummy;
    int dumI = 1;
    std::vector<double> WORK(1);
    int LWORK = -1;
    int INFO;

    dgeev(JOBVL,JOBVR,N,&(A[0]),N, &(eigR[0]), &(eigI[0]), 
             &dummy,dumI,&(V[0]),N,&(WORK[0]),LWORK,&INFO );

    LWORK = int(WORK[0]);
    std::cout<<"Optimal LWORK used in dgeev: " <<LWORK <<std::endl;
    WORK.resize(LWORK);

    dgeev(JOBVL,JOBVR,N,&(A[0]),N, &(eigR[0]), &(eigI[0]), 
             &dummy,dumI,&(V[0]),N,&(WORK[0]),LWORK,&INFO );

    if(INFO != 0) {
      std::cerr<<"Error in solveEigenvalueProblem: FAIL != 0.\n";
      std::cerr<<"INFO: " <<INFO <<std::endl;
      SQCAbort("Error in solveEigenvalueProblem: FAIL != 0.\n");
    }

#ifdef _TIMER_
    TIMER.stop(TTAG11);
    if(rnk==0) std::cout<<" -- Time to solve eigenvalue problem: " <<TIMER.average(TTAG11) <<"\n";
#endif

  }
*/

}

