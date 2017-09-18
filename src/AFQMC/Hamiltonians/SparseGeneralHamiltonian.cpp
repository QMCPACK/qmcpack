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
#include "AFQMC/Matrix/hdf5_readers.hpp"
#include "AFQMC/Matrix/array_partition.hpp"

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

      if(factorizedHamiltonian ) {
        if(H1.size() != nOne || hij_core.size() != nOne_core) {
          app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl;
          return false;
        }
      } else {
        if(H1.size() != nOne || hij_core.size() != nOne_core  || V2.size()!=nTwo ||  Vijkl_core.size()!=nTwo_core ||  Vijkl_mixed.size()!=nTwo_mixed) {
          app_error()<<"Error after readElementsFromFCIDUMP in readFCIDUMP, wrong number of elements.\n" <<std::endl;
          return false;
        }
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
      app_log()<<" Off-diagonal elements of the Fock matrix larger than " <<5*std::max(cutoff1bar,cutoff1bar) <<std::endl;
      for(IndexType i=0; i<NMO; i++) {
       for(IndexType j=i+1; j<NMO; j++) {
        ValueType fij = H(i,j);
        for(std::vector<IndexType>::iterator it = occup_alpha.begin(); it<occup_alpha.end(); it++)
          fij += H_2bar(i,*it,j,*it);
        for(std::vector<IndexType>::iterator it = occup_beta.begin(); it<occup_beta.end(); it++) 
          fij += H(i,*it,j,*it);
        if( std::abs(fij) > 5*std::max(cutoff1bar,cutoff1bar) )
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

    hdf_write();
    if(rnk == 0) {
      ascii_write();
    }

    MPI_Barrier(TG.getNodeCommLocal()); 
    // generate IJ matrix to speedup table seaches
    if(V2.size()>0)
      generateIJ();

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
         if( std::abs(val)==0 && !spinRestricted ) {
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
         if( std::abs(val)==0 && !spinRestricted ) {
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

    //if(factorizedHamiltonian) 
    //  n_reading_cores=1; // no parallelization yet

    // processor info
    int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
    int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
    int nread = (n_reading_cores<=0)?(ncores):(n_reading_cores);
    int node_rank;
    MPI_Comm_rank(TG.getHeadOfNodesComm(),&node_rank);

    // no hamiltonian distribution yet
    min_i = 0;
    max_i = NMO;
    if(!spinRestricted) max_i *= 2;  // or 4?

    app_log()<<" Initializing Hamiltonian from file: " <<fileName <<std::endl;

    V2.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_V2"),TG.getNodeCommLocal());
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

    hdf_archive dump(myComm);
    // these cores will read from hdf file
    if( coreid < nread ) {

      if(!dump.open(fileName,H5F_ACC_RDONLY)) {
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

    }

    int int_blocks,nvecs;
    std::vector<int> Idata(8);
    std::vector<OrbitalType> ivec;    
    std::vector<ValueType> vvec;    

    if(myComm->rank() == 0) 
      if(!dump.read(Idata,"dims")) {
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading dims. \n"; 
        return false;
      } 
   
    myComm->bcast(Idata);

    int_blocks = Idata[2];
    if(NMO < 0) NMO = Idata[3];
    if(NAEA < 0) NAEA = Idata[4];
    if(NAEB < 0) NAEB = Idata[5];
    if(Idata[3] != NMO) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n"; 
      return false;
    }
    if(Idata[4] != NAEA) {
      app_error()<<" ERROR: NAEA differs from value in integral file. \n"; 
      return false;
    }
    if(Idata[5] != NAEB) {
      app_error()<<" ERROR: NAEA differs from value in integral file. \n"; 
      return false;
    }
    spinRestricted = (Idata[6]==0)?(true):(false);
    factorizedHamiltonian = (Idata[7]>0); 
    nvecs = Idata[7];

    H1.resize(Idata[0]);
    if(myComm->rank() == 0 && distribute_Ham && number_of_TGs > NMO) 
      APP_ABORT("Error: number_of_TGs > NMO. \n\n\n");

    occup_alpha.resize(NAEA);
    occup_beta.resize(NAEB);

    if(myComm->rank() == 0) { 
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
    }
    
    myComm->bcast(occup_alpha);
    myComm->bcast(occup_beta);
    myComm->bcast(NuclearCoulombEnergy);
    myComm->bcast(FrozenCoreEnergy);

    if(myComm->rank() == 0) { 

      ivec.resize(2*H1.size());
      if(!dump.read(ivec,"H1_indx")) {
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading H1_indx. \n";
        return false;
      } 
      for(int i=0, j=0; i<H1.size(); i++, j+=2)        
        H1[i] = std::make_tuple(ivec[j],ivec[j+1],0);  

      vvec.resize(H1.size());
      if(!dump.read(vvec,"H1")) {
        app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading H1.  \n";
        return false;
      }

      for(int i=0; i<H1.size(); i++) {
        // keep i<=j by default
        if(std::get<0>(H1[i]) <= std::get<1>(H1[i]))
            std::get<2>(H1[i]) = vvec[i];
        else {
            std::swap(std::get<0>(H1[i]),std::get<1>(H1[i]));
            std::get<2>(H1[i]) = myconj(vvec[i]);
        }
      }

      std::sort (H1.begin(), H1.end(),mySort);
    } 
    myComm->bcast<char>(reinterpret_cast<char*>(H1.data()),H1.size()*sizeof(s2D<ValueType>));
     
    // now read the integrals
    if(skip_V2) {
      // must anything be done???
    } else if(factorizedHamiltonian) {

      app_log()<<" Reading factorized hamiltonian. \n";  

      min_i = 0;
      max_i = nvecs;

      int nrows = NMO*NMO;
      if(!spinRestricted) nrows *= 2;

      if(rank()==0) {

        afqmc::simple_matrix_partition<afqmc::TaskGroup,IndexType,RealType> split(nrows,nvecs,cutoff1bar);
        std::vector<IndexType> counts;
        // count dimensions of sparse matrix
        afqmc::count_entries_hdf5_SpMat(dump,split,std::string("V2fact"),int_blocks,false,counts,TG,true,nread);

        std::vector<IndexType> sets;
        split.partition_over_TGs(TG,false,counts,sets);

        if(distribute_Ham) {
          app_log()<<" Partitioning of (factorized) Hamiltonian Vectors: ";
          for(int i=0; i<=TG.getNumberOfTGs(); i++)
            app_log()<<sets[i] <<" ";
          app_log()<<std::endl;
          app_log()<<" Number of terms in each partitioning: ";
          for(int i=0; i<TG.getNumberOfTGs(); i++)
            app_log()<<accumulate(counts.begin()+sets[i],counts.begin()+sets[i+1],0) <<" ";
          app_log()<<std::endl;

          MPI_Bcast(sets.data(), sets.size(), MPI_INT, 0, myComm->getMPI());
          min_i = sets[TG.getTGNumber()];
          max_i = sets[TG.getTGNumber()+1];
          //for(int i=0; i<nnodes_per_TG; i++)
          //  nCholVec_per_node[i] = sets[i+1]-sets[i];
          //GlobalSpvnSize = std::accumulate(counts.begin(),counts.end(),0);
        }

        // resize Spvn
        int sz = std::accumulate(counts.begin()+min_i,counts.begin()+max_i,0);
        MPI_Bcast(&sz, 1, MPI_INT, 0, TG.getNodeCommLocal());

        V2_fact.setDims(nrows,nvecs);
        V2_fact.allocate(sz);

        // read Spvn    
        afqmc::read_hdf5_SpMat(V2_fact,split,dump,std::string("V2fact"),int_blocks,TG,true,nread);

      } else {

        afqmc::simple_matrix_partition<afqmc::TaskGroup,IndexType,RealType> split(nrows,nvecs,cutoff1bar);
        std::vector<IndexType> counts;
        // count dimensions of sparse matrix
        afqmc::count_entries_hdf5_SpMat(dump,split,std::string("V2fact"),int_blocks,false,counts,TG,true,nread);

        std::vector<IndexType> sets(TG.getNumberOfTGs()+1);
        if(distribute_Ham) {
          MPI_Bcast(sets.data(), sets.size(), MPI_INT, 0, myComm->getMPI());
          min_i = sets[TG.getTGNumber()];
          max_i = sets[TG.getTGNumber()+1];
        }

        int sz;
        if( coreid==0 )
          sz = std::accumulate(counts.begin()+min_i,counts.begin()+max_i,0);
        MPI_Bcast(&sz, 1, MPI_INT, 0, TG.getNodeCommLocal());

        V2_fact.setDims(nrows,nvecs);
        V2_fact.allocate(sz);

        if(n_reading_cores <=0 || coreid < n_reading_cores) 
          split.partition_over_TGs(TG,false,counts,sets);

        // read Spvn    
        afqmc::read_hdf5_SpMat(V2_fact,split,dump,std::string("V2fact"),int_blocks,TG,true,nread);

      } 
      

/*
      // no parallel reading yet
      if(myComm->rank() == 0) {

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
        if(!dump.read(sz,"V2fact_block_sizes")) {
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
            DenseMatrixOperators::product(NMO,NMO,NMO,ValueType(1),Temp.data(),NMO,rotationMatrix.data(),NMO,ValueType(0),Temp1.data(),NMO); 
            DenseMatrixOperators::product_AhB(NMO,NMO,NMO,ValueType(1.0),rotationMatrix.data(),NMO,Temp1.data(),NMO,ValueType(0.0),Temp.data(),NMO); 
            ivec.clear();
            vvec.clear(); 
            for(int j=0; j<NMO; j++) 
             for(int k=0; k<NMO; k++) {
               if(std::abs(Temp(j,k)) > cutoff1bar) {
                 ivec.push_back(j*NMO+k);
                 vvec.push_back(Temp(j,k));
               } 
             } 
            app_log()<<" Change in Chol Vect. number of terms, old:" <<sz[i] <<"  new:" <<ivec.size() <<std::endl; 
 
          }

          Timer.start("Generic4");
          myComm->bcast(ivec.data(),ivec.size(),MPI_COMM_HEAD_OF_NODES);
          myComm->bcast(vvec.data(),vvec.size(),0,MPI_COMM_HEAD_OF_NODES);
          Timer.stop("Generic4");
 
            for(int k=0; k<ivec.size(); k++) {
              if(std::abs(vvec[k]) > cutoff1bar ) 
                V2_fact.add(ivec[k],i,vvec[k],false);
            }

        }

      } else {  // if(myComm->rank() == 0)

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
            myComm->bcast(vvec.data(),vvec.size(),0,MPI_COMM_HEAD_OF_NODES);

              for(int k=0; k<ivec.size(); k++) {
                if(std::abs(vvec[k]) > cutoff1bar )
                  V2_fact.add(ivec[k],i,vvec[k],false);
              }

          }

        }

      }
      Timer.reset("Generic2");
      Timer.start("Generic2");
      V2_fact.compress(TG.getNodeCommLocal());
      Timer.stop("Generic2");

      app_log()<<" -- Average time to read std::vector from h5 file: " <<Timer.average("Generic3") <<"\n";
      app_log()<<" -- Average time to bcast std::vector: " <<Timer.average("Generic4") <<"\n";
      app_log()<<" -- Time to compress factorized Hamiltonian from h5 file: " <<Timer.average("Generic2") <<"\n";
*/
      app_log()<<" Memory used by factorized 2-el integral table: " <<V2_fact.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;
 
    } else {  // factorizedHamiltonian

      std::vector<long> ntpo(NMO,0);
      std::vector<IndexType> indxvec;
      int ntmax;
      std::vector<int> pool_dist;
      if(coreid < nread) {

        Idata.resize(int_blocks);
        // Idata[i]: number of terms per block
        if(!dump.read(Idata,"V2_block_sizes")) {
          app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_block_sizes dataset. \n";
          return false;
        }

        ntmax = *std::max_element(Idata.begin(), Idata.end());
        indxvec.reserve(2*ntmax);
        vvec.reserve(ntmax);

        // divide blocks 
        FairDivide(int_blocks,nread,pool_dist);

        int first_local_block = pool_dist[coreid];
        int last_local_block = pool_dist[coreid+1];

        for(int ib=first_local_block; ib<last_local_block; ib++) {
          if( ib%nnodes != nodeid ) continue;
          if(Idata[ib]==0) continue;
          indxvec.resize(2*Idata[ib]);
          vvec.resize(Idata[ib]);
          if(!dump.read(indxvec,std::string("V2_index_")+std::to_string(ib))) {
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_index_" <<ib <<" dataset. \n";
            return false;
          }
          if(!dump.read(vvec,std::string("V2_vals_")+std::to_string(ib))) {
            app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_vals_" <<ib <<" dataset. \n";
            return false;
          }
            
          std::vector<IndexType>::iterator iti = indxvec.begin();
          for(std::vector<ValueType>::iterator itv = vvec.begin(); itv < vvec.end(); itv++, iti+=2) {
            if(std::abs(*itv) < cutoff1bar )
                continue;
            s4D<ValueType> ijkl = std::make_tuple(  static_cast<OrbitalType>((*iti)/NMO),
                                                    static_cast<OrbitalType>((*(iti+1))/NMO),  
                                                    static_cast<OrbitalType>((*iti)%NMO),
                                                    static_cast<OrbitalType>((*(iti+1))%NMO), ValueType(0));     
            find_smallest_permutation(ijkl); 
            ntpo[std::get<0>(ijkl)]++;
          }   
        }
      }
      myComm->allreduce(ntpo);

      if(distribute_Ham) {
        std::vector<long> sets(TG.getNumberOfTGs()+1);
        if(myComm->rank() == 0) {
          std::vector<long> nv(NMO+1);
          nv[0]=0;
          for(int i=0; i<NMO; i++) 
            nv[i+1]=nv[i]+ntpo[i];
          balance_partition_ordered_set(NMO,nv.data(),sets);
          app_log()<<" Hamiltonian partitioning: \n   Orbitals:        ";
          for(int i=0; i<=TG.getNumberOfTGs(); i++) app_log()<<sets[i] <<" "; 
          app_log()<<std::endl <<"   Terms per block: ";  
          for(int i=0; i<TG.getNumberOfTGs(); i++) app_log()<<nv[sets[i+1]]-nv[sets[i]] <<" ";
          app_log()<<std::endl;  
        }
        MPI_Bcast(sets.data(),sets.size(),MPI_LONG,0,myComm->getMPI());
        min_i = static_cast<int>(sets[TG.getTGNumber()]);
        max_i = static_cast<int>(sets[TG.getTGNumber()+1]);
      }

      long nttot = std::accumulate(ntpo.begin()+min_i,ntpo.begin()+max_i,long(0));
      V2.reserve(nttot);

      if( coreid < nread ) {

        std::vector<IndexType> ivec2;
        std::vector<ValueType> vvec2;
        ivec2.reserve(2*ntmax);
        vvec2.reserve(ntmax);
        int maxv = 10000;
        std::vector<s4D<ValueType>> vals;
        vals.reserve(maxv);

        int first_local_block = pool_dist[coreid];
        int last_local_block = pool_dist[coreid+1];
        int nbtot = last_local_block-first_local_block;
        int niter = nbtot/nnodes + std::min(nbtot%nnodes,1);

        for(int iter=0; iter<niter; iter++) {
          int first_block = first_local_block + nnodes*iter;
          int last_block = std::min(first_block+nnodes,last_local_block);
          int myblock_number = first_block + nodeid;
          if(myblock_number < last_local_block && Idata[myblock_number] > 0) {
            indxvec.resize(2*Idata[myblock_number]);
            vvec.resize(Idata[myblock_number]);
            if(!dump.read(indxvec,std::string("V2_index_")+std::to_string(myblock_number))) {
              app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_index_" <<myblock_number <<" dataset. \n";
              return false;
            }
            if(!dump.read(vvec,std::string("V2_vals_")+std::to_string(myblock_number))) {
              app_error()<<" Error in SparseGeneralHamiltonian::initFromHDF5(): Problems reading V2_vals_" <<myblock_number <<" dataset. \n";
              return false;
            }
          }
         
          for(int k=first_block,ipr=0; k<last_block; k++,ipr++) {
            if(Idata[k] == 0) continue;
            ivec2.resize(2*Idata[k]);
            vvec2.resize(Idata[k]);
            if(ipr==node_rank) {
              assert(myblock_number==k);
              std::copy(indxvec.begin(),indxvec.end(),ivec2.begin());
              std::copy(vvec.begin(),vvec.end(),vvec2.begin());
            }

            MPI_Bcast(ivec2.data(), ivec2.size(), MPI_INT, ipr, TG.getHeadOfNodesComm() );
#if defined(QMC_COMPLEX)
            MPI_Bcast(vvec2.data(),2*vvec2.size(), MPI_DOUBLE, ipr, TG.getHeadOfNodesComm() );
#else
            MPI_Bcast(vvec2.data(),vvec2.size(), MPI_DOUBLE, ipr, TG.getHeadOfNodesComm() );
#endif

            std::vector<IndexType>::iterator iti = ivec2.begin(); 
            for(std::vector<ValueType>::iterator itv = vvec2.begin(); itv < vvec2.end(); itv++, iti+=2) { 
              if(std::abs(*itv) < cutoff1bar )
                  continue;
              s4D<ValueType> ijkl = std::make_tuple(  static_cast<OrbitalType>((*iti)/NMO),
                                                      static_cast<OrbitalType>((*(iti+1))/NMO),  
                                                      static_cast<OrbitalType>((*iti)%NMO),
                                                      static_cast<OrbitalType>((*(iti+1))%NMO), *itv);     
              find_smallest_permutation(ijkl); 
              if( std::get<0>(ijkl) >= min_i && std::get<0>(ijkl) < max_i) { 
                vals.push_back(ijkl);
                if(vals.size()==maxv) {
                  V2.push_back(vals,nread>1);
                  vals.clear();
                }
              }
            }
          }   
        }
        if(vals.size() > 0) 
          V2.push_back(vals,nread>1);
      }

      Timer.reset("Generic2");
      Timer.start("Generic2");
      V2.sort (mySort, TG.getNodeCommLocal(),inplace);
      Timer.stop("Generic2");
      app_log()<<" -- Time to compress Hamiltonian from h5 file: " <<Timer.average("Generic2") <<"\n";
      app_log()<<" Memory used by 2-el integral table: " <<V2.memoryUsage()/1024.0/1024.0 <<" MB. " <<std::endl;

    }

    if( coreid < nread ) {
      dump.pop();
      dump.pop();
      dump.close();
    }
    myComm->barrier();

    // generate IJ matrix to speedup table seaches
    if(!skip_V2 && V2.size()>0)
      generateIJ();

    Timer.stop("Generic");
    app_log()<<" -- Time to initialize Hamiltonian from h5 file: " <<Timer.average("Generic") <<"\n";

    if(!skip_V2) hdf_write();
    if(rank() == 0 && !skip_V2) {
      ascii_write();
    }

/*   Meant for debugging integral codes, remove for released code
    if(rank() == 0 ) {
      std::vector<ValueType> eig(NMO);
      for(IndexType i=0; i<NMO; i++) {
        eig[i] = H(i,i);
        for(std::vector<IndexType>::iterator it = occup_alpha.begin(); it<occup_alpha.end(); it++)
          eig[i] += H_2bar(i,*it,i,*it);
        for(std::vector<IndexType>::iterator it = occup_beta.begin(); it<occup_beta.end(); it++)
          eig[i] += H(i,*it,i,*it);
        app_log()<<i <<"   " <<eig[i] <<std::endl;
      }
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
      ValueType Emp2=0;  
      std::vector<IndexType>::iterator iti,itj,ita,itb;
      for(iti=occup_alpha.begin(); iti<occup_alpha.end(); iti++)
      for(itj=occup_alpha.begin(); itj<occup_alpha.end(); itj++)
      for(ita=virtual_alpha.begin(); ita<virtual_alpha.end(); ita++)
      for(itb=virtual_alpha.begin(); itb<virtual_alpha.end(); itb++)
        Emp2 += H(*iti,*itj,*ita,*itb)*(2.0*H(*ita,*itb,*iti,*itj) - H(*ita,*itb,*itj,*iti))/(eig[*iti]+eig[*itj]-eig[*ita]-eig[*itb]);
      app_log()<<" Emp2: " <<Emp2 <<std::endl;
    }
*/
    return true;

  }

  void SparseGeneralHamiltonian::ascii_write() {

    if(ascii_write_file == std::string("")) return;

    if(factorizedHamiltonian) {
      app_log()<<" factorizedHamiltonian output not yet implemented " <<std::endl;
      return;
    }
 
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
   for(ComplexSMSpMat::intType i=0; i<nvecs; i++) {
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
      return;
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

    if(skip_V2) 
      return; 

    if(hdf_write_file == std::string("")) return;

    bool writeFactorized = factorizedHamiltonian;
    if(hdf_write_type!="default") {
      if(factorizedHamiltonian) writeFactorized = (hdf_write_type == "factorized");
      else writeFactorized = (hdf_write_type != "integrals");
    }

    // pseudo-hack to avoid a lot of indenting 
    if(myComm->rank() > 0) {


      // only parallelized path
      if(!writeFactorized && factorizedHamiltonian) {
        std::vector<OrbitalType> ivec; 
        std::vector<ValueType> vvec;
        for(int k=0; k<NMO; k++) {
          ivec.clear();
          vvec.clear();
          SparseHamiltonianFromFactorization(k,ivec,vvec,cutoff1bar);
        }
      }
  
      return; 
    }

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

    int V2tot=0;
    std::vector<int> Idata(8);
/* writing at the end since I don't have all the information necessarily
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
*/

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

    if(writeFactorized) {

      if(!factorizedHamiltonian) {
        app_error()<<" Error: Can only write factorized hamiltonian if input hamiltonian is already factorized. \n\n\n";
        APP_ABORT(" Error: Can only write factorized hamiltonian if input hamiltonian is already factorized. \n\n\n");
      }

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
      std::vector<ValueType> vvec;
      ivec.reserve(nmax);
      vvec.reserve(nmax);
      for(ComplexSMSpMat::intType i=0; i<nvecs; i++) {

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

      if(factorizedHamiltonian) {

        for(int k=0; k<NMO; k++) {
          ivec.clear();
          vvec.clear();
          SparseHamiltonianFromFactorization(k,ivec,vvec,cutoff1bar); 
          Idata[k]= (k==0?0:Idata[k-1]+vvec.size());
          V2tot+=vvec.size();
          dump.write(ivec,std::string("V2_indx_")+std::to_string(k));
          dump.write(vvec,std::string("V2_")+std::to_string(k));
        }

        dump.write(Idata,"V2_block_sizes");

      } else {
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

        dump.write(Idata,"V2_block_sizes");

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
          V2tot+=vvec.size();
          dump.write(vvec,std::string("V2_")+std::to_string(k));
        }
 
      }

    }

    Idata.resize(8);
    Idata[0]=H1.size();
    Idata[1]=V2tot;
    Idata[2]=0;
    Idata[3]=NMO;
    Idata[4]=NAEA;
    Idata[5]=NAEB;
    Idata[6]=spinRestricted?(0):(1);
    Idata[7]=0;
    if(writeFactorized && factorizedHamiltonian)
      Idata[7] = V2_fact.cols();
    dump.write(Idata,"dims");

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
    //    Now also includes contribution from vn0 = -0.5* sum_{i,l,sigma} (sum_j <i_sigma,j_sigma|j_sigma,l_sigma> ) c+i_sigma cl_sigma
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


      // adding contribution from transformation of H2 into quadratic form
      // 3.  -0.5* sum_{i,l,sigma} (sum_j <i_sigma,j_sigma|j_sigma,l_sigma> ) c+i_sigma cl_sigma 
      //   Calculated with HSPotential and included in Hadd!!!

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
       DenseMatrixOperators::product(NMO,NAEA,NMO,ComplexType(1),P.data(),NMO,S.data(),NAEA,ComplexType(0),Snew.data(),NAEA); 

     Timer.stop("Generic2");
     app_log()<<" -- Average time for dense x dense MM product: " <<Timer.average("Generic2")/nterms <<"\n";

     Timer.reset("Generic2");
     Timer.start("Generic2");

     for(int i=0; i<nterms; i++)
       SparseMatrixOperators::product_SD(NAEA,Pkin.data(),ntermsalpha,S.data(),NAEA,Snew.data(),NAEA);

     Timer.stop("Generic2");
     app_log()<<" -- Average time for sparse x dense MM product: " <<Timer.average("Generic2")/nterms <<"\n";


  }

  void SparseGeneralHamiltonian::calculateHSPotentials_Diagonalization(RealType cut, RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node,  bool sparse, bool parallel)
  {

    if(skip_V2)
      APP_ABORT("Error: Calling SparseGeneralHamiltonian routines with skip_V2=yes. \n\n\n");

    if(factorizedHamiltonian) {
      calculateHSPotentials_FactorizedHam(cut,dt,vn0,Spvn,Dvn,TGprop,nvec_per_node,sparse,parallel);
      return;
    }

    if(parallel) {
      APP_ABORT("Error: calculateHSPotentials_Diagonalization not implemented in parallel. \n");
    }

    // 3. vn0 = -0.5* sum_{i,l,sigma} (sum_j <i_sigma,j_sigma|j_sigma,l_sigma> ) c+i_sigma cl_sigma 
    if(myComm->rank()==0) 
    {  
      if(distribute_Ham) {
        APP_ABORT("Error: calculateHSPotentials_Diagonalization not implemented with distributed hamiltonian. \n\n\n"); 
      }  
      if(spinRestricted) 
        assert(vn0.rows() >= NMO);
      else
        assert(vn0.rows() >= 2*NMO);
      assert(vn0.cols() == NMO);
      std::vector<s4D<ValueType> > v2sym; 
      vn0 = ValueType(0);
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
              vn0(i,Index2Col(l)) -= 0.5*V;
        }
      }  
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
        
      // for dense storage, setting cut to 1e-8 mainly to avoid storing empty vectors
      // with real integrals
      if(!sparse) cut = 1e-8; 
      else if(cut < 1e-12) cut=1e-12;
      int cnt1=0;
      int cnt2=0;
      int cntn=0;
      for(int i=0; i<NMO2; i++) { 
       if(std::abs(eigVal[i]) > std::abs(cutoff_cholesky)) { 
#ifndef QMC_COMPLEX
         if(eigVal[i] < 0) {
           app_log()<<" WARNING: Found negative eigenvalue in REAL build. Ignoring it: " <<i <<" " <<eigVal[i] <<"\n";
           continue;
         }
#endif
         ValueType eig = std::sqrt( 0.25*dt*eigVal[i] );
         int cnt3=0;
         for(int j=0; j<NMO; j++) 
           for(int k=0; k<NMO; k++) { 
             IndexType jk = j*NMO+k; 
             IndexType kj = k*NMO+j; 
             ValueType V = eig*(eigVec(jk,i) + myconj(eigVec(kj,i)));
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
#if defined(QMC_COMPLEX)
         cnt3=0;
         for(int j=0; j<NMO; j++) 
          for(int k=0; k<NMO; k++) { 
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ValueType V = eig*(eigVec(jk,i) - myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) cnt3++;
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) - myconj(eigVec(NMO+NMO+kj,i)));
              if(std::abs(V) > cut) cnt3++;
            }
          }
         if(cnt3 > 0) {
           cntn++;
           cnt1++;
           cnt2 += cnt3;
         } 
#endif
        }
      } 

      if(sparse) {
        Spvn.setDims(NMO2,cnt1);
        Spvn.allocate_serial(cnt2);
      } else {
        Dvn.allocate_serial(NMO2*cnt1);
      }

#if defined(QMC_COMPLEX)
      if(cntn>0) 
        APP_ABORT("Found real Cholesky vectors with real integrals. This is not allowed with dense cholesky vectors. Run with sparse vectors or compile with complex integrals. \n\n\n");
#endif 
 
      int ncols=cnt1;
      cnt1=0;
      cntn=0;
      for(int i=0; i<NMO2; i++) {
       if(std::abs(eigVal[i]) > std::abs(cutoff_cholesky)) {
#ifndef QMC_COMPLEX
         if(eigVal[i] < 0) continue; 
#endif
         ValueType eig = std::sqrt( 0.25*dt*eigVal[i] );
         int cnt3=0;
         for(int j=0; j<NMO; j++)
          for(int k=0; k<NMO; k++) {
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ValueType V = eig*(eigVec(jk,i) + myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) {
              cnt3++;
              //V*=isqrtdt;
              if(sparse) Spvn.add(jk,cnt1,static_cast<SPValueType>(V));
              else Dvn[jk*ncols+cnt1]=static_cast<SPValueType>(V);
            }
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) + myconj(eigVec(NMO*NMO+kj,i)));
              if(std::abs(V) > cut) {
                cnt3++;
                //V*=isqrtdt;
                if(sparse) Spvn.add(NMO*NMO+jk,cnt1,static_cast<SPValueType>(V));
                else Dvn[(jk+NMO*NMO)*ncols+cnt1]=static_cast<SPValueType>(V);
              }
            }
          }
         if(cnt3 > 0) {
           cnt1++; 
         }
#if defined(QMC_COMPLEX)
         eig = ComplexType(0.0,1.0)*std::sqrt( 0.25*dt*eigVal[i] );
         cnt3=0;
         for(int j=0; j<NMO; j++)
          for(int k=0; k<NMO; k++) {
            IndexType jk = j*NMO+k;
            IndexType kj = k*NMO+j;
            ValueType V = eig*(eigVec(jk,i) - myconj(eigVec(kj,i)));
            if(std::abs(V) > cut) {
              cnt3++;
              //V*=sqrtdt;
              if(sparse) Spvn.add(jk,cnt1,static_cast<SPValueType>(V));
              else Dvn[jk*ncols+cnt1]=static_cast<SPValueType>(V); 
            }
            if(!spinRestricted) {
              V = eig*(eigVec(NMO*NMO+jk,i) - myconj(eigVec(NMO*NMO+kj,i)));
              if(std::abs(V) > cut) {
                cnt3++;
                //V*=sqrtdt;
                if(sparse) Spvn.add(NMO*NMO+jk,cnt1,static_cast<SPValueType>(V));
                else Dvn[(jk+NMO*NMO)*ncols+cnt1]=static_cast<SPValueType>(V); 
              }
            }
          }
         if(cnt3 > 0) {
           cnt1++; 
           cntn++;
         }
#endif
       }
      }

      app_log()<<"Number of positive and negative Cholesky std::vectors: " <<cnt1-cntn <<" " <<cntn <<std::endl; 

      app_log()<<"Number of HS potentials: " <<ncols <<std::endl;
      if(sparse) {
        app_log()<<"Number of terms in sparse representation of HS potentials: " <<Spvn.size() <<std::endl;
        app_log()<<"Compressing Spvn. \n";
        Spvn.compress();
        app_log()<<"Done Compressing Spvn. \n";
     }

     if(test_breakup && sparse) {

      if(rnk==0) app_log()<<" -- Testing Hamiltonian factorization. \n";
      Timer.reset("Generic");
      Timer.start("Generic");

      int* cols = Spvn.column_data();
      int* rows = Spvn.row_data();
      int* indx = Spvn.row_index();
      SPValueType* vals = Spvn.values(); 

      SPComplexMatrix v2prod(NMO2);
      for(int i=0; i<NMO2; i++) 
       for(int j=0; j<NMO2; j++) 
        v2prod(i,j) = SparseMatrixOperators::product_SpVSpV<SPValueType>(indx[i+1]-indx[i],cols+indx[i],vals+indx[i],indx[j+1]-indx[j],cols+indx[j],vals+indx[j]);   

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
           ComplexType v2 = dt*H(i,j,k,l);
           ComplexType v2c = static_cast<ComplexType>(v2prod( i*NMO+Index2Col(k) , j*NMO+Index2Col(l)));
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

  void SparseGeneralHamiltonian::calculateHSPotentials_FactorizedHam(RealType cut, RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool parallel)
  {

    if(skip_V2)
      APP_ABORT("Error: Calling SparseGeneralHamiltonian routines with skip_V2=yes. \n\n\n");

    if(distribute_Ham) {
      APP_ABORT("Error: calculateHSPotentials_FactorizedHam not implemented with distributed Hamiltonian. \n");
    }

    assert(vn0.rows() >= NMO);
    assert(vn0.cols() == NMO);
    vn0 = ValueType(0);
    for(int i=0; i<NMO; i++)
     for(int l=i; l<NMO; l++) {
       ValueType vl = ValueType(0);
       for(int j=0; j<NMO; j++)
         vl += H(i,j,j,l);
       vn0(i,l) -= 0.5*vl;
       if(i!=l) vn0(l,i) -= 0.5*myconj(vl);
       if(!spinRestricted) {
         vl=ValueType(0);
         for(int j=NMO; j<2*NMO; j++)
           vl += H(i+NMO,j,j,l+NMO);
         vn0(i+NMO,Index2Col(l+NMO)) -= 0.5*vl;
         if(i!=l) vn0(l+NMO,Index2Col(l+NMO)) -= 0.5*myconj(vl);
       }
     }

    /********************************************************************
    *  You get 2 potentials per Cholesky std::vector   
    *
    *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
    ********************************************************************/

    ValueType sqrtdt = std::sqrt(dt)*0.5;
    int NMO2 = NMO*NMO;
    if(!spinRestricted) NMO2*=2; 

    int n3Vect = V2_fact.cols();
    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

    int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
    if(!parallel) {
      ncores=1;
      coreid=0;
    }

    // transpose V2_fact here and again at the end
    Timer.reset("Generic");
    Timer.start("Generic");
    if(parallel) V2_fact.transpose(TGprop.getNodeCommLocal());
    else if(head_of_nodes) V2_fact.transpose(); 
    Timer.stop("Generic");
    app_log()<<" Time to transpose V2_fact inside calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

    int* cols = V2_fact.column_data();
    int* rows = V2_fact.row_data();
    int* indx = V2_fact.row_index();
    ValueType* vals = V2_fact.values();
    ValueMatrix Ln;

    Ln.resize((spinRestricted?(NMO):(2*NMO)) ,NMO);  

    std::vector<int> cnt_per_vec(2*n3Vect);
    int nnodes = TGprop.getNNodesPerTG();
    int cv0=0,cvN=2*n3Vect;
    std::vector<int> sets;
    std::vector<int> sz_per_node(nnodes);
    int node_number = TGprop.getLocalNodeNumber();
    int rk=0,npr=1;
    if(parallel) {
      npr = myComm->size();
      rk = myComm->rank();
    }

    if(!sparse) cut=1e-8;
    else if(cut < 1e-12) cut=1e-12;
    int cnt=0;
    int cntn=0;
    // count terms and distribute std::vectors 
    // generate sparse version
    for(int n=0; n<n3Vect; n++) { 
      if(n%npr != rk) continue;
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
      if(nm>0) cntn++;
    }

    if(parallel) {
      myComm->allreduce(cnt_per_vec);
      int nvec = std::count_if(cnt_per_vec.begin(),cnt_per_vec.end(),
               [] (int i) { return i>0; } );
      if(sparse) Spvn.setDims(NMO2,nvec); 

      nvec_per_node.resize(nnodes);
      if(nnodes==1) {
        cv0=0;
        cvN=2*n3Vect;
        cnt = std::accumulate(cnt_per_vec.begin(),cnt_per_vec.end(),0);
        if(sparse) Spvn.reserve(cnt);
        else {
          Dvn.allocate(NMO2*nvec);
          Dvn.resize(NMO2*nvec);
        }
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
            if(sparse) cnt+=cnt_per_vec[i];
            else cnt+= (cnt_per_vec[i]>0)?1:0;
            blocks[i+1] = cnt;
          }
          balance_partition_ordered_set(2*n3Vect,blocks.data(),sets);
          myComm->bcast(sets);
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
          myComm->bcast(sets);
          cv0 = sets[node_number];
          cvN = sets[node_number+1];
        } 
        myComm->bcast(nvec_per_node);
        myComm->bcast(sz_per_node);
        if(sparse) Spvn.reserve(sz_per_node[node_number]);
        else {
          Dvn.allocate(NMO2*nvec_per_node[node_number]);
          Dvn.resize(NMO2*nvec_per_node[node_number]);
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
        Spvn.setDims(NMO2,nvec);
        Spvn.allocate_serial(cnt);
      } else {
        Dvn.allocate_serial(NMO2*nvec);
        Dvn.resize_serial(NMO2*nvec);
      }
    }            

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 
    int ncols = nvec_per_node[node_number];

#ifndef QMC_COMPLEX
    if(cntn>0)
      APP_ABORT("Found real Cholesky vectors with real integrals. This is not allowed with dense cholesky vectors. Run with sparse vectors or compile with complex integrals. \n\n\n");
#endif

    int nmax=1000000;
    std::vector<std::tuple<SPValueSMSpMat::intType,SPValueSMSpMat::intType,SPValueType>> vikn;
    if(sparse)
      vikn.reserve(nmax);

// also no need to do this serially, FIX FIX FIX
    cnt=std::accumulate(nvec_per_node.begin(),nvec_per_node.begin()+node_number,0);
    for(int n=0; n<n3Vect; n++) { 
       if( cnt_per_vec[2*n]==0 && cnt_per_vec[2*n+1]==0 ) continue;
       if( !(2*n>=cv0 && 2*n<cvN) && !((2*n+1)>=cv0 && (2*n+1)<cvN)  ) continue;
       if( n%ncores != coreid ) { 
         if(2*n>=cv0 && 2*n<cvN && cnt_per_vec[2*n]>0) cnt++;
         if((2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) cnt++;
         continue;  
       }
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
           ValueType V = (Ln(i,k) + myconj(Ln(k,i))); 
           if(std::abs(V) > cut) { 
             V*=sqrtdt;
             if(sparse) { //Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(V));
               vikn.push_back(std::make_tuple(i*NMO+k,cnt,static_cast<SPValueType>(V)));
               if(vikn.size()==nmax) {  
                 Spvn.add(vikn,true);
                 vikn.clear(); 
               }  
             } else Dvn[(i*NMO+k)*ncols+cnt]=static_cast<SPValueType>(V);
             ++np;
           } 
           if(!spinRestricted) {
             V = (Ln(i+NMO,k) + myconj(Ln(k+NMO,i)));
             if(std::abs(V) > cut) {
               V*=sqrtdt;
               if(sparse) { // Spvn.add(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(V));
                 vikn.push_back(std::make_tuple(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(V)));
                 if(vikn.size()==nmax) {  
                   Spvn.add(vikn,true);
                   vikn.clear(); 
                 }  
               } else Dvn[(i*NMO+k+NMO*NMO)*ncols+cnt]=static_cast<SPValueType>(V);
               ++np;
             }
           }
         }
         ++cnt;
         if(np==0)
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n");
       }
       if((2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) {
#if defined(QMC_COMPLEX) 
         np=0;
         // v-
         for(IndexType i=0; i<NMO; i++)
          for(IndexType k=0; k<NMO; k++) {
           ComplexType V = (Ln(i,k) - myconj(Ln(k,i))); 
           if(std::abs(V) > cut) {
             V*=ComplexType(0.0,1.0)*sqrtdt;
             if(sparse) { // Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(V));
               vikn.push_back(std::make_tuple(i*NMO+k,cnt,static_cast<SPValueType>(V)));
               if(vikn.size()==nmax) {  
                 Spvn.add(vikn,true);
                 vikn.clear(); 
               }  
             } else Dvn[(i*NMO+k)*ncols+cnt]=static_cast<SPValueType>(V);
             ++np;
           }
           if(!spinRestricted) {
             V = (Ln(i+NMO,k) - myconj(Ln(k+NMO,i)));
             if(std::abs(V) > cut) {
               V*=ComplexType(0.0,1.0)*sqrtdt;
               if(sparse) { //Spvn.add(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(V));
                 vikn.push_back(std::make_tuple(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(V)));
                 if(vikn.size()==nmax) {
                   Spvn.add(vikn,true);
                   vikn.clear();
                 }
               } else Dvn[(i*NMO+k+NMO*NMO)*ncols+cnt]=static_cast<SPValueType>(V);
               ++np;
             }
           }
         }
         ++cnt;
         if(np==0)
           APP_ABORT("Error: This should not happen. Found empty cholesky std::vector. \n");
#else
         APP_ABORT("Error: THIS SHOULD NOT HAPPEN. Found complex cholesky vector in real build. \n\n\n");
#endif
       }
    }
    if(vikn.size()>0) {  
      Spvn.add(vikn,true);
      vikn.clear(); 
    }  

    if(rank()==0 && nnodes>1 && parallel) {
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

      if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

      Timer.reset("Generic");
      Timer.start("Generic");
      if(parallel) Spvn.compress(TGprop.getNodeCommLocal());
      else if(head_of_nodes) Spvn.compress();
      Timer.stop("Generic");
      app_log()<<" Time to compress Spvn in calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;
    }

    Timer.reset("Generic");
    Timer.start("Generic");
    if(parallel) V2_fact.transpose(TGprop.getNodeCommLocal());
    else if(head_of_nodes) V2_fact.transpose(); 
    Timer.stop("Generic");
    app_log()<<" Time to transpose V2_fact inside calculateHSPotentials_FactorizedHam(): " <<Timer.average("Generic") <<std::endl;

    if(parallel) MPI_Barrier(TGprop.getNodeCommLocal()); 

  }

  void SparseGeneralHamiltonian::calculateHSPotentials(RealType cut, RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool parallel)
  {

    if(skip_V2)
      APP_ABORT("Error: Calling SparseGeneralHamiltonian routines with skip_V2=yes. \n\n\n");

    if(factorizedHamiltonian) {
      calculateHSPotentials_FactorizedHam(cut,dt,vn0,Spvn,Dvn,TGprop,nvec_per_node,sparse,parallel);
      return;
    }

    if(distribute_Ham && !parallel) {
        APP_ABORT("Error: Distributed hamiltonian requires parallel Cholesky factorization. \n\n\n");
    }

    int rnk=0;
    rnk = rank();
    cholesky_residuals.clear();
    cholesky_residuals.reserve(2*NMO*NMO);

    if(spinRestricted) {


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
          npr = static_cast<long>(myComm->size());
          rk = static_cast<long>(myComm->rank());
        }
        for(s4Dit it = V2.begin(); it != V2.end(); it++, nt++) {
          if( nt%npr != rk) continue;
          // dumb and slow, but easy for now
          find_equivalent_OneBar_for_integral_list(*it,v2sym); 
          for(int n=0; n<v2sym.size(); n++) {   
            IndexType i,j,k,l;
            ValueType V;
            std::tie (i,j,k,l,V) = v2sym[n];   
            int sector = getSpinSector(i,j,k,l);
            // <i,j | j,l> -> v(i,l)  
            if( (sector==0 || sector==3) && j==k) 
              vn0(i,Index2Col(l)) -= 0.5*V;
          }
        }  
        if(parallel) {
          // since allreduce(vn0) somehow doesn't work
          //myComm->allreduce(vn0);
          std::vector<ComplexType> g(NMO*NMO);
          std::copy(vn0.begin(),vn0.begin()+NMO*NMO,g.begin());
          myComm->allreduce(g);
          std::copy(g.begin(),g.end(),vn0.begin());
        }
      }
             

      /********************************************************************
      *               Calculate Cholesky decomposition 
      *
      *   1. The mapping of the 2-el repulsion integrals to a 2-D matrix
      *      is done as follows:
      *         V(i,j,k,l) -->  V( i*NMO+k,  l*NMO+j )
      ********************************************************************/

     Timer.reset("Generic");
     Timer.start("Generic");

     // used to split (i,k) space over all processors for the construction and storage of cholesky vectors
     int npr = myComm->size(), rk = myComm->rank(); 
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
     std::vector<s2D<RealType> > IKLmax(npr); 
     s2D<RealType> mymax;
     int maxloc=0;

     std::vector<ValueType> Lcomm;
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
             app_log()<<" WARNING: Found negative/complex Duv: " <<i <<" " <<k <<" " <<Lcomm[nt] <<std::endl;
             if(zero_bad_diag_2eints)
               Lcomm[nt] = ValueType(0);
           }
#else
           if(Lcomm[nt] < ValueType(0)) {
             app_log()<<" WARNING: Found negative Duv: " <<i <<" " <<k <<" " <<Lcomm[nt] <<std::endl;
             if(zero_bad_diag_2eints)
               Lcomm[nt] = ValueType(0);
           }
#endif
         }
         if(nt>=tg_ik1) break;
       }
       myComm->allreduce(Lcomm); 
       std::copy( Lcomm.begin()+ik0, Lcomm.begin()+ik1, Duv.begin() );
     } else {
       for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
         for(IndexType k=0; k<NMO; k++,nt++) {
           if(nt<ik0 || nt>=ik1) continue;
           Duv[ik] = H(i,k,k,i);  
#ifndef QMC_COMPLEX
           if(Duv[ik] < ValueType(0)) {
             app_log()<<" WARNING: Found negative Duv: " <<i <<" " <<k <<" " <<Duv[ik] <<std::endl;  
             if(zero_bad_diag_2eints) 
               Duv[ik] = ValueType(0);
           }
#else
           if(Duv[ik].imag() > 1e-10 || Duv[ik].real() < RealType(0)) {
             app_log()<<" WARNING: Found negative/complex Duv: " <<i <<" " <<k <<" " <<Duv[ik] <<std::endl;
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
     mymax = std::make_tuple(-1,-1,0);
     for(IndexType i=0, nt=0, ik=0; i<NMO; i++) {
      for(IndexType k=0; k<NMO; k++,nt++) {
        if(nt<ik0 || nt>=ik1) continue;
        if( std::abs(Duv[ik]) > max) {
          max = std::get<2>(mymax) =std::abs(Duv[ik]);  
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

     if(test_2eint && !parallel && !distribute_Ham) {

       /* <ij|kl> <--> M_(ik)_(lj) where M must be positive definite.
        * --> |M_(ik)_(lj)| <= sqrt( |M_(ik)_(ik)| |M_(lj)_(lj)| )
        * --> |<ij|kl>| <= sqrt( |<ik|ki>| |<lj|jl| )   
        *  
        */
       for(s4Dit it = V2.begin(); it != V2.end(); it++) {
         IndexType i,j,k,l;
         ValueType w1,w2,w3;
         std::tie (i,j,k,l,w1) = *it;  
 
         w2 = H(i,k,k,i);
         w3 = H(l,j,j,l);
         ValueType w4 = H(j,l,l,j);
         if( std::abs(w1) > std::sqrt(std::abs(w2*w3)) ) {
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

       cholesky_residuals.push_back(std::abs(max));

       // calculate new cholesky std::vector based on (ii,kk)
       L.push_back(std::vector<ValueType>(maxnterms));  
       std::vector<ValueType>& Ln = L.back();
       std::vector<ValueType>::iterator it = Ln.begin();

       if(rk==maxloc) {
         for(int n=0; n<L.size()-1; n++)
           Lnmax[n] = L[n][ii*NMO+kk-ik0]; 
       }
       if(parallel && L.size()>1) {
         myComm->bcast(Lnmax.data(),L.size()-1,maxloc,myComm->getMPI()); 
       }

       if(distribute_Ham) {
         std::fill(Lcomm.begin(),Lcomm.end(),ValueType(0));
         s4D<ValueType> s;
         for(IndexType i=0, nt=0; i<NMO; i++) {
           for(IndexType k=0; k<NMO; k++, nt++) {
             if(nt<tg_ik0 || nt>=tg_ik1) continue;
             s = std::make_tuple(i,kk,k,ii,ValueType(0));
             bool cjgt = find_smallest_permutation(s);
             if( std::get<0>(s) >= min_i && std::get<0>(s) < max_i ) {
               s4Dit it_ = std::lower_bound( V2.begin(), V2.end(), s, mySort);
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
         myComm->allreduce(Lcomm);
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
       mymax = std::make_tuple(-1,-1,0);
       for(IndexType i=0,ik=0,nt=0; i<NMO; i++) {
        for(IndexType k=0; k<NMO; k++,nt++) {
         if(nt<ik0 || nt>=ik1) continue;
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
        if(nt>=ik1) break;
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
       if(rnk==0) app_log()<<" -- Time to test Cholesky factorization: " <<Timer.average("Generic") <<"\n";

     }

      /********************************************************************
      *  You get 2 potentials per Cholesky std::vector   
      *
      *    vn(+-)_{i,k} = sum_n 0.5*( L^n_{i,k} +- conj(L^n_{k,i}) )            
      ********************************************************************/

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
#if defined(QMC_COMPLEX)
          MPI_Gather(L[n].data(),2*maxnterms,MPI_DOUBLE,Lcomm.data(),2*maxnterms,MPI_DOUBLE,0,myComm->getMPI());
#else
          MPI_Gather(L[n].data(),maxnterms,MPI_DOUBLE,Lcomm.data(),maxnterms,MPI_DOUBLE,0,myComm->getMPI());
#endif
          if(rank()==0) Ls = Lcomm.data();
        } else {
          Ls = L[n].data();
        } 
        if(rank()==0) {
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
          myComm->barrier();
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
          if(rank()==0) {
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
      if(parallel) myComm->barrier();

      Timer.stop("Generic2");
      if(rnk==0) app_log()<<"     -- setup: " <<Timer.average("Generic2") <<"\n";

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
         //MPI_Allgather(L[n].data(),2*maxnterms,MPI_DOUBLE,Lcomm.data(),2*maxnterms,MPI_DOUBLE,myComm->getMPI());
         MPI_Gather(L[n].data(),2*maxnterms,MPI_DOUBLE,Lcomm.data(),2*maxnterms,MPI_DOUBLE,0,myComm->getMPI());
         if(head_of_nodes)
           MPI_Bcast(Lcomm.data(),2*Lcomm.size(),MPI_DOUBLE,0,MPI_COMM_HEAD_OF_NODES);   
#else
         MPI_Allgather(L[n].data(),maxnterms,MPI_DOUBLE,Lcomm.data(),maxnterms,MPI_DOUBLE,myComm->getMPI());
#endif
         Ls = Lcomm.data();
       } else {
         Ls = L[n].data();
       }
       Timer.stop("Generic2");
       Timer.start("Generic3");
       if(head_of_nodes && 2*n>=cv0 && 2*n<cvN && cnt_per_vec[2*n]>0) {
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
       if(head_of_nodes && (2*n+1)>=cv0 && (2*n+1)<cvN && cnt_per_vec[2*n+1]>0) {
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
         myComm->barrier(); 
       Timer.stop("Generic3");
      }
      app_log()<<"     -- av comm time: " <<Timer.average("Generic2") <<"\n";
      app_log()<<"     -- av insert time: " <<Timer.average("Generic3") <<"\n";

      Timer.stop("Generic");
      app_log()<<" -- Time to assemble Cholesky Matrix: " <<Timer.average("Generic") <<"\n";

      if(rank()==0 && nnodes>1 && parallel) { 
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

        if(parallel) Spvn.compress(TG.getNodeCommLocal());
        else if(head_of_nodes) Spvn.compress();
        //if(head_of_nodes) Spvn.compress();

        Timer.stop("Generic");
        app_log()<<"Done Compressing Spvn. \n";
        if(rnk==0) app_log()<<" -- Time to Compress Cholesky Matrix: " <<Timer.average("Generic") <<"\n";

      }

      if(parallel) myComm->barrier();

    } else {

     // way behind in development, needs parallel and dense options 
     APP_ABORT("Error: Cholesky factorization for unrestricted disabled. Finish implementation. \n\n\n");

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

      ValueType sqrtdt = std::sqrt(dt)*0.5;

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

#ifndef QMC_COMPLEX
      if(cntm>0) 
        APP_ABORT("Error: Found negative cholesky vectors in REAL build. \n");
#endif

      cnt=0;
      for(int n=0; n<L.size(); n++) { 
       int np=0;
       // v+
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) { 
         if(std::abs( (L[n][i*NMO+k] + myconj(L[n][k*NMO+i])) ) > cut) { 
           Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(sqrtdt*(L[n][i*NMO+k] + myconj(L[n][k*NMO+i]))));
           ++np;
         }
         if(std::abs( (L[n][NMO*NMO+i*NMO+k] + myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) {
           Spvn.add(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(sqrtdt*(L[n][NMO*NMO+i*NMO+k] + myconj(L[n][NMO*NMO+k*NMO+i]))));
           ++np;
         }
        }
       if(np>0) ++cnt;
#if defined(QMC_COMPLEX)
       np=0;
       // v-
       for(IndexType i=0; i<NMO; i++)
        for(IndexType k=0; k<NMO; k++) { 
         if(std::abs( (L[n][i*NMO+k] - myconj(L[n][k*NMO+i])) ) > cut) { 
           Spvn.add(i*NMO+k,cnt,static_cast<SPValueType>(ComplexType(0,1.0)*sqrtdt*(L[n][i*NMO+k] - myconj(L[n][k*NMO+i]))));
           ++np;
         }
         if(std::abs( (L[n][NMO*NMO+i*NMO+k] - myconj(L[n][NMO*NMO+k*NMO+i])) ) > cut) {
           Spvn.add(NMO*NMO+i*NMO+k,cnt,static_cast<SPValueType>(ComplexType(0,1.0)*sqrtdt*(L[n][NMO*NMO+i*NMO+k] - myconj(L[n][NMO*NMO+k*NMO+i]))));
           ++np;
         }
        }
       if(np>0) ++cnt;
#endif
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
      SPValueType* vals = Spvn.values();

      int NMO2 = spinRestricted?(NMO*NMO):(2*NMO*NMO);

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
      double scl = spinRestricted?(1.0):(4.0);
      app_log()<<"\n ********************************************\n Average error due to truncated eigenvalue factorization (in units of cutoff), max error : " <<s/cutoff_cholesky/NMO/NMO/NMO/NMO/scl <<" " <<max <<" \n ********************************************\n"<<std::endl;

       Timer.stop("Generic");
       if(rnk==0) app_log()<<" -- Time to test eigenvalue factorization: " <<Timer.average("Generic") <<"\n";
     }

  }

  bool SparseGeneralHamiltonian::parse(xmlNodePtr cur)
  {

    app_log()<<"\n\n --------------- Parsing Hamiltonian input ------------------ \n\n";

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
    std::string str5("yes"); 
    std::string str6("no"); 
    ParameterSet m_param;
    m_param.add(order,"orderStates","std::string");    
    m_param.add(cutoff1bar,"cutoff_1bar","double");
    m_param.add(cutoff_cholesky,"cutoff_decomp","double");
    m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
    m_param.add(cutoff_cholesky,"cutoff_factorization","double");
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(filetype,"filetype","std::string");
    m_param.add(filename,"filename","std::string");
    m_param.add(hdf_write_file,"hdf_write_file","std::string");
    m_param.add(hdf_write_type,"hdf_write_type","std::string");
    m_param.add(ascii_write_file,"ascii_write_file","std::string");
    m_param.add(rotation,"rotation","std::string");
    m_param.add(number_of_TGs,"nblocks","int");
    m_param.add(bkp,"test_breakup","std::string");
    m_param.add(str1,"printEig","std::string");
    m_param.add(str2,"test_2eint","std::string");
    m_param.add(str3,"fix_2eint","std::string");
    m_param.add(str4,"test_algo","std::string");
    m_param.add(str5,"inplace","std::string");
    m_param.add(str6,"skip_V2","std::string");
    m_param.add(n_reading_cores,"num_io_cores","int");
    m_param.put(cur);

    orderStates=false;
    std::transform(order.begin(),order.end(),order.begin(),(int (*)(int))tolower);
    std::transform(filetype.begin(),filetype.end(),filetype.begin(),(int (*)(int))tolower);
    std::transform(hdf_write_type.begin(),hdf_write_type.end(),hdf_write_type.begin(),(int (*)(int))tolower);
    std::transform(bkp.begin(),bkp.end(),bkp.begin(),(int (*)(int))tolower);
    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int))tolower);
    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int))tolower);
    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int))tolower);
    std::transform(str4.begin(),str4.end(),str4.begin(),(int (*)(int))tolower);
    std::transform(str5.begin(),str5.end(),str5.begin(),(int (*)(int))tolower);
    std::transform(str6.begin(),str6.end(),str6.begin(),(int (*)(int))tolower);
    if(order == "yes" || order == "true") orderStates = true;  
    if(bkp == "yes" || bkp == "true") test_breakup = true;  
    if(str1 == "yes" || str1 == "true") printEig = true;  
    if(str2 == "yes" || str2 == "true") test_2eint = true;  
    if(str3 == "yes" || str3 == "true") zero_bad_diag_2eints = true;  
    if(str4 == "no" || str4 == "false") test_algo = false;  
    if(str5 == "no" || str5 == "false") inplace = false;  
    if(str6 == "yes" || str6 == "true") skip_V2 = true;  

    if(skip_V2) 
      app_log()<<" Skipping 2 electron integrals. Only correct if other objects are initialized from hdf5 files. \n";
   
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
  //  For complex, there are 3 extra non-symmetric terms:
  //  J1a = <il|kj>
  //  J2a = <ik|lj>
  //  J3a = <il|jk> 
  //    In this case, make sure you use the fact that: <ij|kl> = conj( <kl|ij> )
  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_closed_shell(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v)
  {
    // simple algorithm for now
    // 1. add all contributions blindly
    // 2. check for repeated and remove
    // 3. apply symmetry elimination
   
    if(aa_only)
      APP_ABORT(" Error in SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_closed_shell(). aa_only is not allowed. \n\n\n");

    v.reserve(24);  

#ifndef QMC_COMPLEX 
    J1a=J1;
    J2a=J2;
    J3a=J3;
#endif

    // <ij||kl> -> (ijkl)
    // Symmetries:
    //   (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)* 

    ValueType J1J2 = 4.0*J1-2.0*J2;
    if(std::abs(J1J2) > cut) {
      push_ijkl(i,j,k,l,J1J2,v);
      push_ijkl(k,l,i,j,myconj(J1J2),v);
    }
    ValueType J1J3 = 4.0*J1a-2.0*J3a;
    if(std::abs(J1J3) > cut) {
      push_ijkl(i,l,k,j,J1J3,v);
      push_ijkl(j,k,l,i,myconj(J1J3),v);
    }

    ValueType J2J1 = 4.0*J2-2.0*J1;
    if(std::abs(J2J1) > cut) {
      push_ijkl(i,j,l,k,J2J1,v);
      push_ijkl(k,l,j,i,myconj(J2J1),v);
    }
    ValueType J2J3 = 4.0*J2a-2.0*J3;
    if(std::abs(J2J3) > cut) {
      push_ijkl(i,k,l,j,J2J3,v);
      push_ijkl(j,l,k,i,myconj(J2J3),v);
    }

    ValueType J3J1 = 4.0*J3a-2.0*J1a;
    if(std::abs(J3J1) > cut) {
      push_ijkl(i,l,j,k,J3J1,v);
      push_ijkl(k,j,l,i,myconj(J3J1),v);
    }
    ValueType J3J2 = 4.0*J3-2.0*J2a;
    if(std::abs(J3J2) > cut) {
      push_ijkl(i,k,j,l,J3J2,v);
      push_ijkl(l,j,k,i,myconj(J3J2),v);
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

  // For a given quartet ijkl and the 3 (for real) different permutations 
  // in the reduced list, add all possible contributions to the Hamiltonian
  // taking into account cutoff and (ik)<->(jl) symmetry  
  // J1 = <ij|kl>
  // J2 = <ij|lk>
  // J3 = <ik|jl> 
  //  For complex, there are 3 extra non-symmetric terms:
  //  J1a = <il|kj>
  //  J2a = <ik|lj>
  //  J3a = <il|jk> 
  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_spinRestricted(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v) 
  { 
    // simple algorithm for now
    // 1. add all contributions blindly
    // 2. check for repeated and remove
    // 3. apply symmetry elimination
    // 4. if aa_only, only alpha/alpha component is added 
    //
   
    v.reserve(48);  
    v.clear();
  
    int i2 = i+NMO; 
    int j2 = j+NMO; 
    int k2 = k+NMO; 
    int l2 = l+NMO; 
#ifndef QMC_COMPLEX
    J1a=J1;
    J2a=J2;
    J3a=J3;
#else
      APP_ABORT("Not yet working for std::complex \n\n\n");
#endif
    // 2bar terms
    ValueType J1J2 = J1-J2;
    if(std::abs(J1J2) > cut) {
      push_ijkl(i,j,k,l,J1J2,v);
      push_ijkl(k,l,i,j,J1J2,v);  
      if(!aa_only) push_ijkl(i2,j2,k2,l2,J1J2,v);
      if(!aa_only) push_ijkl(k2,l2,i2,j2,J1J2,v);  
    }
    ValueType J1J3 = J1-J3;
    if(std::abs(J1J3) > cut) {
      push_ijkl(i,l,k,j,J1J3,v);  
      push_ijkl(j,k,l,i,J1J3,v);  
      if(!aa_only) push_ijkl(i2,l2,k2,j2,J1J3,v);  
      if(!aa_only) push_ijkl(j2,k2,l2,i2,J1J3,v);  
    }

    ValueType J2J1 = J2-J1;
    if(std::abs(J2J1) > cut) {
      push_ijkl(i,j,l,k,J2J1,v);
      push_ijkl(k,l,j,i,J2J1,v);
      if(!aa_only) push_ijkl(i2,j2,l2,k2,J2J1,v);
      if(!aa_only) push_ijkl(k2,l2,j2,i2,J2J1,v);
    }
    ValueType J2J3 = J2-J3;
    if(std::abs(J2J3) > cut) {
      push_ijkl(i,k,l,j,J2J3,v);
      push_ijkl(j,l,k,i,J2J3,v);
      if(!aa_only) push_ijkl(i2,k2,l2,j2,J2J3,v);
      if(!aa_only) push_ijkl(j2,l2,k2,i2,J2J3,v);
    }
       
    ValueType J3J1 = J3-J1;
    if(std::abs(J3J1) > cut) {
      push_ijkl(i,l,j,k,J3J1,v);
      push_ijkl(k,j,l,i,J3J1,v);
      if(!aa_only) push_ijkl(i2,l2,j2,k2,J3J1,v);
      if(!aa_only) push_ijkl(k2,j2,l2,i2,J3J1,v);
    }    
    ValueType J3J2 = J3-J2;
    if(std::abs(J3J2) > cut) {
      push_ijkl(i,k,j,l,J3J2,v);
      push_ijkl(l,j,k,i,J3J2,v);
      if(!aa_only) push_ijkl(i2,k2,j2,l2,J3J2,v);
      if(!aa_only) push_ijkl(l2,j2,k2,i2,J3J2,v);
    }    

    // 1bar terms
    if(std::abs(J1) > cut && !aa_only) {
      push_ijkl(i,j2,k,l2,J1,v);
      push_ijkl(k,l2,i,j2,J1,v);
      push_ijkl(i,l2,k,j2,J1,v);
      push_ijkl(j,k2,l,i2,J1,v);
      push_ijkl(j,i2,l,k2,J1,v);
      push_ijkl(l,k2,j,i2,J1,v);
      push_ijkl(l,i2,j,k2,J1,v);
      push_ijkl(k,j2,i,l2,J1,v);
    }       
    if(std::abs(J2) > cut && !aa_only) {
      push_ijkl(i,j2,l,k2,J2,v);
      push_ijkl(k,l2,j,i2,J2,v);
      push_ijkl(i,k2,l,j2,J2,v);
      push_ijkl(j,l2,k,i2,J2,v);
      push_ijkl(j,i2,k,l2,J2,v);
      push_ijkl(l,k2,i,j2,J2,v);
      push_ijkl(k,i2,j,l2,J2,v);
      push_ijkl(l,j2,i,k2,J2,v);
    }    
    if(std::abs(J3) > cut && !aa_only) {
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

  // For a given quartet ijkl and the 3 (for real) different permutations 
  // in the reduced list, add all possible contributions to the Hamiltonian
  // taking into account cutoff and (ik)<->(jl) symmetry  
  // J1 = <ij|kl>
  // J2 = <ij|lk>
  // J3 = <ik|jl> 
  //  For complex, there are 3 extra non-symmetric terms:
  //  J1a = <il|kj>
  //  J2a = <ik|lj>
  //  J3a = <il|jk> 
  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_ghf(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v) 
  { 
    // simple algorithm for now
    // 1. add all contributions blindly
    // 2. check for repeated and remove
    // 3. apply symmetry elimination
    //
    // GHF: all terms from _spinRestricted routine plus: M_(i-alpha k-beta)_(j-beta l-alpha) = <ij|lk> 
  
    if(aa_only)
      APP_ABORT(" Error in SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_ghf(). aa_only is not allowed. \n\n\n");
 
    v.reserve(48);  
    v.clear();
  
    int i2 = i+NMO; 
    int j2 = j+NMO; 
    int k2 = k+NMO; 
    int l2 = l+NMO; 
#ifndef QMC_COMPLEX
    J1a=J1;
    J2a=J2;
    J3a=J3;
#else
      APP_ABORT("Not yet working for std::complex \n\n\n");
#endif
    // 2bar terms
    ValueType J1J2 = J1-J2;
    if(std::abs(J1J2) > cut) {
      push_ijkl(i,j,k,l,J1J2,v,true);
      push_ijkl(k,l,i,j,J1J2,v,true);  
      push_ijkl(i2,j2,k2,l2,J1J2,v,true);
      push_ijkl(k2,l2,i2,j2,J1J2,v,true);  
    }
    ValueType J1J3 = J1-J3;
    if(std::abs(J1J3) > cut) {
      push_ijkl(i,l,k,j,J1J3,v,true);  
      push_ijkl(j,k,l,i,J1J3,v,true);  
      push_ijkl(i2,l2,k2,j2,J1J3,v,true);  
      push_ijkl(j2,k2,l2,i2,J1J3,v,true);  
    }

    ValueType J2J1 = J2-J1;
    if(std::abs(J2J1) > cut) {
      push_ijkl(i,j,l,k,J2J1,v,true);
      push_ijkl(k,l,j,i,J2J1,v,true);
      push_ijkl(i2,j2,l2,k2,J2J1,v,true);
      push_ijkl(k2,l2,j2,i2,J2J1,v,true);
    }
    ValueType J2J3 = J2-J3;
    if(std::abs(J2J3) > cut) {
      push_ijkl(i,k,l,j,J2J3,v,true);
      push_ijkl(j,l,k,i,J2J3,v,true);
      push_ijkl(i2,k2,l2,j2,J2J3,v,true);
      push_ijkl(j2,l2,k2,i2,J2J3,v,true);
    }
       
    ValueType J3J1 = J3-J1;
    if(std::abs(J3J1) > cut) {
      push_ijkl(i,l,j,k,J3J1,v,true);
      push_ijkl(k,j,l,i,J3J1,v,true);
      push_ijkl(i2,l2,j2,k2,J3J1,v,true);
      push_ijkl(k2,j2,l2,i2,J3J1,v,true);
    }    
    ValueType J3J2 = J3-J2;
    if(std::abs(J3J2) > cut) {
      push_ijkl(i,k,j,l,J3J2,v,true);
      push_ijkl(l,j,k,i,J3J2,v,true);
      push_ijkl(i2,k2,j2,l2,J3J2,v,true);
      push_ijkl(l2,j2,k2,i2,J3J2,v,true);
    }    

    // 1bar terms: <alpha,beta|alpha,beta>
    if(std::abs(J1) > cut) {
      push_ijkl(i,j2,k,l2,J1,v,true);
      push_ijkl(k,l2,i,j2,J1,v,true);
      push_ijkl(i,l2,k,j2,J1,v,true);
      push_ijkl(j,k2,l,i2,J1,v,true);
      push_ijkl(j,i2,l,k2,J1,v,true);
      push_ijkl(l,k2,j,i2,J1,v,true);
      push_ijkl(l,i2,j,k2,J1,v,true);
      push_ijkl(k,j2,i,l2,J1,v,true);
    }       
    if(std::abs(J2) > cut) {
      push_ijkl(i,j2,l,k2,J2,v,true);
      push_ijkl(k,l2,j,i2,J2,v,true);
      push_ijkl(i,k2,l,j2,J2,v,true);
      push_ijkl(j,l2,k,i2,J2,v,true);
      push_ijkl(j,i2,k,l2,J2,v,true);
      push_ijkl(l,k2,i,j2,J2,v,true);
      push_ijkl(k,i2,j,l2,J2,v,true);
      push_ijkl(l,j2,i,k2,J2,v,true);
    }    
    if(std::abs(J3) > cut) {
      push_ijkl(i,l2,j,k2,J3,v,true);
      push_ijkl(k,j2,l,i2,J3,v,true);
      push_ijkl(i,k2,j,l2,J3,v,true);
      push_ijkl(l,j2,k,i2,J3,v,true);
      push_ijkl(l,i2,k,j2,J3,v,true);
      push_ijkl(j,k2,i,l2,J3,v,true);
      push_ijkl(k,i2,l,j2,J3,v,true);
      push_ijkl(j,l2,i,k2,J3,v,true);
    }

    // GHF extra exchage terms: <alpha,beta|beta,alpha>
    // M_iakb_jbla = <iajb|lakb> 
    // J1 = <ij|kl>
    // J2 = <ij|lk> 
    // J3 = <ik|jl> 
    J1*=RealType(-1); 
    J2*=RealType(-1); 
    J3*=RealType(-1); 
    if(std::abs(J1) > cut) {
      push_ijkl(i,j2,l2,k,J1,v,true);
      push_ijkl(k,j2,l2,i,J1,v,true);
      push_ijkl(i,l2,j2,k,J1,v,true);
      push_ijkl(k,l2,j2,i,J1,v,true);
      push_ijkl(j,i2,k2,l,J1,v,true);
      push_ijkl(l,i2,k2,j,J1,v,true);
      push_ijkl(j,k2,i2,l,J1,v,true);
      push_ijkl(l,k2,i2,j,J1,v,true);
    }       
    if(std::abs(J2) > cut) {
      push_ijkl(i,j2,k2,l,J2,v,true);
      push_ijkl(l,j2,k2,i,J2,v,true);
      push_ijkl(i,k2,j2,l,J2,v,true);
      push_ijkl(l,k2,j2,i,J2,v,true);
      push_ijkl(j,i2,l2,k,J2,v,true);
      push_ijkl(k,i2,l2,j,J2,v,true);
      push_ijkl(j,l2,i2,k,J2,v,true);
      push_ijkl(k,l2,i2,j,J2,v,true);
    }    
    if(std::abs(J3) > cut) {
      push_ijkl(i,k2,l2,j,J3,v,true);
      push_ijkl(j,k2,l2,i,J3,v,true);
      push_ijkl(i,l2,k2,j,J3,v,true);
      push_ijkl(j,l2,k2,i,J3,v,true);
      push_ijkl(k,i2,j2,l,J3,v,true);
      push_ijkl(l,i2,j2,k,J3,v,true);
      push_ijkl(k,j2,i2,l,J3,v,true);
      push_ijkl(l,j2,i2,k,J3,v,true);
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

  // For a given quartet ijkl and the 3 (for real) different permutations 
  // in the reduced list, add all possible contributions to the Hamiltonian
  // taking into account cutoff and (ik)<->(jl) symmetry  
  // J1 = <ij|kl>
  // J2 = <ij|lk>
  // J3 = <ik|jl> 
  //  For complex, there are 3 extra non-symmetric terms:
  //  J1a = <il|kj>
  //  J2a = <ik|lj>
  //  J3a = <il|jk> 
  void SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_general(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v) 
  {
    APP_ABORT("Finsigh implementation. \n\n\n");

    if(aa_only)
      APP_ABORT(" Error in SparseGeneralHamiltonian::find_all_contributions_to_hamiltonian_general(). aa_only is not allowed. \n\n\n");
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
      if( klij != ijkl && klij != jilk ) v.push_back(std::make_tuple(k,l,i,j,myconj(V)));   
      if( lkji != ijkl && lkji != jilk && lkji != klij ) v.push_back(std::make_tuple(l,k,j,i,myconj(V)));   
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
          v.push_back(std::make_tuple(k,l,i,j,myconj(static_cast<RealType>(2.0)*V)));   
          if(ijkl == lkji)  {
            std::cerr<<" Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl == lkji): " <<i <<" " <<j <<" " <<k <<" " <<l <<"  " <<V  <<std::endl;
            APP_ABORT("Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl == lkji) \n"); 
          }
        } else 
          v.push_back(std::make_tuple(k,l,i,j,myconj(V)));

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

  bool SparseGeneralHamiltonian::find_smallest_permutation(s4D<ValueType>& val) {

#if defined(QMC_COMPLEX)
    // jl < ik
    if(  std::forward_as_tuple(std::get<1>(val),std::get<3>(val) ) < std::forward_as_tuple(std::get<0>(val),std::get<2>(val) )  ) {
        std::swap(std::get<0>(val),std::get<1>(val));
        std::swap(std::get<2>(val),std::get<3>(val));
    }
    // kl < ij
    if(  std::forward_as_tuple(std::get<2>(val),std::get<3>(val) ) < std::forward_as_tuple(std::get<0>(val),std::get<1>(val) )  ) {
      std::swap(std::get<0>(val),std::get<2>(val));
      std::swap(std::get<1>(val),std::get<3>(val));
      std::get<4>(val) = std::conj(std::get<4>(val));
      // jl < ik again since ij<->kl swap occured  
      if(  std::forward_as_tuple(std::get<1>(val),std::get<3>(val) ) < std::forward_as_tuple(std::get<0>(val),std::get<2>(val) )  ) {
        std::swap(std::get<0>(val),std::get<1>(val));
        std::swap(std::get<2>(val),std::get<3>(val));
      }
      return true;
    } else {
      // only possibility is that l < i, since I know that the current i is smaller than j and k
      // 
      if( std::forward_as_tuple(std::get<3>(val),std::get<2>(val) ) < std::forward_as_tuple(std::get<0>(val),std::get<1>(val) )  ) { 
        std::swap(std::get<0>(val),std::get<3>(val));
        std::swap(std::get<2>(val),std::get<1>(val));
        std::get<4>(val) = std::conj(std::get<4>(val));
        return true;
      }
      return false; 
    }
#else
    // i < k
    if( std::get<2>(val) < std::get<0>(val) ) std::swap(std::get<0>(val),std::get<2>(val) );    
    // j < l
    if( std::get<3>(val) < std::get<1>(val) ) std::swap(std::get<1>(val),std::get<3>(val) );    
    // ik < jl
    if( std::get<1>(val) < std::get<0>(val) ) { 
        std::swap(std::get<0>(val),std::get<1>(val) ); 
        std::swap(std::get<2>(val),std::get<3>(val) ); 
    } else if( (std::get<1>(val) == std::get<0>(val)) && ( std::get<3>(val) < std::get<2>(val) )) 
        std::swap(std::get<2>(val),std::get<3>(val) );    
    return false;
#endif

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

  // generates all symmetry inequivalent integrals V2(indx,j,k,l) for a given indx.
  // Vijkl is only meaningful for myComm->rank()==0.
  bool SparseGeneralHamiltonian::SparseHamiltonianFromFactorization( int indx, std::vector<OrbitalType>& jkl, std::vector<ValueType>& intgs, const RealType cut)
  {

    if(!spinRestricted) {
      app_error()<<" Error: SparseGeneralHamiltonian::FullHamiltonianFromFactorization only implemented for spin restricted integrals. \n";
      APP_ABORT("Error: SparseGeneralHamiltonian::FullHamiltonianFromFactorization only implemented for spin restricted integrals. \n");
    }

    if(!factorizedHamiltonian) {
      app_error()<<" Error:  Calling SparseGeneralHamiltonian::FullHamiltonianFromFactorization without factorized hamiltonian. \n"; 
      APP_ABORT("Error:  Calling SparseGeneralHamiltonian::FullHamiltonianFromFactorization without factorized hamiltonian. \n");
    }

    ValueType zero = ValueType(0);

#if defined(QMC_COMPLEX)
    APP_ABORT(" Error: SparseGeneralHamiltonian::SparseHamiltonianFromFactorization doesn't yet work for std::complex matrix elements. \n");
#endif

    intgs.reserve(NMO*NMO);  // reasonable guess
    intgs.clear();
    jkl.reserve(3*NMO*NMO);
    jkl.clear();

    long cnter=0, npr = myComm->size(), rk = myComm->rank();
    ValueType J1;

    if(DiagHam.size() == 0) {
      DiagHam.resize(NMO,NMO);
      for(OrbitalType i=0; i<NMO; i++)
      for(OrbitalType k=i; k<NMO; k++, cnter++) {
        if( cnter%npr != rk ) continue;
        DiagHam(i,k) =  H(i,k,k,i);
        DiagHam(k,i) = DiagHam(i,k); 
#if defined(QMC_COMPLEX)
        if(DiagHam(i,k).imag() > 1e-8) {
            app_error()<<" Error: Found complex diagonal on hamiltonian. " <<i <<" " <<k <<" " <<DiagHam(i,k) <<std::endl;
            APP_ABORT("Error: Found complex diagonal on hamiltonian.");
        }
#endif
      }
      myComm->allreduce(DiagHam);      
    }  

#if defined(QMC_COMPLEX)
// don't know how to do this quickly in complex case
    auto ijkl = std::make_tuple(indx,0,0,0); 
    for(OrbitalType j=indx; j<NMO; j++)  {
    for(OrbitalType k=indx; k<NMO; k++)  {
    for(OrbitalType l=indx; l<NMO; l++, cnter++)  {
      if( cnter%npr != rk ) continue;

      // |<ij|kl>| <= sqrt( <ik|ki> * <lj|jl>  )  
      if( std::sqrt( std::abs(DiagHam(indx,k)*DiagHam(l,j)) ) > cut ) {
        std::get<1>(ijkl)=j;
        std::get<2>(ijkl)=k;
        std::get<3>(ijkl)=l;
        // am I the smallest equivalent permutation?
        if( ijkl <= std::forward_as_tuple(j,indx,l,k)  &&
            ijkl <= std::forward_as_tuple(k,l,indx,j) &&
            ijkl <= std::forward_as_tuple(l,k,j,indx) ) { 
          J1 = H(indx,j,k,l);
          if(std::abs(J1) > cut) {
             jkl.push_back(j);
             jkl.push_back(k);
             jkl.push_back(l);
             intgs.push_back(J1);
          }
        }
      }

    }
    }
    }
#else
    for(OrbitalType j=indx; j<NMO; j++)  {
    for(OrbitalType k=indx; k<NMO; k++)  {
      OrbitalType l0 = (j==indx)?k:j;  
      for(OrbitalType l=l0; l<NMO; l++, cnter++)  {
        if( cnter%npr != rk ) continue;

        // |<ij|kl>| <= sqrt( <ik|ki> * <lj|jl>  )  
        if( std::sqrt( std::abs(DiagHam(indx,k)*DiagHam(l,j)) ) > cut ) {
          J1 = H(indx,j,k,l);
          if(std::abs(J1) > cut) {
             jkl.push_back(j);
             jkl.push_back(k);
             jkl.push_back(l);
             intgs.push_back(J1);
          }
        }
      
      }
    }
    }
#endif

    std::vector<int> lnint(1);
    std::vector<int> gnint, disp;
    lnint[0] = intgs.size();
    if(myComm->rank()==0) {
      gnint.resize(npr);
    }
    myComm->gather(lnint,gnint,0);

    if(myComm->rank() == 0) {

      int ntot = 0;
      for(int i=0; i<gnint.size(); i++) ntot += gnint[i];

      disp.resize(npr);
      disp[0]=disp[1]=0;
      for(int i=2; i<npr; i++) 
        disp[i] = disp[i-1] + gnint[i-1]; 

      int nint0 = gnint[0];
      gnint[0]=0; // no need to gather with myself, this way I avoid the need for a second vector
      intgs.resize(ntot);
      myComm->gatherv(intgs.data(), intgs.data()+nint0, 0, gnint, disp, 0, myComm->getMPI());
      for(int i=2; i<npr; i++) disp[i] *= 3;
      for(int i=1; i<npr; i++) gnint[i] *= 3;
      jkl.resize(3*ntot);
      myComm->gatherv(jkl.data(), jkl.data()+3*nint0, 0, gnint, disp, 0, myComm->getMPI());

    } else {
      myComm->gatherv(intgs.data(), intgs.data(), intgs.size(), gnint, disp, 0, myComm->getMPI());
      myComm->gatherv(jkl.data(), jkl.data(), jkl.size(), gnint, disp, 0, myComm->getMPI());
      intgs.clear();
      jkl.clear();
    }
  }

  bool SparseGeneralHamiltonian::createHamiltonianForPureDeterminant(int walker_type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b, std::vector<s1D<ValueType> >& hij, SPValueSMSpMat& Vijkl, const RealType cut)
  {

    // walker_type: 0-closed_shell density matrix, 1-ROHF/UHF density matrix, 2-GHF density matrix
   

    //  For alpha-alpha and beta-beta store two-bar integrals directly
    //  For alpha-beta and beta-alpha, store one bar integrals 
    //
    // Symmetries for real orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij> = <lk||ji>  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji> = -<lk||ij> 
    //
    //                                    
    //  For one-bar integrals: <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
    //                                 = <ji|lk> = <li|jk> = <jk|li> = <lk|ji> 
    //
    // Symmetries for Complex orbitals:
    //  For two-bar integrals: <ij||kl> = <ji||lk> = <kl||ij>* = <lk||ji>*  
    //                                  = -<ij||lk> = -<ji||kl> = -<kl||ji>* = -<lk||ij>* 
    //  For one-bar integrals: <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
    //  Notice that in this case:   <ij|kl> != <kj|il> and other permutations            
    // 

    if(skip_V2) 
      APP_ABORT("Error: Calling SparseGeneralHamiltonian routines with skip_V2=yes. \n\n\n");

#ifdef AFQMC_DEBUG
    app_log()<<" In SparseGeneralHamiltonian :: createHamiltonianForPureDeterminant." <<std::endl; 
#endif

    if(!spinRestricted && walker_type==2) {
      APP_ABORT("Error: GHF density matrix only implemented with spinRestricted integrals. \n");
    }

    hij.clear();

    bool closed_shell = walker_type==0;

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
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+i , ValueType(2.0)*V );
          } else  {  // ij/ji both (alpha/alpha) and (beta/beta)
            if(occ_a[i]) *ith++ = std::make_tuple( i*NMO+j , ValueType(2.0)*V );
            if(occ_a[j]) *ith++ = std::make_tuple( j*NMO+i , ValueType(2.0)*myconj(V) );
          }
        } else {
          if( i == j ) {  // ii alpha / ii beta 
            if(occ_a[i]) *ith++ = std::make_tuple( Index2Mat(i,i,walker_type==2) , V );
            if(occ_b[i+NMO]) *ith++ = std::make_tuple(Index2Mat(i+NMO,i+NMO,walker_type==2), V );
          } else  {  // ij/ji both (alpha/alpha) and (beta/beta)
            if(occ_a[i]) *ith++ = std::make_tuple( Index2Mat(i,j,walker_type==2), V );
            if(occ_a[j]) *ith++ = std::make_tuple( Index2Mat(j,i,walker_type==2), myconj(V) );
            if(occ_b[i+NMO]) *ith++ = std::make_tuple( Index2Mat(i+NMO,j+NMO,walker_type==2) , V );
            if(occ_b[j+NMO]) *ith++ = std::make_tuple( Index2Mat(j+NMO,i+NMO,walker_type==2) , myconj(V) );
          }
        }
      } else {
        if( i == j ) {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j,walker_type==2), V );
        } else {
          if(occ_a[i] || occ_b[i]) *ith++ = std::make_tuple( Index2Mat(i,j,walker_type==2), V );
          if(occ_a[j] || occ_b[j]) *ith++ = std::make_tuple( Index2Mat(j,i,walker_type==2), myconj(V) );
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

        Timer.reset("Generic");
        Timer.start("Generic");

        ValueType zero = ValueType(0);

        std::vector<s4D<ValueType> > vs4D;
        vs4D.reserve(24);

        cnt2=0; 
        long cnter=0, npr = myComm->size(), rk = myComm->rank();
        OrbitalType i,j,k,l,j1,k1,l1,j2,k2,l2;
        ValueType J1,J2,J3,J1a=zero,J2a=zero,J3a=zero,fct;

        if(!spinRestricted)
          APP_ABORT("Error: Need to finish implementation for !spinRestricted. \n\n\n");

        DiagHam.resize(NMO,NMO); 
        for(IndexType i=0; i<NMO; i++) 
        for(IndexType k=i; k<NMO; k++, cnter++) { 
          if( cnter%npr != rk ) continue;
          DiagHam(i,k) =  H(i,k,k,i);
          DiagHam(k,i) = DiagHam(i,k); 
#if defined(QMC_COMPLEX)
          if(DiagHam(i,k).imag() > 1e-8) {
              app_error()<<" Error: Found complex diagonal on hamiltonian. " <<i <<" " <<k <<" " <<DiagHam(i,k) <<std::endl;
              APP_ABORT("Error: Found complex diagonal on hamiltonian.");
          }
#endif
        }
        myComm->allreduce(DiagHam);
/*
        if( isComplex(Diag(0,0)) )
          MPI_Allreduce(MPI_IN_PLACE, Diag.data(),  2*Diag.size(), MPI_DOUBLE, MPI_SUM,
                myComm->getMPI());
        else
          MPI_Allreduce(MPI_IN_PLACE, Diag.data(),  Diag.size(), MPI_DOUBLE, MPI_SUM,
                myComm->getMPI());
*/

        int occi, occj, occk, occl; 
        cnter=0;
        // Approximation: (similar to non factorized case)
        //   - if <ij|kl>  <  cutoff1bar, set to zero 
        for(IndexType i=0; i<NMO; i++)  {
        for(IndexType j=i; j<NMO; j++,cnter++)  {
         if( cnter%npr != rk ) continue;
         occi = (occ_a[i]||occ_b[i])?1:0;
         occj = (occ_a[j]||occ_b[j])?1:0;
        for(IndexType k=j; k<NMO; k++)  {
         occk = (occ_a[k]||occ_b[k])?1:0;
        for(IndexType l=k; l<NMO; l++)  {

            occl = (occ_a[l]||occ_b[l])?1:0;
            if( occi + occj + occk + occl < 2 ) continue; 

            // NOTE NOTE NOTE:
            // <ik|ik> can get a complex component due to truncation error, eliminate it here

            J1=J2=J3=zero; 
            // J1 = <ij|kl>   
            // J1 < sqrt( <ik|ki> * <lj|jl> )  
            if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) { 
              J1 = H(i,j,k,l);  
#if defined(QMC_COMPLEX)
              if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
            }

            // J2 = <ij|lk> 
            if(i==j || l==k) {
              J2=J1; 
            } else { 
              if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) {
                J2 = H(i,j,l,k);  
#if defined(QMC_COMPLEX)
                if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
              }  
            }

            // J3 = <ik|jl> 
            if(j==k) {
              J3=J1; 
            } else if(i==l) {
              J3 = myconj(J1);
            } else { 
              if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) { 
                J3 = H(i,k,j,l);  
#if defined(QMC_COMPLEX)
                if(i==j && k==l) J3 = ValueType(J3.real(),0); 
#endif
              }
            }

#if defined(QMC_COMPLEX)
            //  J2a = <ik|lj>
            if(l==j) {
              J2a=J3;
            } else if(i==k) {
              J2a=J3;
            } else if(k==j) {
              J2a=J2;
            } else if(i==l) {
              J2a=std::conj(J2);
            } else {
              if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
                J2a = H(i,k,l,j);
            }
            
            //  J3a = <il|jk> 
            if(l==j) {
              J3a=J2;
            } else if(i==k) {
              J3a=std::conj(J2);
            } else if(k==l) {
              J3a=J3;
            } else if(i==j) {
              J3a=std::conj(J3);
            } else {
              if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
                J3a = H(i,l,j,k);
            }

            //  For complex, there are 3 extra non-symmetric terms:
            //  J1a = <il|kj>
            if(k==l) {
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
              if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
                J1a = H(i,l,k,j);
            }
#endif



            if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue; 

            vs4D.clear();
            if(walker_type==0) {
              find_all_contributions_to_hamiltonian_closed_shell(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D); 
            } else if(walker_type==1) {
              if(spinRestricted) 
                find_all_contributions_to_hamiltonian_spinRestricted(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D); 
              else 
                find_all_contributions_to_hamiltonian_general(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D); 
            } else if(walker_type==2) {
                find_all_contributions_to_hamiltonian_ghf(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D); 
            } else {
              APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n"); 
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
        //   - if <ij|kl>  <  cutoff1bar, set to zero 
        for(IndexType i=0; i<NMO; i++)  {
        for(IndexType j=i; j<NMO; j++,cnter++)  {
         if( cnter%npr != rk ) continue;
         occi = (occ_a[i]||occ_b[i])?1:0;
         occj = (occ_a[j]||occ_b[j])?1:0;
        for(IndexType k=j; k<NMO; k++)  {
         occk = (occ_a[k]||occ_b[k])?1:0;
        for(IndexType l=k; l<NMO; l++)  {

            occl = (occ_a[l]||occ_b[l])?1:0;
            if( occi + occj + occk + occl < 2 ) continue; 

            J1=J2=J3=zero; 
            // J1 = <ij|kl>   
            // J1 < sqrt( <ik|ki> * <lj|jl> )  
            if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(l,j)) ) > cutoff1bar ) { 
              J1 = H(i,j,k,l);  
#if defined(QMC_COMPLEX)
              if(i==k && j==l) J1 = ValueType(J1.real(),0);
#endif
            }

            // J2 = <ij|lk> 
            if(i==j || l==k) {
              J2=J1; 
            } else { 
              if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(k,j)) ) > cutoff1bar ) { 
                J2 = H(i,j,l,k); 
#if defined(QMC_COMPLEX)
                if(i==l && j==k) J2 = ValueType(J2.real(),0);
#endif
              } 
            }

            // J3 = <ik|jl> 
            if(j==k) {
              J3=J1; 
            } else if(i==l) {
              J3 = myconj(J1);
            } else { 
              if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(l,k)) ) > cutoff1bar ) {
                J3 = H(i,k,j,l);  
#if defined(QMC_COMPLEX)
                if(i==j && k==l) J3 = ValueType(J3.real(),0);
#endif
              }
            }

#if defined(QMC_COMPLEX)
            //  J2a = <ik|lj>
            if(l==j) {
              J2a=J3;
            } else if(i==k) {
              J2a=J3;
            } else if(k==j) {
              J2a=J2;
            } else if(i==l) {
              J2a=std::conj(J2);
            } else {
              if( std::sqrt( std::abs(DiagHam(i,l)*DiagHam(j,k)) ) > cutoff1bar )
                J2a = H(i,k,l,j);
            }
            
            //  J3a = <il|jk> 
            if(l==j) {
              J3a=J2;
            } else if(i==k) {
              J3a=std::conj(J2);
            } else if(k==l) {
              J3a=J3;
            } else if(i==j) {
              J3a=std::conj(J3);
            } else {
              if( std::sqrt( std::abs(DiagHam(i,j)*DiagHam(k,l)) ) > cutoff1bar )
                J3a = H(i,l,j,k);
            }

            //  For complex, there are 3 extra non-symmetric terms:
            //  J1a = <il|kj>
            if(k==l) {
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
              if( std::sqrt( std::abs(DiagHam(i,k)*DiagHam(j,l)) ) > cutoff1bar )
                J1a = H(i,l,k,j);
            }
            
#endif

            if( std::abs(J1)<cutoff1bar && std::abs(J2)<cutoff1bar && std::abs(J3)<cutoff1bar && std::abs(J1a)<cutoff1bar && std::abs(J2a)<cutoff1bar && std::abs(J3a)<cutoff1bar) continue; 

            vs4D.clear();
            if(walker_type==0) {
              find_all_contributions_to_hamiltonian_closed_shell(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
            } else if(walker_type==1) {
              if(spinRestricted)
                find_all_contributions_to_hamiltonian_spinRestricted(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
              else
                find_all_contributions_to_hamiltonian_general(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
            } else if(walker_type==2) {
                find_all_contributions_to_hamiltonian_ghf(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
            } else {
              APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
            }
            cnt2+=add_allowed_terms(vs4D,occ_a,occ_b, Vijkl, true, walker_type==2);                    
        }
        }
        }
        }
        myComm->barrier();
        Timer.stop("Generic");
        app_log()<<"Time to generate full Hamiltonian from factorized form: " <<Timer.total("Generic") <<std::endl; 
 

    } else {  //factorizedHamiltonian 

// for parallel algorithm with distributed matrices, count the number of terms
// in the hamiltonian and divide evenly among nodes on the TG.
// Then the cores on the local node work on the local segment 

      Timer.reset("Generic");
      Timer.start("Generic");

      std::vector<s4D<ValueType> > vs4D;  
      vs4D.reserve(48);
      long N = spinRestricted?NMO:2*NMO;
      // this should already exist, just in case
      //if(IJ.size() == 0) 
        //generateIJ();

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
      long npr = myComm->size(), rk = myComm->rank();
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
            find_all_contributions_to_hamiltonian_closed_shell(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==1) {
            if(spinRestricted)
              find_all_contributions_to_hamiltonian_spinRestricted(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
            else
              find_all_contributions_to_hamiltonian_general(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==2) {
              find_all_contributions_to_hamiltonian_ghf(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else {
            APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
          }
          cnt2+=count_allowed_terms(vs4D,occ_a,occ_b);                    
        }
      }

      myComm->allreduce(cnt2);
long tmp_ = cnt2;
      Vijkl.allocate(cnt2+1000);
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
            find_all_contributions_to_hamiltonian_closed_shell(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==1) {
            if(spinRestricted)
              find_all_contributions_to_hamiltonian_spinRestricted(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
            else
              find_all_contributions_to_hamiltonian_general(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else if(walker_type==2) {
              find_all_contributions_to_hamiltonian_ghf(aa_only,i,j,k,l,J1,J2,J3,J1a,J2a,J3a,cut,vs4D);
          } else {
            APP_ABORT(" Error: Unknown walker type in createHamiltonianForPureDeterminant. \n");
          }
          cnt2+=add_allowed_terms(vs4D,occ_a,occ_b, Vijkl, true, walker_type==2);                    
        }
      }

      myComm->barrier();
      Timer.stop("Generic");
      app_log()<<"Time to generate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;
  
    } // factrorizedHamiltonnian

    Timer.reset("Generic");
    Timer.start("Generic");
    if(!communicate_Vijkl(Vijkl)) return false;
    Timer.stop("Generic");
    app_log()<<"Time to communicate Hamiltonian: " <<Timer.total("Generic") <<std::endl; 

    Timer.reset("Generic");
    Timer.start("Generic");

    if(!Vijkl.remove_repeated_and_compress(TG.getNodeCommLocal())) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

    Timer.stop("Generic");
    app_log()<<"Time to remove_repeated_and_compress Hamiltonian: " <<Timer.total("Generic") <<std::endl;

    return true;  

  }

  bool SparseGeneralHamiltonian::communicate_Vijkl(SPValueSMSpMat& Vijkl) 
  {

    myComm->barrier();

    if(head_of_nodes) {

      int rk,npr;
      long ptr,n0;
      MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk);
      MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr);
      n0 = Vijkl.size(); // my number of terms, always from zero to n0
      ptr = n0; // position to copy elements to 
      std::vector<long> size(npr);
      size[rk] = n0;
      myComm->gsum(size,MPI_COMM_HEAD_OF_NODES);
      long ntot = 0;
      for(int i=0; i<npr; i++) ntot+=size[i];

app_log()<<"size: " <<ntot <<" " <<n0 <<std::endl;

      if(ntot > Vijkl.capacity()) {
        app_error()<<" Problems gathering hamiltonian. Capacity of std::vector is not sufficient: " <<ntot <<" " <<Vijkl.capacity() <<" \n";
        return false;
      }

app_log()<<" before resize: " <<std::endl;
      Vijkl.resize_serial(ntot);
app_log()<<" after resize: " <<std::endl;

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
app_log()<<" after bcasts: " <<std::endl;
    }
    myComm->barrier();
    return true; 
  }

#ifndef QMC_COMPLEX
  bool SparseGeneralHamiltonian::communicate_Vijkl(SPComplexSMSpMat& Vijkl)
  {

    myComm->barrier();

    if(head_of_nodes) {

      int rk,npr;
      long ptr,n0;
      MPI_Comm_rank(MPI_COMM_HEAD_OF_NODES,&rk);
      MPI_Comm_size(MPI_COMM_HEAD_OF_NODES,&npr);
      n0 = Vijkl.size(); // my number of terms, always from zero to n0
      ptr = n0; // position to copy elements to 
      std::vector<long> size(npr);
      size[rk] = n0;
      myComm->gsum(size,MPI_COMM_HEAD_OF_NODES);
      long ntot = 0;
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
#endif

  bool SparseGeneralHamiltonian::createHamiltonianForGeneralDeterminant (int walker_type, const ComplexMatrix& A,std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vijkl, const RealType cut) 
  {

    ComplexMatrix M,N;
    const ComplexType one = ComplexType(1.0);
    const ComplexType zero = ComplexType(0.0);
    int npr = myComm->size(), rk = myComm->rank();

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
        M(i,j) = ValueType(2.0)*V; //
        if( i!=j ) M(j,i) = ValueType(2.0)*myconj(V);
      }
      DenseMatrixOperators::product_AhB(NAEA,NMO,NMO,one,A.data(),NAEA,M.data(),NMO,zero,N.data(),NMO);
      int cnt=0;
/*
      for(int i=0; i<NMO; i++)
       for(int j=0; j<NMO; j++)
         std::cout<<i <<" " <<j <<" " <<M(i,j) <<"\n";
      for(int i=0; i<NAEA; i++)
       for(int j=0; j<NMO; j++)
         std::cout<<i <<" " <<j <<" " <<N(i,j) <<"\n";
      for(int i=0; i<NMO; i++)
       for(int j=0; j<NAEA; j++)
         std::cout<<i <<" " <<j <<" " <<A(i,j) <<"\n";
*/
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
        APP_ABORT("Error: TODO: createHamiltonianForGeneralDeterminant not implemented for GHF walker type yet. \n\n\n");
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

      int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
      int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
  
      // <ab||kl> = sum_n Qk(a,n) * Rl(b,n) - Rl(a,n)*Qk(b,n),
      // where:
      //   Qk(a,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
      //   Rl(a,n) = sum_i conj(Amat(i,a)) * conj(V2_fact(li,n))
      // For real build, Qk=Rk
      //
      // For parallelization, distribute (k,l) pairs over nodes.
      // Build ahead of time Qk/Rl matrices in shared memory to reduce memory/setup time.
      // Assemble integrals in parallel and fully distributed.
      // Collect on all nodes.
      //    - For distributed hamiltonians, you do not need to worry about keeping contiguous
      //    segments of the hamiltonian. Only that the distribution over nodes is roughly equal.
      //    Write a simple algorithm that balances the number of terms in a TG. 
      //  

      Timer.reset("Generic");
      Timer.start("Generic");
//      SPComplexSMVector::pointer Qkptr;
//      SPComplexSMVector::pointer Rlptr;
//      SPComplexSMVector tQk;    
//      tQk.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_tQk"),TG.getNodeCommLocal());  
//      SPComplexSMVector Qk;    
//      Qk.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_Qk"),TG.getNodeCommLocal());  
      SPComplexSMVector Rl;    
      Rl.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_Rl"),TG.getNodeCommLocal());  

      SPComplexSMSpMat SptQk;
      SptQk.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_SptQk"),TG.getNodeCommLocal());
      SPComplexSMSpMat SpQk;
      SpQk.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_SpQk"),TG.getNodeCommLocal());
      SPComplexSMSpMat SpRl;
      SpRl.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_SpRl"),TG.getNodeCommLocal());

      if(useSpMSpM) {
        APP_ABORT(" Error: Not fully implemented. \n\n\n");
      }  

      int NMO2 = (walker_type == 0)?NMO:2*NMO;
      int NAEA2 = (walker_type == 0)?NAEA:2*NAEA;
      int ngrp = std::min(NMO2,nnodes);
      std::vector<int> M_split(ngrp+1);

      // split the orbitals among processors 
      FairDivide(NMO2,ngrp,M_split);

      // Construct your set of Q(k,a,m), R(l,b,m)
      int l0 = (nodeid<ngrp)?(M_split[nodeid]):(-1);
      int lN = (nodeid<ngrp)?(M_split[nodeid+1]):(-1);
      bool amIAlpha = true; // for simplicity, subset of bands must be either elpha or beta, not mixed
      if( l0 < NMO && (lN-1) < NMO )
        amIAlpha = true;
      else if( l0 >= NMO && lN >= NMO )
        amIAlpha = false;
      else {
        std::cerr<<"l0, lN, nodeid, ngrp, NMO: " <<l0 <<" " <<lN <<" " <<nodeid <<" " <<ngrp <<" " <<NMO <<std::endl;  
        std::cerr<<" Error: Current algorithm requires an even number of processors. \n\n\n";  
        APP_ABORT(" Error: Current algorithm requires an even number of processors. \n\n\n");  
      }
      int norb = lN-l0;
      int maxnorb = 0;
      for(int i=0; i<ngrp; i++) maxnorb = std::max(maxnorb,M_split[i+1]-M_split[i]);  
      int nvec = V2_fact.cols();
      // Rl(k, a, m), k:[0:NMO2], a:[0:NAEA], m:[0:nvec] 
      int mat_size = norb * NAEA * nvec; 
      if(nodeid < ngrp ) {
        //Qk.resize(mat_size);  
        if(!useSpMSpM)
          Rl.resize(mat_size);  
      }
      if(coreid==0 && !useSpMSpM) std::fill(Rl.begin(),Rl.end(),SPComplexType(0.0));
//      if(coreid==0) std::fill(Qk.begin(),Qk.end(),SPComplexType(0.0));

      int bl0=-1, blN=-1;
      int nwork = std::min(norb*NAEA,ncores);  
      if(coreid <nwork)
        if(amIAlpha)
          std::tie(bl0, blN) = FairDivideBoundary(coreid,norb*NAEA,nwork);
        else
          std::tie(bl0, blN) = FairDivideBoundary(coreid,norb*NAEB,nwork);
      // right now this is an upper bound
      int nak = blN-bl0;

      SpQk.setDims(norb*NAEA,nvec);
      if(useSpMSpM) 
        SpRl.setDims(nvec,norb*NAEA);

      std::vector<SPComplexType> vec(nvec);
      std::vector<std::tuple<SPValueSMSpMat::intType,SPValueSMSpMat::intType,SPComplexType>> abkl;
      int nmax = 1000000;
      abkl.reserve(nmax);  

      if(distribute_Ham) {
        APP_ABORT(" Finish THIS (43)!!! \n\n\n");
      } else {

/*
        //   Q(k,a,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
        //   R(l,a,n) = sum_i conj(Amat(i,a)) * conj(V2_fact(li,n))
        if(walker_type == 0) {

          long sz=0;
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            int n_ = (l<NMO)?0:NMO;
            int NEL = (l<NMO)?NAEA:NAEB;
            for(int a=0; a<NEL; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
              for(int i=0; i<NMO; i++) {
                int il = i*NMO+l-n_;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(il));
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(il+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                for(int n=0; n<nv; n++)
                  vec[itLcol[n]] += Aiac * itLval[n];
              }
              for(int n=0; n<nvec; n++)
                if(std::abs(vec[n])>cut) sz++;
            }
          }
          {
            long c_ = sz;
            MPI_Allreduce(&c_,&sz,1,MPI_LONG,MPI_SUM,TG.getNodeCommLocal());
          }
          SpQk.allocate(sz);
          SpQk.barrier();
          abkl.clear();
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            int n_ = (l<NMO)?0:NMO;
            int NEL = (l<NMO)?NAEA:NAEB;
            for(int a=0; a<NEL; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
              for(int i=0; i<NMO; i++) {
                int il = i*NMO+l-n_;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(il));
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(il+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                for(int n=0; n<nv; n++)
                  vec[itLcol[n]] += Aiac * itLval[n];
              }
              for(int n=0; n<nvec; n++)
                if(std::abs(vec[n])>cut) {
                  abkl.push_back( std::make_tuple(cnt, n, vec[n]) );
                  if(abkl.size()==nmax) {
                    SpQk.add(abkl,true);
                    abkl.clear();
                  }
                }
            }
          }
          if(abkl.size()>0) {
            SpQk.add(abkl,true);
            abkl.clear();
          }
          SpQk.compress(TG.getNodeCommLocal())

#if defined(QMC_COMPLEX)
// Think about how to calculate Rl quickly since now it is (nvec,norb*NAEA)
// If it is not possible to evaluate quickly, then calculate it in Qk with (norb*NAEA,nvec)
// form and then transpose it into Rl
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            for(int a=0; a<NAEA; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
              for(int i=0; i<NMO; i++) {
                int li = l*NMO+i;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(li));  
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(li+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;  
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol =  V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i,a));
                for(int n=0; n<nv; n++)   
                  vec[itLcol[n]] += Aiac * std::conj(itLval[n]); 
                //Rlptr = Rl.values() + nt*NAEA*nvec + a*nvec; 
                //for(int n=0; n<nv; n++)   
                //  Rlptr[itLcol[n]] += Aiac * std::conj(itLval[n]); 
              }
              BLAS::axpy(nvec,SPComplexType(1.0),vec.data(),1,Rl.values()+cnt,norb*NAEA);
            }    
          }
#else
          if(useSpMSpM) {
            SpRl.allocate(SpQk.size());
            if(coreid==0) {
              std::copy(SpQk.values(),SpQk.values()+SpQk.size(),SpRl.values());
              std::copy(SpQk.column_data(),SpQk.column_data()+SpQk.size(),SpRl.row_data());
              std::copy(SpQk.row_data(),SpQk.row_data()+SpQk.size(),SpRl.column_data());
            }
            SpRl.compress(TG.getNodeCommLocal());
            SpRl.barrier();
          } else {
            if(coreid==0)
              std::fill(Rl.begin(),Rl.end(),SPComplexType(0));
            Rl.barrier();
            int n0_,n1_;
            assert(SpQk.size() >= ncores);
            std::tie(n0_, n1_) = FairDivideBoundary(coreid,SpQk.size(),ncores);
            SPComplexSMSpMat::iterator val = SpQk.vals_begin()+n0_;
            SPComplexSMSpMat::iterator vend = SpQk.vals_begin()+n1_;
            SPComplexSMSpMat::int_iterator col = SpQk.cols_begin()+n0_;
            SPComplexSMSpMat::int_iterator row = SpQk.rows_begin()+n0_;
            int ncol = SpQk.rows();
            while(val != vend)
              Rl[ (*(col++))*ncol + (*(row++)) ] = *(val++);
            Rl.barrier();
          }
#endif
        } else if(walker_type == 1) {  

          assert(spinRestricted);
*/
/*
          // Construct Qk[k,n,nvec]
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            int n_ = (l<NMO)?0:NMO;
            int NEL = (l<NMO)?NAEA:NAEB;
            for(int a=0; a<NEL; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              for(int i=0; i<NMO; i++) {
                int il = i*NMO+l-n_;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(il));  
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(il+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;  
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                Qkptr = Qk.values() + nt*NEL*nvec + a*nvec; 
                for(int n=0; n<nv; n++)   
                  Qkptr[itLcol[n]] += Aiac * itLval[n]; 
              }
            }    
          }
*/ 

        if(walker_type == 0 || walker_type == 1) {
          // Construct SpQk[k,n,nvec]
          long sz=0;
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            int n_ = (l<NMO)?0:NMO;
            int NEL = (l<NMO)?NAEA:NAEB;
            for(int a=0; a<NEL; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
              for(int i=0; i<NMO; i++) {
                int il = i*NMO+l-n_;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(il));
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(il+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                for(int n=0; n<nv; n++)
                  vec[itLcol[n]] += Aiac * itLval[n];
              }
              for(int n=0; n<nvec; n++)
                if(std::abs(vec[n])>cut) sz++;     
            }
          }
          {
            long c_ = sz; 
            MPI_Allreduce(&c_,&sz,1,MPI_LONG,MPI_SUM,TG.getNodeCommLocal());  
          }
          SpQk.allocate(sz);
          SpQk.barrier();
          abkl.clear();
          for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
            int n_ = (l<NMO)?0:NMO;
            int NEL = (l<NMO)?NAEA:NAEB;
            for(int a=0; a<NEL; a++, cnt++) {
              if( cnt%ncores != coreid ) continue;
              std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
              for(int i=0; i<NMO; i++) {
                int il = i*NMO+l-n_;
                SPValueSMSpMat::intType n0 = *(V2_fact.row_index(il));
                SPValueSMSpMat::intType n1 = *(V2_fact.row_index(il+1));
                if(n1==n0) continue;
                SPValueSMSpMat::intType nv = n1-n0;
                const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                for(int n=0; n<nv; n++)
                  vec[itLcol[n]] += Aiac * itLval[n];
              }
              for(int n=0; n<nvec; n++)
                if(std::abs(vec[n])>cut) { 
                  abkl.push_back( std::make_tuple(cnt, n, vec[n]) );
                  if(abkl.size()==nmax) {
                    SpQk.add(abkl,true);
                    abkl.clear();
                  }
                }
            }
          }
          if(abkl.size()>0) {
            SpQk.add(abkl,true);
            abkl.clear();
          }
          SpQk.compress(TG.getNodeCommLocal());

#if defined(QMC_COMPLEX)
          if(!useSpMSpM) {
            for(int k=l0, nt=0, cnt=0; k<lN; k++, nt++) {
              int n_ = (k<NMO)?0:NMO;
              int NEL = (k<NMO)?NAEA:NAEB;
              for(int a=0; a<NEL; a++, cnt++) {
                if( cnt%ncores != coreid ) continue;
                std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
                for(int i=0; i<NMO; i++) {
                  int ki = (k-n_)*NMO+i;
                  SPValueSMSpMat::intType n0 = *(V2_fact.row_index(ki));  
                  SPValueSMSpMat::intType n1 = *(V2_fact.row_index(ki+1));
                  if(n1==n0) continue;
                  int nv = n1-n0;  
                  const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                  const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                  const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                  for(int n=0; n<nv; n++)   
                    vec[itLcol[n]] += Aiac * std::conj(itLval[n]); 
                }
                BLAS::axpy(nvec,SPComplexType(1.0),vec.data(),1,Rl.values()+cnt,norb*NEL);
              }    
            }
            Rl.barrier();
          } else {  

            // Construct SpRl[nvec,k,n]
            sz=0;
            for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
              int n_ = (l<NMO)?0:NMO;
              int NEL = (l<NMO)?NAEA:NAEB;
              for(int a=0; a<NEL; a++, cnt++) {
                if( cnt%ncores != coreid ) continue;
                std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
                for(int i=0; i<NMO; i++) {
                  int li = (l-n_)*NMO+i;
                  SPValueSMSpMat::intType n0 = *(V2_fact.row_index(li));
                  SPValueSMSpMat::intType n1 = *(V2_fact.row_index(li+1));
                  if(n1==n0) continue;
                  SPValueSMSpMat::intType nv = n1-n0;
                  const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                  const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                  const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                  for(int n=0; n<nv; n++)
                    vec[itLcol[n]] += Aiac * std::conj(itLval[n]);
                }
                for(int n=0; n<nvec; n++)
                  if(std::abs(vec[n])>cut) sz++;
              }
            }
            {
              long c_ = sz;
              MPI_Allreduce(&c_,&sz,1,MPI_LONG,MPI_SUM,TG.getNodeCommLocal());
            }
            SpRl.allocate(sz);
            SpRl.barrier();
            abkl.clear();
            for(int l=l0, nt=0, cnt=0; l<lN; l++, nt++) {
              int n_ = (l<NMO)?0:NMO;
              int NEL = (l<NMO)?NAEA:NAEB;
              for(int a=0; a<NEL; a++, cnt++) {
                if( cnt%ncores != coreid ) continue;
                std::fill(vec.begin(),vec.end(),SPComplexType(0,0));
                for(int i=0; i<NMO; i++) {
                  int li = (l-n_)*NMO+i;
                  SPValueSMSpMat::intType n0 = *(V2_fact.row_index(li));
                  SPValueSMSpMat::intType n1 = *(V2_fact.row_index(li+1));
                  if(n1==n0) continue;
                  SPValueSMSpMat::intType nv = n1-n0;
                  const SPValueSMSpMat::pointer itLval = V2_fact.values(n0);
                  const SPValueSMSpMat::intPtr itLcol = V2_fact.column_data(n0);
                  const ComplexType Aiac = std::conj(A(i+n_,a));  // alpha/beta
                  for(int n=0; n<nv; n++)
                    vec[itLcol[n]] += Aiac * std::conj(itLval[n]);
                }
                for(int n=0; n<nvec; n++)
                  if(std::abs(vec[n])>cut) {
                    abkl.push_back( std::make_tuple(n, cnt, vec[n]) );
                    if(abkl.size()==nmax) {
                      SpRl.add(abkl,true);
                      abkl.clear();
                    }
                  }
              }
            }
            if(abkl.size()>0) {
              SpRl.add(abkl,true);
              abkl.clear();
            }
            SpRl.compress(TG.getNodeCommLocal());
          }
#else
          if(useSpMSpM) {
            SpRl.allocate(SpQk.size());
            if(coreid==0) {
              std::copy(SpQk.values(),SpQk.values()+SpQk.size(),SpRl.values());
              std::copy(SpQk.column_data(),SpQk.column_data()+SpQk.size(),SpRl.row_data());
              std::copy(SpQk.row_data(),SpQk.row_data()+SpQk.size(),SpRl.column_data());
            }
            SpRl.compress(TG.getNodeCommLocal());
            SpRl.barrier();
          } else {
            if(coreid==0)
              std::fill(Rl.begin(),Rl.end(),SPComplexType(0));
            Rl.barrier();
            int n0_,n1_,sz_=SpQk.size();
            assert(SpQk.size() >= ncores); 
            std::tie(n0_, n1_) = FairDivideBoundary(coreid,sz_,ncores);   
            SPComplexSMSpMat::iterator val = SpQk.vals_begin()+n0_;
            SPComplexSMSpMat::iterator vend = SpQk.vals_begin()+n1_;
            SPComplexSMSpMat::int_iterator col = SpQk.cols_begin()+n0_; 
            SPComplexSMSpMat::int_iterator row = SpQk.rows_begin()+n0_; 
            int ncol = SpQk.rows();
            while(val != vend) 
              Rl[ (*(col++))*ncol + (*(row++)) ] = *(val++);
            Rl.barrier();   
          }
#endif
        } else if(walker_type == 2) {  
          APP_ABORT(" Error in createHamiltonianForGeneralDeterminant: GHF not implemented. \n\n\n");
        } else {
          APP_ABORT(" Error in createHamiltonianForGeneralDeterminant: Unknown walker_type. \n\n\n");
        }

      }  
      SpQk.barrier();

      // let the maximum message be 1 GB
      int maxnt = std::max(1,static_cast<int>(std::floor(1024.0*1024.0*1024.0/sizeof(SPComplexType))));

      std::vector<int> nkbounds;          // local bounds for communication  
      std::vector<int> Qknum(nnodes);     // number of blocks per node
      std::vector<int> Qksizes;           // number of terms and number of k vectors in block for all nodes
      if(head_of_nodes) {
        int n_=0, ntcnt=0, n0=0; 
        int NEL = (amIAlpha)?NAEA:NAEB;
        for(int i=0; i<norb; i++) {
          int ntt = *SpQk.row_index((i+1)*NEL) - *SpQk.row_index(i*NEL); 
          assert(ntt < maxnt);  
          if(ntcnt+ntt > maxnt) {
            nkbounds.push_back(ntcnt);
            nkbounds.push_back(i-n0);
            ntcnt=ntt;
            n0=i;
            n_++;
          } else {
            ntcnt+=ntt;
          }   
        }
        if(ntcnt > 0) {
          // push last block
          n_++;
          nkbounds.push_back(ntcnt);
          nkbounds.push_back(norb-n0);
        }

        MPI_Allgather(&n_,1,MPI_INT,Qknum.data(),1,MPI_INT,MPI_COMM_HEAD_OF_NODES);

        int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0); 
        Qksizes.resize(2*ntt);

        std::vector<int> cnts(nnodes);
        std::vector<int> disp(nnodes);
        int cnt=0;
        for(int i=0; i<nnodes; i++) {
          cnts[i] = Qknum[i]*2;
          disp[i]=cnt;
          cnt+=cnts[i];  
        }
        MPI_Allgatherv(nkbounds.data(),nkbounds.size(),MPI_INT,Qksizes.data(),cnts.data(),disp.data(),MPI_INT,MPI_COMM_HEAD_OF_NODES );

      }

      MPI_Bcast(Qknum.data(),nnodes,MPI_INT,0,TG.getNodeCommLocal());
      int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0); 
      if(!head_of_nodes)  
        Qksizes.resize(2*ntt);
      MPI_Bcast(Qksizes.data(),Qksizes.size(),MPI_INT,0,TG.getNodeCommLocal());

// store {nterms,nk} for all nodes 
// use it to know communication pattern 

      int maxnk = 0;        // maximum number of k vectors in communication block 
      long maxqksize = 0;   // maximum size of communication block
      for(int i=0; i<ntt; i++) {
        if(Qksizes[2*i] > maxqksize) maxqksize = Qksizes[2*i];
        if(Qksizes[2*i+1] > maxnk) maxnk = Qksizes[2*i+1];
      }    

      SPComplexSMVector Ta;
      Ta.setup(head_of_nodes,std::string("SparseGeneralHamiltonian_Ta"),TG.getNodeCommLocal());
      Ta.resize(norb*NAEA*maxnk*NAEA);
      Ta.barrier();
       
      // setup working sparse matrix  
      SptQk.setDims(maxnk * NAEA, nvec);
      SptQk.allocate(maxqksize+1000);

      abkl.clear();  
      std::vector<SPComplexSMSpMat::intType> rowI;
      rowI.reserve( (maxnk*NAEA/ncores) + 1 );  

      Timer.stop("Generic");
      app_log()<<"Time to construct distributed rotated Cholesky vectors: " <<Timer.total("Generic") <<std::endl;

      Timer.reset("Generic");
      Timer.start("Generic");

      std::size_t sz=0;  
      std::size_t sz1=0, sz2=0, sz3=0;  

      // now calculate fully distributed matrix elements
      for(int nn=0, nb=0, nkcum=0; nn<ngrp; nn++) {

        // just checking
        assert(nkcum==M_split[nn]);
        if(M_split[nn+1]==M_split[nn]) continue;
        int nblk = Qknum[nn]; 
        long ntermscum=0;
        for( int bi = 0; bi < nblk; bi++, nb++) {
          int nterms = Qksizes[2*nb];      // number of terms in block 
          int nk = Qksizes[2*nb+1];        // number of k-blocks in block
          int k0 = nkcum;                  // first value of k in block
          nkcum+=nk;                       
          int kN = nkcum;                  // last+1 value
          int NEL0 = (k0<NMO)?NAEA:NAEB;   // number of electrons in this spin block
          assert(nk > 0 && nk <= maxnk );  // just checking

/*
Timer.reset("Generic2");
Timer.start("Generic2");
          if(head_of_nodes) {
            if(nn == nodeid) {
              std::copy(Qk.values()+bi*maxnk*NEL0*nvec,Qk.values()+(bi*maxnk+nk)*NEL0*nvec,tQk.values());
            }  
#if defined(AFQMC_SP)
            MPI_Bcast(tQk.values(),nk*NEL0*nvec*2,MPI_FLOAT,nn,MPI_COMM_HEAD_OF_NODES);
#else
            MPI_Bcast(tQk.values(),nk*NEL0*nvec*2,MPI_DOUBLE,nn,MPI_COMM_HEAD_OF_NODES);
#endif
          } 
          Qk.barrier();

Timer.stop("Generic2");
app_log()<<"Time to Bcast tQk: " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;
*/

Timer.reset("Generic2");
Timer.start("Generic2");

          if(head_of_nodes) {
            if(nn == nodeid) {
              long n0 = *SpQk.row_index( (k0-M_split[nn])*NEL0 );
              long n1 = *SpQk.row_index( (k0-M_split[nn]+nk)*NEL0 );
              int nt_ = static_cast<int>(n1-n0);
              assert(ntermscum==n0);
              assert(nt_==nterms);
              SptQk.resize_serial(nterms);
              std::copy(SpQk.values()+n0,SpQk.values()+n1,SptQk.values());
              std::copy(SpQk.column_data()+n0,SpQk.column_data()+n1,SptQk.column_data());
              for(int i=0, j=(k0-M_split[nn])*NEL0; i<=nk*NEL0; i++, j++)
                SptQk.row_index()[i] = SpQk.row_index()[j]-n0;
            }
            MPI_Bcast(&nterms,1,MPI_LONG,nn,MPI_COMM_HEAD_OF_NODES);
            SptQk.resize_serial(nterms);
#if defined(AFQMC_SP)
            MPI_Bcast(SptQk.values(),nterms*2,MPI_FLOAT,nn,MPI_COMM_HEAD_OF_NODES);
#else
            MPI_Bcast(SptQk.values(),nterms*2,MPI_DOUBLE,nn,MPI_COMM_HEAD_OF_NODES);
#endif
            MPI_Bcast(SptQk.column_data(),nterms,MPI_INT,nn,MPI_COMM_HEAD_OF_NODES);
            MPI_Bcast(SptQk.row_index(),nk*NEL0+1,MPI_INT,nn,MPI_COMM_HEAD_OF_NODES);
          }
          SptQk.setCompressed();
          SptQk.barrier();

          // for safety, keep track of sum
          ntermscum += static_cast<long>(nterms);

Timer.stop("Generic2");
app_log()<<"Time to Bcast SptRl: " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;

Timer.reset("Generic2");
Timer.start("Generic2");

          if(walker_type == 0) {

            // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
            // Qa(ka,lb) = Q(k,a,:)*R(l,b,:)
            //DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEA,SPComplexType(0),Qa.values()+bl0,norb*NAEA);
            //Qa.barrier();
            SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),SptQk.values(),SptQk.column_data(), SptQk.row_index(), Rl.values()+bl0,norb*NAEA,SPComplexType(0),Ta.values()+bl0,norb*NAEA);
            SpQk.barrier();
            SPComplexType four = SPComplexType(4.0);
            SPComplexType two = SPComplexType(2.0);
            for(int k=k0, ka=0; k<kN; k++) {
              for(int a=0; a<NAEA; a++, ka++) {
                for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b  
                  int b = lb%NAEA;
                  if(b<a) continue;
                  int l = lb/NAEA+l0;
                  int la = (l-l0)*NAEA+a;  
                  int kb = (k-k0)*NAEA+b;  
                  SPComplexType qkalb = *(Ta.values() + ka*norb*NAEA + lb);  // Qa(ka,lb)   
                  SPComplexType qlakb = *(Ta.values() + kb*norb*NAEA + la);  // Qa(kb,la)   
                  if(std::abs( four*qkalb - two*qlakb ) > cut) sz++;
                }
              }
            }

          } else if(walker_type == 1) {

            if( M_split[nn] < NMO && (M_split[nn+1]-1) < NMO ) {
              // k is alpha
            
              if(amIAlpha) {  
                // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
                // Qa(ka,lb) = Q(k,a,:)*R(l,b,:)
//                DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEA,SPComplexType(0),Qa.values()+bl0,norb*NAEA);
//                Qa.barrier();
                SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),SptQk.values(),SptQk.column_data(), SptQk.row_index(), Rl.values()+bl0,norb*NAEA,SPComplexType(0),Ta.values()+bl0,norb*NAEA);
                Ta.barrier();
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEA; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b  
                      int b = lb%NAEA;
                      if(b<=a) continue;
                      int l = lb/NAEA+l0;
                      int la = (l-l0)*NAEA+a;
                      int kb = (k-k0)*NAEA+b;
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEA + lb);  // Qa(ka,lb)   
                      SPComplexType qlakb = *(Ta.values() + kb*norb*NAEA + la);  // Qa(kb,la)   
                      if(std::abs( qkalb - qlakb ) > cut) sz++;
                      if(std::abs( qkalb - qlakb ) > cut) sz1++;
                    }
                  }
                }
              } else {
                // <a,b | k,l> = Qa(ka,lb) = Q(k,a,:)*R(l,b,:) 
                //DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEB,SPComplexType(0),Qa.values()+bl0,norb*NAEB);
                SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),SptQk.values(),SptQk.column_data(), SptQk.row_index(), Rl.values()+bl0,norb*NAEB,SPComplexType(0),Ta.values()+bl0,norb*NAEB);
                Ta.barrier();
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEA; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEB+b  
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEB + lb);  // Qa(ka,lb)   
                      if(std::abs( qkalb ) > cut) sz++;
                      if(std::abs( qkalb ) > cut) sz2++;
                    }
                  }
                }
              }  

            } else if(M_split[nn] >= NMO && M_split[nn+1] >= NMO) {
              // k is beta           
              if(!amIAlpha) {
                // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
                // Qa(ka,lb) = Q(k,a,:)*R(l,b,:)                
                //DenseMatrixOperators::product(nk*NAEB,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEB,SPComplexType(0),Qa.values()+bl0,norb*NAEB);
                //Qa.barrier();
                SparseMatrixOperators::product_SpMatM(nk*NAEB,int(blN-bl0),nvec,SPComplexType(1.0),SptQk.values(),SptQk.column_data(), SptQk.row_index(), Rl.values()+bl0,norb*NAEB,SPComplexType(0),Ta.values()+bl0,norb*NAEB);
                Ta.barrier();
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEB; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEB+b  
                      int b = lb%NAEB;
                      if(b<=a) continue;
                      int l = lb/NAEB+l0;
                      int la = (l-l0)*NAEB+a;
                      int kb = (k-k0)*NAEB+b;
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEB + lb);  // Qa(ka,lb)   
                      SPComplexType qlakb = *(Ta.values() + kb*norb*NAEB + la);  // Qa(kb,la)   
                      if(std::abs( qkalb - qlakb ) > cut) sz++;
                      if(std::abs( qkalb - qlakb ) > cut) sz3++;
                    }
                  }
                }
              }  
            } else {
              APP_ABORT(" Error: This should not happen. \n\n\n");
            } 

          } else if(walker_type == 2) {
            APP_ABORT(" Error in createHamiltonianForGeneralDeterminant: GHF not implemented. \n\n\n");
          }
MPI_Barrier(myComm->getMPI()); // to measure actual time
Timer.stop("Generic2");
app_log()<<"Time to calculate Qab: " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;

        }

      }

      std::size_t sz_local=sz;
      MPI_Allreduce(&sz_local,&sz,1,MPI_UNSIGNED_LONG,MPI_SUM,myComm->getMPI()); 
      // not yet distributed, later on decide distribution here based on FairDivide
      // and allocate appropriately          
      app_log()<<"Number of terms in Vijkl: " <<sz <<" " <<sz*(sizeof(ValueType)+sizeof(int)*2)/1024.0/1024.0 <<" MB" <<std::endl;
      Vijkl.allocate(sz);

      std::size_t sz_1=0, sz_2=0, sz_3=0;

      // now calculate fully distributed matrix elements
      for(int nn=0, nb=0, nkcum=0; nn<ngrp; nn++) {

        assert(nkcum==M_split[nn]);
        if(M_split[nn+1]==M_split[nn]) continue;
        int nblk = Qknum[nn];
        long ntermscum=0;
        for( int bi = 0; bi < nblk; bi++, nb++) {
          int nterms = Qksizes[2*nb];      // number of terms in block 
          int nk = Qksizes[2*nb+1];        // number of k-blocks in block
          int k0 = nkcum;                  // first value of k in block
          nkcum+=nk;
          int kN = nkcum;                  // last+1 value
          int NEL0 = (k0<NMO)?NAEA:NAEB;   // number of electrons in this spin block
          assert(nk > 0 && nk <= maxnk );  // just checking

/*
          if(head_of_nodes) {
            if(nn == nodeid) {
              std::copy(Qk.values()+bi*maxnk*NEL0*nvec,Qk.values()+(bi*maxnk+nk)*NEL0*nvec,tQk.values());
            }  
#if defined(AFQMC_SP)
            MPI_Bcast(tQk.values(),nk*NEL0*nvec*2,MPI_FLOAT,nn,MPI_COMM_HEAD_OF_NODES);
#else
            MPI_Bcast(tQk.values(),nk*NEL0*nvec*2,MPI_DOUBLE,nn,MPI_COMM_HEAD_OF_NODES);
#endif
          } 
          tQk.barrier();

*/

Timer.reset("Generic2");
Timer.start("Generic2");

          assert(nblk==1);
          if(head_of_nodes) {
            if(nn == nodeid) {
              long n0 = *SpQk.row_index( (k0-M_split[nn])*NEL0 );
              long n1 = *SpQk.row_index( (k0-M_split[nn]+nk)*NEL0 );
              int nt_ = static_cast<int>(n1-n0);
              assert(ntermscum==n0);
              assert(nt_==nterms);
              SptQk.resize_serial(nterms);
              std::copy(SpQk.values()+n0,SpQk.values()+n1,SptQk.values());
              std::copy(SpQk.column_data()+n0,SpQk.column_data()+n1,SptQk.column_data());
              for(int i=0, j=(k0-M_split[nn])*NEL0; i<=nk*NEL0; i++, j++)
                SptQk.row_index()[i] = SpQk.row_index()[j]-n0;
            }
            MPI_Bcast(&nterms,1,MPI_LONG,nn,MPI_COMM_HEAD_OF_NODES);
            SptQk.resize_serial(nterms);
#if defined(AFQMC_SP)
            MPI_Bcast(SptQk.values(),nterms*2,MPI_FLOAT,nn,MPI_COMM_HEAD_OF_NODES);
#else
            MPI_Bcast(SptQk.values(),nterms*2,MPI_DOUBLE,nn,MPI_COMM_HEAD_OF_NODES);
#endif
            MPI_Bcast(SptQk.column_data(),nterms,MPI_INT,nn,MPI_COMM_HEAD_OF_NODES);
            MPI_Bcast(SptQk.row_index(),nk*NEL0+1,MPI_INT,nn,MPI_COMM_HEAD_OF_NODES);
          }
          SptQk.setCompressed();
          SptQk.barrier();

Timer.stop("Generic2");
app_log()<<"Time to Bcast SptQk: " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;

          // for safety, keep track of sum
          ntermscum += static_cast<long>(nterms);

//*/
          if(walker_type == 0) {

            // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
            // Qa(ka,lb) = 4*Q(k,a,:)*R(l,b,:)
            //DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEA,SPComplexType(0),Qa.values()+bl0,norb*NAEA);
            //Qa.barrier();
            SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),
                     SptQk.values(),SptQk.column_data(), SptQk.row_index(),
                     Rl.values()+bl0,norb*NAEA,SPComplexType(0),Ta.values()+bl0,norb*NAEA);
            Ta.barrier();
            SPComplexType four = SPComplexType(4.0);
            SPComplexType two = SPComplexType(2.0);
            for(int k=k0, ka=0; k<kN; k++) {
              for(int a=0; a<NAEA; a++, ka++) {
                for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b  
                  int b = lb%NAEA;
                  if(b<a) continue;
                  int l = lb/NAEA+l0;
                  int la = (l-l0)*NAEA+a;
                  int kb = (k-k0)*NAEA+b;
                  SPComplexType qkalb = *(Ta.values() + ka*norb*NAEA + lb);  // Qa(ka,lb)   
                  SPComplexType qlakb = *(Ta.values() + kb*norb*NAEA + la);  // Qa(kb,la)   
                  if(std::abs( four*qkalb - two*qlakb ) > cut) {
                    abkl.push_back( std::make_tuple(a*NMO+k, b*NMO+l, 2.0*(four*qkalb - two*qlakb)) );
                    if(abkl.size()==nmax) {
                      Vijkl.add(abkl,true);
                      abkl.clear();
                    } 
                  }
                }
              }
            }

          } else if(walker_type == 1) {

            if( M_split[nn] < NMO && (M_split[nn+1]-1) < NMO ) {
              // k is alpha

              if(amIAlpha) {
                // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
                // Qa(ka,lb) = Q(k,a,:)*R(l,b,:)
              
/*
                Qa.barrier();
Timer.reset("Generic2");
Timer.start("Generic2");

                DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEA,SPComplexType(0),Qa.values()+bl0,norb*NAEA);
                Qa.barrier();


Timer.stop("Generic2");
app_log()<<"Time for dense DGEMM : " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;
*/
/*
              SptQk.toOneBase();
              SpRl.toOneBase();
              if(coreid==0)  
                std::fill(Ta.begin(),Ta.end(),SPComplexType(0));  
              SpRl.barrier();
Timer.reset("Generic2");
Timer.start("Generic2");

#if defined(HAVE_MKL)
              if(coreid < nk*NAEA) {
                // partition on the fly
                int nwork = std::min(nk*NAEA,ncores);
                int ak0,akN;
                std::tie(ak0, akN) = FairDivideBoundary(coreid,nk*NAEA,nwork);

                char TRANS = 'N';
                int M_ = (akN-ak0), N_ = norb*NAEA;
                int ldc_ = nk*NAEA;
                rowI.resize(M_+1);
                SPComplexSMSpMat::intPtr r = SptQk.row_index(ak0);
                int p0 = *r; 
                std::vector<SPComplexSMSpMat::intType>::iterator it=rowI.begin();
                std::vector<SPComplexSMSpMat::intType>::iterator itend=rowI.end();
                int i_=0;
                for(; it!=itend; it++,r++) 
                  *it = *r-p0+1;
                // CAREFUL!!! C matrix in fortran format
                mkl_zcsrmultd (&TRANS, &M_ , &N_ , &nvec , 
                        SptQk.values(p0-1), SptQk.column_data(p0-1) , rowI.data(), //SptQk.row_index(ak0), 
                        SpRl.values()  , SpRl.column_data()  , SpRl.row_index() , 
                        Ta.values()+ak0, &ldc_ );     
              }
              SptQk.barrier();  
#else 
              APP_ABORT(" Error: Requires MKL. Code alternative (Talk to Miguel) \n\n\n");
#endif

Timer.stop("Generic2");
app_log()<<"Time for zcsrmultd : " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;

              if(coreid==0) {
                for(int a=0; a<nk*NAEA; a++)
                  for(int b=0; b<norb*NAEA; b++)
                    if( std::abs(*(Qa.values()+a*norb*NAEA+b) - *(Ta.values()+b*nk*NAEA+a) ) > 5e-5)
                      app_log()<<" Diff: " <<a <<" " <<b <<" " <<*(Qa.values()+a*norb*NAEA+b) <<" " <<*(Ta.values()+b*nk*NAEA+a) <<std::endl;
              }  

              SptQk.toZeroBase();
              SpRl.toZeroBase();
              SpQk.barrier();  
*/
///*
//Timer.reset("Generic2");
//Timer.start("Generic2");
              SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),
                            SptQk.values(),SptQk.column_data(), SptQk.row_index(), 
                            Rl.values()+bl0,norb*NAEA,SPComplexType(0),Ta.values()+bl0,norb*NAEA);
              SpQk.barrier();  
//Timer.stop("Generic2");
//app_log()<<"Time for sparse DGEMM : " <<nn <<"  " <<Timer.total("Generic2") <<std::endl;
//              SpQk.barrier();  

//              if(coreid==0) {
//                for(int a=0; a<nk*NAEA; a++)
//                  for(int b=0; b<norb*NAEA; b++)
//                    if( std::abs(*(Qa.values()+a*norb*NAEA+b) - *(Ta.values()+a*norb*NAEA+b) ) > 1e-6)
//                      app_log()<<" Diff: " <<a <<" " <<b <<" " <<*(Qa.values()+a*norb*NAEA+b) <<" " <<*(Ta.values()+a*norb*NAEA+b) <<std::endl;
//              }

//*/
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEA; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEA+b  
                      int b = lb%NAEA;
                      if(b<=a) continue;
                      int l = lb/NAEA+l0;
                      int la = (l-l0)*NAEA+a;
                      int kb = (k-k0)*NAEA+b;
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEA + lb);  // Qa(ka,lb)   
                      SPComplexType qlakb = *(Ta.values() + kb*norb*NAEA + la);  // Qa(kb,la)
                      if(std::abs( qkalb - qlakb ) > cut) {
                        sz_1++;
                        abkl.push_back( std::make_tuple(a*NMO+k, b*NMO+l, 2.0*(qkalb - qlakb)) );
                        if(abkl.size()==nmax) {
                          Vijkl.add(abkl,true);
                          abkl.clear();
                        }
                      }
                    }
                  }
                }

              } else {
                //DenseMatrixOperators::product(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEB,SPComplexType(0),Qa.values()+bl0,norb*NAEB);
                SparseMatrixOperators::product_SpMatM(nk*NAEA,int(blN-bl0),nvec,SPComplexType(1.0),
                            SptQk.values(),SptQk.column_data(), SptQk.row_index(),
                            Rl.values()+bl0,norb*NAEB,SPComplexType(0),Ta.values()+bl0,norb*NAEB);
                Ta.barrier();
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEA; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEB+b  
                      int b = lb%NAEB;
                      int l = lb/NAEB+l0;
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEB + lb);  // Qa(ka,lb)   
                      if(std::abs( qkalb ) > cut) { 
                        sz_2++;
                        abkl.push_back( std::make_tuple(a*NMO+k, NMO*NMO+b*NMO+l-NMO, 2.0*(qkalb)));
                        if(abkl.size()==nmax) {
                          Vijkl.add(abkl,true);
                          abkl.clear();
                        }
                      }
                    }
                  }
                }
              }  

            } else if(M_split[nn] >= NMO && M_split[nn+1] >= NMO) {
              // k is beta
              if(!amIAlpha) {
                // <a,b | k,l> = Qa(ka,lb) - Qb(la,kb)
                // Qa(ka,lb) = Q(k,a,:)*R(l,b,:)
                //DenseMatrixOperators::product(nk*NAEB,int(blN-bl0),nvec,SPComplexType(1.0),tQk.values(),nvec,Rl.values()+bl0,norb*NAEB,SPComplexType(0),Qa.values()+bl0,norb*NAEB);
                //Qa.barrier(); 
                SparseMatrixOperators::product_SpMatM(nk*NAEB,int(blN-bl0),nvec,SPComplexType(1.0),
                            SptQk.values(),SptQk.column_data(), SptQk.row_index(),
                            Rl.values()+bl0,norb*NAEB,SPComplexType(0),Ta.values()+bl0,norb*NAEB);
                Ta.barrier(); 
                for(int k=k0, ka=0; k<kN; k++) {
                  for(int a=0; a<NAEB; a++, ka++) {
                    for(int lb=bl0; lb<blN; lb++) { // lb = (l-l0)*NAEB+b  
                      int b = lb%NAEB;
                      if(b<=a) continue;
                      int l = lb/NAEB+l0;
                      int la = (l-l0)*NAEB+a;
                      int kb = (k-k0)*NAEB+b;
                      SPComplexType qkalb = *(Ta.values() + ka*norb*NAEB + lb);  // Qa(ka,lb)   
                      SPComplexType qlakb = *(Ta.values() + kb*norb*NAEB + la);  // Qa(kb,la)   
                      if(std::abs( qkalb - qlakb ) > cut) {
                        sz_3++;
                        abkl.push_back( std::make_tuple(NMO*NMO+a*NMO+k-NMO, NMO*NMO+b*NMO+l-NMO, 2.0*(qkalb - qlakb)) );
                        if(abkl.size()==nmax) {
                          Vijkl.add(abkl,true);
                          abkl.clear();
                        }
                      } 
                    }
                  }
                }
              }
            } else {
              APP_ABORT(" Error: This should not happen. \n\n\n");
            }

          } else if(walker_type == 2) {
            APP_ABORT(" Error in createHamiltonianForGeneralDeterminant: GHF not implemented. \n\n\n");
          }
        } // bi < nblk
      }

      if(abkl.size()>0) {
        Vijkl.add(abkl,true);
        abkl.clear();
      }
      SpQk.barrier();
      myComm->barrier();

      //debug debug debug
      assert(sz1==sz_1);
      assert(sz2==sz_2);
      assert(sz3==sz_3);

      Timer.stop("Generic");
      app_log()<<"Time to calculate distributed Hamiltonian: " <<Timer.total("Generic") <<std::endl;

      if(!communicate_Vijkl(Vijkl)) return false;
      myComm->barrier();

      Timer.stop("Generic");
      app_log()<<"Time to generate 2-body Hamiltonian: " <<Timer.total("Generic") <<std::endl;

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
             s1 = fct*ValueType(4.0)*V*(*itAia);
             s2 = fct*ValueType(2.0)*V*(*itAja);
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
               s1 = fct*4.0*V*(*itAia);
               s2 = fct*2.0*V*(*itAja);
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
        APP_ABORT("Error: TODO: createHamiltonianForGeneralDeterminant not implemented for GHF walker type yet. \n\n\n");
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
             s1 = fct*4.0*V*(*itAia);
             s2 = fct*2.0*V*(*itAja);
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
             Vijkl.add( a*NMO+k0, b*NMO+l0, static_cast<SPComplexType>(Qa(a,b)), true);
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
                 Vijkl.add( a*NMO+k0, b*NMO+l0, static_cast<SPComplexType>(Qa(a,b)), true);
               if(std::abs(Qb(a,b)) > cut) 
                 Vijkl.add( NMO*NMO+a*NMO+k0, NMO*NMO+b*NMO+l0, static_cast<SPComplexType>(Qb(a,b)), true);
               if(std::abs(Qab(a,b)) > cut)
                 Vijkl.add( a*NMO+k0, NMO*NMO+b*NMO+l0, static_cast<SPComplexType>(Qab(a,b)), true);
               if(std::abs(Qab2(a,b)) > cut && k0!=l0)
                 Vijkl.add( a*NMO+l0, NMO*NMO+b*NMO+k0, static_cast<SPComplexType>(Qab2(a,b)), true);
             }

          } else {
            app_error()<<" Error: Finish implementation. \n\n\n";
            return false;
            if(k0 < NMO && l0 < NMO) { // aa 
              for(int a=0; a<NAEA; a++, itAia++, itAja++) {
               s1 = fct*4.0*V*(*itAia);
               s2 = fct*2.0*V*(*itAja);
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
        APP_ABORT("Error: TODO: createHamiltonianForGeneralDeterminant not implemented for GHF walker type yet. \n\n\n");
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

    if(!Vijkl.remove_repeated_and_compress(TG.getNodeCommLocal())) {
      APP_ABORT("Error in call to SparseMatrix::remove_repeated(). \n");
    }

#ifdef AFQMC_DEBUG
      app_log()<<" Done compressing sparse hamiltonians. " <<std::endl;
#endif
    
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
            app_error()<<"Error in SparseGeneralHamiltonian :: createHamiltonianForPureDeterminant(). I should never get here. \n"; 
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

