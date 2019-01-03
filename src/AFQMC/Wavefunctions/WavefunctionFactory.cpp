//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#include<random>

#include "io/hdf_archive.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Wavefunctions/NOMSD.hpp"
#include "AFQMC/Wavefunctions/PHMSD.hpp"
#include "AFQMC/HamiltonianOperations/HamOpsIO.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

Wavefunction WavefunctionFactory::fromASCII(TaskGroup_& TGprop, TaskGroup_& TGwfn, 
                                            xmlNodePtr cur, WALKER_TYPES walker_type, Hamiltonian& h, 
                                            RealType cutvn, int targetNW)
{
  if(cur == NULL)
    APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

  std::string type("msd");
  std::string info("info0");
  std::string init_type("");
  std::string name("");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(type,"type");
  oAttrib.add(info,"info");
  oAttrib.add(init_type,"init");
  oAttrib.add(name,"name");
  oAttrib.put(cur);

  std::transform(type.begin(),type.end(),type.begin(),(int (*)(int)) tolower);
  std::transform(init_type.begin(),init_type.end(),init_type.begin(),(int (*)(int)) tolower);

  if(InfoMap.find(info) == InfoMap.end()) {
    app_error()<<"ERROR: Undefined info in WavefunctionFactory. \n";
    APP_ABORT("ERROR: Undefined info in WavefunctionFactory. \n");
  }

  RealType cutv2(0.);   
  int initialDet(1);
  int ndets_to_read(-1); // if not set, read the entire file
  int initial_configuration=0;  
  double randomize_guess(0.0);
  std::string starting_det("");
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  std::string write_trial_density_matrix("");
  ParameterSet m_param;
  m_param.add(filename,"filename","std::string");
  m_param.add(restart_file,"restart_file","std::string");
  m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
  m_param.add(cutv2,"cutoff","double");
  m_param.add(initialDet,"initialDetType","int");
  m_param.add(starting_det,"starting_det","std:string");
  m_param.add(ndets_to_read,"ndet","int");
  m_param.add(initial_configuration,"initial_configuration","int");
  m_param.add(randomize_guess,"randomize_guess","double");
  m_param.put(cur);

  AFQMCInfo& AFinfo = InfoMap[info];
  auto NCE = h.getNuclearCoulombEnergy();

  int NMO = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;

  std::ifstream in;
  in.open(filename.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening file:  " <<filename <<std::endl;
     APP_ABORT("Problems opening file. \n");
  }
  std::string wfn_type = getWfnType(in);

  using Alloc = boost::mpi3::intranode::allocator<ComplexType>;
  if(type=="msd" || type=="nomsd") {

    app_log()<<" Wavefunction type: NOMSD\n";

    std::vector<ComplexType> ci;
    std::vector<PsiT_Matrix> PsiT;
    if(wfn_type == "occ" || wfn_type == "mixed") {
      std::vector<PsiT_Matrix> PsiT_MO; // read_ph_wavefunction returns the full MO matrix
      // careful here!!!!
      // number of terms in PsiT_MO depend on RHF/UHF type, not on walker_type!!!
      ph_excitations<int,ComplexType> abij = read_ph_wavefunction(in,ndets_to_read,walker_type,
                    TGwfn.Node(),NMO,NAEA,NAEB,PsiT_MO);
      assert(abij.number_of_configurations() == ndets_to_read);
      int NEL = (walker_type==NONCOLLINEAR)?(NAEA+NAEB):NAEA;  
      int N_ = (walker_type==NONCOLLINEAR)?2*NMO:NMO;
      ComplexType one(1.0,0.0);
      if(walker_type==COLLINEAR) 
        PsiT.reserve(2*ndets_to_read);  
      else
        PsiT.reserve(ndets_to_read);  
      ci.reserve(ndets_to_read);
      auto refc = abij.reference_configuration();
      // add reference  
      ci.emplace_back(std::get<2>(*abij.configurations_begin()));  
      if(wfn_type == "occ")
        PsiT.emplace_back(PsiT_Matrix({NEL,N_},{0,0},1,Alloc(TGwfn.Node())));
      else 
        PsiT.emplace_back(PsiT_Matrix({NEL,N_},{0,0},
                          get_nnz(PsiT_MO[0],refc,NEL,0),Alloc(TGwfn.Node())));
      if(TGwfn.Node().root()) {
        if(wfn_type == "occ") {
          for(int k=0; k<NEL; k++)
            PsiT.back().emplace_back({k,*(refc+k)},one);
        } else {
          for(int k=0; k<NEL; k++) {
            size_t ki = *(refc+k); // occupied state #k
            auto col = PsiT_MO[0].non_zero_indices2_data(ki);
            auto val = PsiT_MO[0].non_zero_values_data(ki);
            for(size_t ic=0, icend=PsiT_MO[0].num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
              PsiT.back().emplace_back({k,*col},*val);
          }
        }
      }
      if(walker_type==COLLINEAR) {
        if(wfn_type == "occ")
          PsiT.emplace_back(PsiT_Matrix({NAEB,NMO},{0,0},1,Alloc(TGwfn.Node())));
        else
          PsiT.emplace_back(PsiT_Matrix({NAEB,NMO},{0,0},
                            get_nnz(PsiT_MO.back(),refc+NAEA,NAEB,NMO),Alloc(TGwfn.Node())));
        if(TGwfn.Node().root()) {
          if(wfn_type == "occ") {
            for(int k=0; k<NAEB; k++)
              PsiT.back().emplace_back({k,*(refc+NAEA+k)-NMO},one);
          } else {
            for(int k=0; k<NAEB; k++) {
              size_t ki = *(refc+NAEA+k)-NMO; // occupied state #k
              // this should be correct 
              auto col = PsiT_MO.back().non_zero_indices2_data(ki);
              auto val = PsiT_MO.back().non_zero_values_data(ki);
              for(size_t ic=0, icend=PsiT_MO.back().num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
                PsiT.back().emplace_back({k,*col},*val);
            }
          }
        } 
      }
      // work array
      std::vector<int> iwork(NAEA); 
      auto configurations = abij.configurations_begin()+1;
      for(; configurations<abij.configurations_end(); ++configurations) {
        ci.emplace_back(std::get<2>(*configurations));
        abij.get_alpha_configuration( std::get<0>(*configurations) ,iwork); 
        if(wfn_type == "occ")
          PsiT.emplace_back(PsiT_Matrix({NEL,N_},{0,0},1,Alloc(TGwfn.Node())));
        else 
          PsiT.emplace_back(PsiT_Matrix({NEL,N_},{0,0},
                        get_nnz(PsiT_MO[0],iwork.data(),NEL,0),Alloc(TGwfn.Node())));
        if(TGwfn.Node().root()) { 
          // add excited configuration
          if(wfn_type == "occ") {
            for(int k=0; k<NEL; k++) 
              PsiT.back().emplace_back({k,iwork[k]},one); 
          } else {
            for(int k=0; k<NEL; k++) {
              size_t ki = iwork[k]; // occupied state #k
              auto col = PsiT_MO[0].non_zero_indices2_data(ki);
              auto val = PsiT_MO[0].non_zero_values_data(ki);
              for(size_t ic=0, icend=PsiT_MO[0].num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
                PsiT.back().emplace_back({k,*col},*val);
            }
          }
        }
        if(walker_type==COLLINEAR) {
          abij.get_beta_configuration( std::get<1>(*configurations) ,iwork); 
          if(wfn_type == "occ")
            PsiT.emplace_back(PsiT_Matrix({NAEB,NMO},{0,0},1,Alloc(TGwfn.Node())));
          else
            PsiT.emplace_back(PsiT_Matrix({NAEB,NMO},{0,0},
                            get_nnz(PsiT_MO.back(),iwork.data(),NAEB,NMO),Alloc(TGwfn.Node())));
          if(TGwfn.Node().root()) { 
            if(wfn_type == "occ") {
              for(int k=0; k<NAEB; k++) 
                PsiT.back().emplace_back({k,iwork[k]-NMO},one); 
            } else {
              for(int k=0; k<NAEB; k++) {
                size_t ki = iwork[k]-NMO; // occupied state #k
                auto col = PsiT_MO.back().non_zero_indices2_data(ki);
                auto val = PsiT_MO.back().non_zero_values_data(ki);
                for(size_t ic=0, icend=PsiT_MO.back().num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
                  PsiT.back().emplace_back({k,*col},*val);
              }
            }
          } 
        } 
      } 
    } else if(wfn_type == "matrix") {
      read_general_wavefunction(in,ndets_to_read,walker_type,TGwfn.Node(),
                    NMO,NAEA,NAEB,PsiT,ci);
    } else {
      APP_ABORT("Error: Unknown wfn_type in WavefunctionFactory with MSD wavefunction.\n");
    }
    TGwfn.node_barrier();

    // if requested, create restart file
    // Will use phdf5 in the future, for now only head node writes
    hdf_archive dump(TGwfn.Global());
    if(restart_file != "") {
      if(TGwfn.Global().root()) {
        if(!dump.create(restart_file,H5F_ACC_EXCL)) {
          app_error()<<" Error opening restart file in WavefunctionFactory. \n";
          APP_ABORT("");
        }
        dump.push("Wavefunction");
        dump.push("NOMSD");
        std::vector<int> dims{NMO,NAEA,NAEB,walker_type,ndets_to_read};
        dump.write(dims,"dims"); 
        std::vector<ValueType> et{NCE}; 
        dump.write(et,"NCE"); 
        dump.write(ci,"CICOEFFICIENTS"); 
        for(int i=0; i<PsiT.size(); ++i) { 
          dump.push(std::string("PsiT_")+std::to_string(i));
          csr_hdf5::CSR2HDF(dump,PsiT[i]);
          dump.pop();  
        }
      }  
    }
    // multideterminant NOMSD needs 
    auto HOps(h.getHamiltonianOperations(wfn_type == "occ" && walker_type==CLOSED,
			ndets_to_read>1,walker_type,PsiT,cutvn,cutv2,TGprop,TGwfn,dump));  
    if(restart_file != "") 
      if(TGwfn.Global().root()) {
        dump.pop();
        dump.pop();
        dump.close();
      }  
    TGwfn.node_barrier();
    // add initial_guess
    auto guess = initial_guess.find(name);
    if( guess == initial_guess.end() ) {
      auto newg = initial_guess.insert(
                std::make_pair(name,boost::multi::array<ComplexType,3>({2,NMO,NAEA})));
      int iC = (walker_type!=COLLINEAR?initial_configuration:2*initial_configuration);
      if( iC >= PsiT.size() )
        APP_ABORT(" Error: initial_configuration > ndets_to_read \n"); 
      if(!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n"); 
      using std::conj;
      std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
      {  
        auto pbegin = PsiT[iC].pointers_begin(); 
        auto pend = PsiT[iC].pointers_end(); 
        auto p0 = pbegin[0]; 
        auto v0 = PsiT[iC].non_zero_values_data();  
        auto c0 = PsiT[iC].non_zero_indices2_data();  
        for(int i=0; i<PsiT[iC].shape()[0]; i++) 
          for(int ip=pbegin[i]; ip<pend[i]; ip++) { 
            ((newg.first)->second)[0][c0[ip-p0]][i] = conj(v0[ip-p0]);
          }
      }  
      if(walker_type==COLLINEAR) {
        auto pbegin = PsiT[iC+1].pointers_begin(); 
        auto pend = PsiT[iC+1].pointers_end();  
        auto p0 = pbegin[0];                
        auto v0 = PsiT[iC+1].non_zero_values_data();
        auto c0 = PsiT[iC+1].non_zero_indices2_data();
        for(int i=0; i<PsiT[iC+1].shape()[0]; i++) 
          for(int ip=pbegin[i]; ip<pend[i]; ip++) { 
            ((newg.first)->second)[1][c0[ip-p0]][i] = conj(v0[ip-p0]);  
          }
      }  
    } else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n"); 

    //return Wavefunction{}; 
    return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW)); 
  } else if(type == "phmsd") {

    app_log()<<" Wavefunction type: PHMSD\n";

    /* Implementation notes: 
     *  - PsiT: [Nact, NMO] where Nact is the number of active space orbitals, 
     *                     those that participate in the ci expansion
     *  - The half rotation is done with respect to the supermatrix PsiT
     *  - Need to calculate Nact and create a mapping from orbital index to actice space index.
     *    Those orbitals in the corresponding virtual space (not in active) map to -1 as a precaution.
     */   

    // assuming walker_type==COLLINEAR for now, specialize a type for perfect pairing PHMSD
    if(walker_type!=COLLINEAR)
      APP_ABORT("Error: PHMSD requires a COLLINEAR calculation.\n");
    std::vector<PsiT_Matrix> PsiT_MO;
    ph_excitations<int,ComplexType> abij = read_ph_wavefunction(in,ndets_to_read,walker_type,
                  TGwfn.Node(),NMO,NAEA,NAEB,PsiT_MO);
    int NEL = (walker_type==NONCOLLINEAR)?(NAEA+NAEB):NAEA;  
    int N_ = (walker_type==NONCOLLINEAR)?2*NMO:NMO;
    if(wfn_type == "occ") {
      //build PsiT_MO
      ComplexType one(1.0,0.0);
      // wfn_type == "occ" implies a single reference now, since integrals can't be UHF
      PsiT_MO.reserve(1);
      PsiT_MO.emplace_back(PsiT_Matrix({N_,N_},{0,0},1,Alloc(TGwfn.Node())));
      if(TGwfn.Node().root())  
        for(int k=0; k<N_; k++) 
          PsiT_MO.back().emplace_back({k,k},one); 
    } else if(wfn_type == "mixed") {
      // nothing to do
    } else if(wfn_type == "matrix") {
      APP_ABORT("Error: wfn_type=matrix not allowed in WavefunctionFactory with PHMSD wavefunction.\n");
    } else {
      APP_ABORT("Error: Unknown wfn_type in WavefunctionFactory with MSD wavefunction.\n");
    }
    TGwfn.node_barrier();

    // find active space orbitals and create super trial matrix PsiT
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve( PsiT_MO.size() );
    // expect mapped over range [0-2*NMO], but alpha and beta sectors with 0-based active indexes
    std::map<int,int> mo2active(find_active_space(PsiT_MO.size()==1,abij,NMO,NAEA,NAEB));
    std::map<int,int> acta2mo; 
    std::map<int,int> actb2mo; 
    std::vector<int> active_alpha;
    std::vector<int> active_beta;
    std::vector<int> active_combined;
    for(int i=0; i<NMO; i++) {
      if(mo2active[i]>=0) {
        active_alpha.push_back(i);
        acta2mo[mo2active[i]] = i;
      }  
      if(mo2active[i+NMO]>=0) {
        active_beta.push_back(i);
        actb2mo[mo2active[i+NMO]] = i+NMO;
      }  
      if(mo2active[i]>=0 || mo2active[i+NMO]>=0) active_combined.push_back(i);
    }
    if(PsiT_MO.size() == 1) {
      // RHF reference
      PsiT.emplace_back(PsiT_Matrix({active_combined.size(),NMO},{0,0},
                     get_nnz(PsiT_MO[0],active_combined.data(),active_combined.size(),0),
                     Alloc(TGwfn.Node())));
      if(TGwfn.Node().root()) {
        for(int k=0; k<active_combined.size(); k++) {
          size_t ki = active_combined[k]; // occupied state #k
          auto col = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val = PsiT_MO[0].non_zero_values_data(ki);
          for(size_t ic=0, icend=PsiT_MO[0].num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
            PsiT[0].emplace_back({k,*col},*val);
        }
      }
    } else {
      // UHF reference
      PsiT.emplace_back(PsiT_Matrix({active_alpha.size(),NMO},{0,0},
                     get_nnz(PsiT_MO[0],active_alpha.data(),active_alpha.size(),0),
                     Alloc(TGwfn.Node())));
      if(TGwfn.Node().root()) {
        for(int k=0; k<active_alpha.size(); k++) {
          size_t ki = active_alpha[k]; // occupied state #k
          auto col = PsiT_MO[0].non_zero_indices2_data(ki);
          auto val = PsiT_MO[0].non_zero_values_data(ki);
          for(size_t ic=0, icend=PsiT_MO[0].num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
            PsiT.back().emplace_back({k,*col},*val);
        }
      }
      PsiT.emplace_back(PsiT_Matrix({active_beta.size(),NMO},{0,0},
                     get_nnz(PsiT_MO[1],active_beta.data(),active_beta.size(),0),
                     Alloc(TGwfn.Node())));
      if(TGwfn.Node().root()) {
        for(int k=0; k<active_beta.size(); k++) {
          size_t ki = active_beta[k]; // occupied state #k
          auto col = PsiT_MO[1].non_zero_indices2_data(ki);
          auto val = PsiT_MO[1].non_zero_values_data(ki);
          for(size_t ic=0, icend=PsiT_MO[1].num_non_zero_elements(ki); ic<icend; ic++, ++col, ++val)
            PsiT[1].emplace_back({k,*col},*val);
        }
      }
    }
    // now that mappings have been constructed, map indexes of excited state orbitals
    // to the corresponding active space indexes
    if(TGwfn.Node().root()) {
      // map reference  
      auto refc = abij.reference_configuration();
      for(int i=0; i<NAEA+NAEB; i++, ++refc) *refc = mo2active[*refc]; 
      for(int n=1; n<abij.maximum_excitation_number()[0]; n++) {  
        auto it = abij.alpha_begin(n);
        auto ite = abij.alpha_end(n);
        for(; it<ite; ++it) {
          auto exct = (*it)+n; // only need to map excited state indexes
          for(int np=0; np<n; ++np, ++exct)
            *exct = mo2active[*exct];  
        }       
      }
      for(int n=1; n<abij.maximum_excitation_number()[1]; n++) {
        auto it = abij.beta_begin(n);
        auto ite = abij.beta_end(n);
        for(; it<ite; ++it) {
          auto exct = (*it)+n; // only need to map excited state indexes
          for(int np=0; np<n; ++np, ++exct)
            *exct = mo2active[*exct]; 
        } 
      }
    }
    TGwfn.node_barrier();
    

    // if requested, create restart file
    // Will use phdf5 in the future, for now only head node writes
    hdf_archive dump(TGwfn.Global());
/* no restart yet!
    if(restart_file != "") {
      if(TGwfn.Global().root()) {
        if(!dump.create(restart_file,H5F_ACC_EXCL)) {
          app_error()<<" Error opening restart file in WavefunctionFactory. \n";
          APP_ABORT("");
        }
        dump.push("Wavefunction");
        dump.push("NOMSD");
        std::vector<int> dims{NMO,NAEA,NAEB,walker_type,ndets_to_read};
        dump.write(dims,"dims"); 
        std::vector<ValueType> et{NCE}; 
        dump.write(et,"NCE"); 
        dump.write(ci,"CICOEFFICIENTS"); 
        for(int i=0; i<PsiT.size(); ++i) { 
          dump.push(std::string("PsiT_")+std::to_string(i));
          csr_hdf5::CSR2HDF(dump,PsiT[i]);
          dump.pop();  
        }
      }  
    }
*/
    WALKER_TYPES reference_type = (PsiT.size()==1?CLOSED:COLLINEAR);
    // is PureSD actually faster??? CHECK!!!
    // NOTE: For UHF reference, treat HOps as a 2 determinant wavefunction of a
    //        CLOSED walker type. This way you can treat alpha/beta sectors independently in 
    //        HOps through the index corresponding 0/1.			
    // never add coulomb to half rotated v2 tensor in PHMSD	
    //auto HOps(h.getHamiltonianOperations(wfn_type == "occ",false, 
    auto HOps(h.getHamiltonianOperations(false,false,
                                         CLOSED,PsiT,cutvn,cutv2,TGprop,TGwfn,dump));  
/*
    if(restart_file != "") 
      if(TGwfn.Global().root()) {
        dump.pop();
        dump.pop();
        dump.close();
      }  
*/
    TGwfn.node_barrier();
    // add initial_guess
    // when propagating Nact states, change this here
    auto guess = initial_guess.find(name);
    if( guess == initial_guess.end() ) {
      auto newg = initial_guess.insert(
                std::make_pair(name,boost::multi::array<ComplexType,3>({2,NMO,NAEA})));
      if(!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n"); 
      auto& Psi0((newg.first)->second);
      randomize_guess = std::abs(randomize_guess);
      if(randomize_guess > 1e-12) 
        app_log()<<" Randomizing initial guess with uniform distribution: " <<randomize_guess <<std::endl;
      std::default_random_engine generator(777);
      std::uniform_real_distribution<double> distribution(-randomize_guess,randomize_guess);  
      int iC = initial_configuration;
      if( iC >= abij.number_of_configurations() )
        APP_ABORT(" Error: initial_configuration > ndets \n");
      using std::conj;
      std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
      //auto refc = abij.reference_configuration();
      {  
        std::vector<int> alphaC(NAEA);
        abij.get_alpha_configuration(std::get<0>(*(abij.configurations_begin()+iC)),alphaC);
        auto pbegin = PsiT[0].pointers_begin(); 
        auto pend = PsiT[0].pointers_end(); 
        auto p0 = pbegin[0]; 
        auto v0 = PsiT[0].non_zero_values_data();  
        auto c0 = PsiT[0].non_zero_indices2_data();  
        // only takinf NAEA states, increase later if super SM is needed
        for(int i=0; i<NAEA; i++) { 
          //int ik = *(refc+i);
          int ik = alphaC[i]; 
          for(int ip=pbegin[ik]; ip<pend[ik]; ip++) 
            Psi0[0][c0[ip-p0]][i] = conj(v0[ip-p0]);
        }
        if(randomize_guess > 1e-12)
          for(int i=0; i<NMO; i++)
            for(int a=0; a<NAEA; a++)
              Psi0[0][i][a] += distribution(generator);
      }  
      if(walker_type==COLLINEAR) {
        std::vector<int> betaC(NAEB);
        abij.get_beta_configuration(std::get<1>(*(abij.configurations_begin()+iC)),betaC);
        auto pbegin = PsiT.back().pointers_begin(); 
        auto pend = PsiT.back().pointers_end();  
        auto p0 = pbegin[0];                
        auto v0 = PsiT.back().non_zero_values_data();
        auto c0 = PsiT.back().non_zero_indices2_data();
        // only takinf NAEB states, increase later if super SM is needed
        for(int i=0; i<NAEB; i++) { 
          //int ik = *(refc+NAEA+i);
          int ik = betaC[i]; 
          for(int ip=pbegin[ik]; ip<pend[ik]; ip++) 
            Psi0[1][c0[ip-p0]][i] = conj(v0[ip-p0]);  
        }
        if(randomize_guess > 1e-12)
          for(int i=0; i<NMO; i++)
            for(int a=0; a<NAEB; a++)
              Psi0[1][i][a] += distribution(generator);
      }  
    } else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n"); 

    // setup configuration coupligs
    using index_aos = ma::sparse::array_of_sequences<int,int,
                                                   boost::mpi3::intranode::allocator<int>,
                                                   boost::mpi3::intranode::is_root>;
//    std::allocator<ComplexType> alloc_{}; //boost::mpi3::intranode::allocator<ComplexType>;
    boost::mpi3::intranode::allocator<int> alloc_{TGwfn.Node()};
    
    // alpha
    std::vector<int> counts_alpha(abij.number_of_unique_excitations()[0]);
    std::vector<int> counts_beta(abij.number_of_unique_excitations()[1]);
    if(TGwfn.Node().root()) {
      for(auto it=abij.configurations_begin(); it<abij.configurations_end(); ++it) {
        ++counts_alpha[std::get<0>(*it)];
        ++counts_beta[std::get<1>(*it)];
      }
    }
    TGwfn.Node().broadcast_n(counts_alpha.begin(),counts_alpha.size());
    TGwfn.Node().broadcast_n(counts_beta.begin(),counts_beta.size());
    index_aos beta_coupled_to_unique_alpha(counts_alpha.size(),counts_alpha,alloc_);  
    index_aos alpha_coupled_to_unique_beta(counts_beta.size(),counts_beta,alloc_);  
    if(TGwfn.Node().root()) {
      int ni=0;
      for(auto it=abij.configurations_begin(); it<abij.configurations_end(); ++it, ++ni) {
        beta_coupled_to_unique_alpha.emplace_back(std::get<0>(*it),ni);
        alpha_coupled_to_unique_beta.emplace_back(std::get<1>(*it),ni);
      }
    }
    TGwfn.Node().barrier();

    //return Wavefunction{}; 
    return Wavefunction(PHMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(acta2mo),
                        std::move(actb2mo),std::move(abij),std::move(beta_coupled_to_unique_alpha),
                        std::move(alpha_coupled_to_unique_beta),std::move(PsiT),
                        walker_type,NCE,targetNW)); 
  } else if(type == "generalmsd") {
    app_error()<<" Error: Wavefunction type GeneralMSD not yet implemented. \n"; 
    APP_ABORT(" Error: Wavefunction type GeneralMSD not yet implemented. \n");
    return Wavefunction{}; 
  } else {
    app_error()<<" Error: Unknown wave-function type: " <<type <<std::endl;
    APP_ABORT(" Error: Unknown wave-function type. \n");
    return Wavefunction{}; 
  }
}

Wavefunction WavefunctionFactory::fromHDF5(TaskGroup_& TGprop, TaskGroup_& TGwfn, 
                                           xmlNodePtr cur, WALKER_TYPES walker_type,  
                                            RealType cutvn, int targetNW)
{
  if(cur == NULL)
    APP_ABORT("Error: NULL xml pointer in HamiltonianFactory::parse(). \n");

  std::string type("MSD");
  std::string info("info0");
  std::string init_type("");
  std::string name("");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(type,"type");
  oAttrib.add(info,"info");
  oAttrib.add(init_type,"init");
  oAttrib.add(name,"name");
  oAttrib.put(cur);

  std::transform(type.begin(),type.end(),type.begin(),(int (*)(int)) tolower);
  std::transform(init_type.begin(),init_type.end(),init_type.begin(),(int (*)(int)) tolower);

  if(InfoMap.find(info) == InfoMap.end()) {
    app_error()<<"ERROR: Undefined info in WavefunctionFactory. \n";
    APP_ABORT("ERROR: Undefined info in WavefunctionFactory. \n");
  }

  RealType cutv2(0.);
  int initialDet(1);
  int ndets_to_read(-1); // if not set, read the entire file
  std::string starting_det("");
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  std::string write_trial_density_matrix("");
  ParameterSet m_param;
  m_param.add(filename,"filename","std::string");
  m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
  m_param.add(cutv2,"cutoff","double");
  m_param.add(initialDet,"initialDetType","int");
  m_param.add(starting_det,"starting_det","std:string");
  m_param.add(ndets_to_read,"ndet","int");
  m_param.put(cur);

  AFQMCInfo& AFinfo = InfoMap[info];
  ValueType NCE;

  int NMO = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;

  std::vector<ComplexType> ci;
  std::vector<PsiT_Matrix> PsiT;
  std::vector<int> excitations;

  using Alloc = boost::mpi3::intranode::allocator<ComplexType>;
  // HOps, ci, PsiT, NCE
  hdf_archive dump(TGwfn.Global());
  if(!dump.open(filename,H5F_ACC_RDONLY)) {
    app_error()<<" Error hdf5 file in WavefunctionFactory. \n";
    APP_ABORT("");
  }
  if(!dump.push("Wavefunction",false)) {
    app_error()<<" Error in WavefunctionFactory: Group Wavefunction not found. \n";
    APP_ABORT("");
  }
  if(!dump.push("NOMSD",false)) {
    app_error()<<" Error in WavefunctionFactory: Group NOMSD not found. \n";
    APP_ABORT("");
  }

    // check for consistency in parameters
  std::vector<int> dims(5); //{NMO,NAEA,NAEB,walker_type,ndets_to_read};
  if(TGwfn.Global().root()) {
    if(!dump.read(dims,"dims")) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Problems reading dims. \n";
      APP_ABORT("");
    }
    if(NMO!=dims[0]) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Inconsistent NMO . \n";
      APP_ABORT("");
    }
    if(NAEA!=dims[1]) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Inconsistent  NAEA. \n";
      APP_ABORT("");
    }
    if(NAEB!=dims[2]) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Inconsistent  NAEB. \n";
      APP_ABORT("");
    }
    if(walker_type!=dims[3]) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Inconsistent  walker_type. \n";
      APP_ABORT("");
    }
    if(ndets_to_read < 1) ndets_to_read=dims[4];  
    if(ndets_to_read > dims[4]) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Inconsistent  ndets_to_read. \n";
      APP_ABORT("");
    }
    if(!dump.read(ci,"CICOEFFICIENTS")) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Problems reading CICOEFFICIENTS. \n";
      APP_ABORT("");
    }
    ci.resize(ndets_to_read);
    std::vector<ValueType> dum;  
    if(!dump.read(dum,"NCE")) {
      app_error()<<" Error in WavefunctionFactory::fromHDF5(): Problems reading NCE. \n";
      APP_ABORT("");
    }
    NCE = dum[0];  
  }
  TGwfn.Global().broadcast_n(dims.data(),dims.size());
  if(ndets_to_read < 1) ndets_to_read=dims[4];  
  ci.resize(ndets_to_read);
  TGwfn.Global().broadcast_n(ci.data(),ci.size());
  TGwfn.Global().broadcast_value(NCE);

  if(type=="msd" || type=="nomsd") {

    int nd = (walker_type==COLLINEAR?2*ndets_to_read:ndets_to_read);
    PsiT.reserve(nd); 
    using Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    for(int i=0; i<nd; ++i) {
      if(!dump.push(std::string("PsiT_")+std::to_string(i),false)) {
        app_error()<<" Error in WavefunctionFactory: Group PsiT not found. \n";
        APP_ABORT("");
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix,Alloc>(dump,TGwfn.Node())); //,Alloc(TGwfn.Node())));
      dump.pop();
    }

    // add initial_guess
    auto guess = initial_guess.find(name);
    if( guess == initial_guess.end() ) {
      auto newg = initial_guess.insert(
                std::make_pair(name,boost::multi::array<ComplexType,3>({2,NMO,NAEA})));
      if(!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n");
      using std::conj;
      std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
      {
        auto pbegin = PsiT[0].pointers_begin();
        auto pend = PsiT[0].pointers_end();
        auto p0 = pbegin[0];
        auto v0 = PsiT[0].non_zero_values_data();
        auto c0 = PsiT[0].non_zero_indices2_data();
        for(int i=0; i<PsiT[0].shape()[0]; i++)
          for(int ip=pbegin[i]; ip<pend[i]; ip++)
            ((newg.first)->second)[0][c0[ip-p0]][i] = conj(v0[ip-p0]);
      }
      if(walker_type==COLLINEAR) {
        auto pbegin = PsiT[1].pointers_begin();
        auto pend = PsiT[1].pointers_end();
        auto p0 = pbegin[0];
        auto v0 = PsiT[1].non_zero_values_data();
        auto c0 = PsiT[1].non_zero_indices2_data();
        for(int i=0; i<PsiT[1].shape()[0]; i++)
          for(int ip=pbegin[i]; ip<pend[i]; ip++)
            ((newg.first)->second)[1][c0[ip-p0]][i] = conj(v0[ip-p0]);
      }
    } else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");

    HamiltonianOperations HOps(loadHamOps(dump,walker_type,NMO,NAEA,NAEB,PsiT,TGprop,TGwfn,cutvn,cutv2));

    return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW));
  } else if(type == "phmsd") {
/*
    int nd = (walker_type==COLLINEAR?2*ndets_to_read:ndets_to_read);
    PsiT.reserve(nd); 
    using Alloc = boost::mpi3::intranode::allocator<ComplexType>;
    for(int i=0; i<nd; ++i) {
      if(!dump.push(std::string("PsiT_")+std::to_string(i),false)) {
        app_error()<<" Error in WavefunctionFactory: Group PsiT not found. \n";
        APP_ABORT("");
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix,Alloc>(dump,TGwfn.Node())); //,Alloc(TGwfn.Node())));
      dump.pop();
    }

    // add initial_guess
    auto guess = initial_guess.find(name);
    if( guess == initial_guess.end() ) {
      auto newg = initial_guess.insert(
                std::make_pair(name,boost::multi::array<ComplexType,3>({2,NMO,NAEA})));
      if(!newg.second)
        APP_ABORT(" Error: Problems adding new initial guess. \n");
      using std::conj;
      std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
      {
        auto pbegin = PsiT[0].pointers_begin();
        auto pend = PsiT[0].pointers_end();
        auto p0 = pbegin[0];
        auto v0 = PsiT[0].non_zero_values_data();
        auto c0 = PsiT[0].non_zero_indices2_data();
        for(int i=0; i<PsiT[0].shape()[0]; i++)
          for(int ip=pbegin[i]; ip<pend[i]; ip++)
            ((newg.first)->second)[0][c0[ip-p0]][i] = conj(v0[ip-p0]);
      }
      if(walker_type==COLLINEAR) {
        auto pbegin = PsiT[1].pointers_begin();
        auto pend = PsiT[1].pointers_end();
        auto p0 = pbegin[0];
        auto v0 = PsiT[1].non_zero_values_data();
        auto c0 = PsiT[1].non_zero_indices2_data();
        for(int i=0; i<PsiT[1].shape()[0]; i++)
          for(int ip=pbegin[i]; ip<pend[i]; ip++)
            ((newg.first)->second)[1][c0[ip-p0]][i] = conj(v0[ip-p0]);
      }
    } else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");

    HamiltonianOperations HOps(loadHamOps(dump,false,walker_type,NMO,NAEA,NAEB,PsiT,TGprop,TGwfn,cutvn,cutv2));
    return Wavefunction(PHMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW));
*/
    app_error()<<" Error: Wavefunction type PHMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type PHMSD not yet implemented. \n");
    return Wavefunction();
  } else if(type == "generalmsd") {
    app_error()<<" Error: Wavefunction type GeneralMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type GeneralMSD not yet implemented. \n");
    return Wavefunction{};
  } else {
    app_error()<<" Error: Unknown wave-function type: " <<type <<std::endl;
    APP_ABORT(" Error: Unknown wave-function type. \n");
    return Wavefunction{};
  }
}

}

}
