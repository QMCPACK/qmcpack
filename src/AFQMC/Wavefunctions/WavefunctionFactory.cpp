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
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
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
  int ndets_to_read(-1); // if not set, read the entire file
  int initial_configuration=0;  
  double randomize_guess(0.0);
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  std::string write_trial_density_matrix("");
  ParameterSet m_param;
  m_param.add(filename,"filename","std::string");
  m_param.add(restart_file,"restart_file","std::string");
  m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
  m_param.add(cutv2,"cutoff","double");
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

  using Alloc = shared_allocator<ComplexType>;
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
        PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL,N_},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
      else
        PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL,N_},tp_ul_ul{0,0},
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
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB,NMO},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
        else
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB,NMO},tp_ul_ul{0,0},
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
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL,N_},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
        else 
          PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NEL,N_},tp_ul_ul{0,0},
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
            PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB,NMO},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
          else
            PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{NAEB,NMO},tp_ul_ul{0,0},
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
        dump.write(ci,"ci_coeffs");
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
      using ma::conj;
      std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
      {  
        auto pbegin = PsiT[iC].pointers_begin(); 
        auto pend = PsiT[iC].pointers_end(); 
        auto p0 = pbegin[0]; 
        auto v0 = PsiT[iC].non_zero_values_data();  
        auto c0 = PsiT[iC].non_zero_indices2_data();  
        for(int i=0; i<PsiT[iC].size(0); i++) 
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
        for(int i=0; i<PsiT[iC+1].size(0); i++) 
          for(int ip=pbegin[i]; ip<pend[i]; ip++) { 
            ((newg.first)->second)[1][c0[ip-p0]][i] = conj(v0[ip-p0]);  
          }
      }  
    } else
      APP_ABORT(" Error: Problems adding new initial guess, already exists. \n"); 

    if(TGwfn.TG_local().size() > 1) {
      SlaterDetOperations SDetOp( SlaterDetOperations_shared<ComplexType>(
                        ((walker_type!=NONCOLLINEAR)?(NMO):(2*NMO)),
                        ((walker_type!=NONCOLLINEAR)?(NAEA):(NAEA+NAEB)) ));
      return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(SDetOp),std::move(HOps),
                        std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW)); 
    } else 
    {
      SlaterDetOperations SDetOp( SlaterDetOperations_serial<device_allocator<ComplexType>>(
                        ((walker_type!=NONCOLLINEAR)?(NMO):(2*NMO)),
                        ((walker_type!=NONCOLLINEAR)?(NAEA):(NAEA+NAEB)) ));
      return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(SDetOp),std::move(HOps),
                        std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW));
    }

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
      PsiT_MO.emplace_back(PsiT_Matrix(tp_ul_ul{N_,N_},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_combined.size(),NMO},tp_ul_ul{0,0},
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_alpha.size(),NMO},tp_ul_ul{0,0},
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_beta.size(),NMO},tp_ul_ul{0,0},
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
      using ma::conj;
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
                                                   shared_allocator<int>,
                                                   ma::sparse::is_root>;
//    std::allocator<ComplexType> alloc_{}; //shared_allocator<ComplexType>;
    shared_allocator<int> alloc_{TGwfn.Node()};
    
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
                                           Hamiltonian& h, RealType cutvn, int targetNW)
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
  int ndets_to_read(-1); // if not set, read the entire file
  std::string str("false");
  std::string filename("");
  std::string restart_file("");
  ParameterSet m_param;
  m_param.add(filename,"filename","std::string");
  m_param.add(restart_file,"restart_file","std::string");
  m_param.add(cutv2,"cutoff","double");
  m_param.add(ndets_to_read,"ndet","int");
  m_param.put(cur);

  AFQMCInfo& AFinfo = InfoMap[info];
  ValueType NCE;

  int NMO = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;

  std::vector<int> excitations;

  using Alloc = shared_allocator<ComplexType>;
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
  hdf_archive restart(TGwfn.Global());
  if(restart_file != "") {
    if(type == "phmsd") {
      app_log() << " Restarting from PHMSD trial wavefunction not implemented. \n";
    } else {
      if(!restart.create(restart_file,H5F_ACC_EXCL)) {
        app_error()<<" Error opening restart file in WavefunctionFactory. \n";
        APP_ABORT("");
      }
    }
  }

  if(type=="msd" || type=="nomsd") {

    app_log()<<" Wavefunction type: NOMSD\n";
    if(!dump.push("NOMSD",false)) {
      app_error()<<" Error in WavefunctionFactory: Group NOMSD not found.\n";
      APP_ABORT("");
    }
    std::vector<ComplexType> ci;
    ValueType NCE;

    // Read common trial wavefunction input options.
    getCommonInput(dump, NMO, NAEA, NAEB, ndets_to_read, ci,
                   walker_type, TGwfn.Global().root());
    if(restart_file == "") {
      NCE = h.getNuclearCoulombEnergy();
    } else {
      std::vector<ValueType> tmp;
      if(!dump.readEntry(tmp, "NCE")) {
        app_error()<<" Error in WavefunctionFactory::fromHDF5(): Problems reading NCE.\n";
        APP_ABORT("");
      }
      NCE = tmp[0];
    }

    TGwfn.Global().broadcast_n(ci.data(),ci.size());
    TGwfn.Global().broadcast_value(NCE);

    // Create Trial wavefunction.
    int nd = (walker_type==COLLINEAR?2*ndets_to_read:ndets_to_read);
    std::vector<PsiT_Matrix> PsiT;
    PsiT.reserve(nd);
    using Alloc = shared_allocator<ComplexType>;
    for(int i=0; i<nd; ++i) {
      if(!dump.push(std::string("PsiT_")+std::to_string(i),false)) {
        app_error()<<" Error in WavefunctionFactory: Group PsiT not found. \n";
        APP_ABORT("");
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix,Alloc>(dump,TGwfn.Node())); //,Alloc(TGwfn.Node())));
      dump.pop();
    }

    // Set initial walker's Slater matrix.
    getInitialGuess(dump, name, NMO, NAEA, NAEB, walker_type);

    // Restore Hamiltonian Operations Object from file if it exists.
    bool read_ham_op = restart.is_group("HamiltonianOperations");
    auto HOps(getHamOps(read_ham_op,restart,walker_type,NMO,NAEA,NAEB,PsiT,TGprop,TGwfn,cutvn,cutv2,ndets_to_read,h));

    if(TGwfn.TG_local().size() > 1) {
      SlaterDetOperations SDetOp( SlaterDetOperations_shared<ComplexType>(
                        ((walker_type!=NONCOLLINEAR)?(NMO):(2*NMO)),
                        ((walker_type!=NONCOLLINEAR)?(NAEA):(NAEA+NAEB)) ));
      return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(SDetOp),std::move(HOps),
                        std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW));
    } else
    {
      SlaterDetOperations SDetOp( SlaterDetOperations_serial<device_allocator<ComplexType>>(
                        ((walker_type!=NONCOLLINEAR)?(NMO):(2*NMO)),
                        ((walker_type!=NONCOLLINEAR)?(NAEA):(NAEA+NAEB)) ));
      return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(SDetOp),std::move(HOps),
                        std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW));
    }

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
    std::string wfn_type;
    if(!dump.push("PHMSD",false)) {
      app_error()<<" Error in WavefunctionFactory: Group PHMSD not found. \n";
      APP_ABORT("");
    }
    ph_excitations<int,ComplexType> abij = read_ph_wavefunction_hdf(dump,ndets_to_read,walker_type,
                                                                    TGwfn.Node(),NMO,NAEA,NAEB,PsiT_MO,
                                                                    wfn_type);
    int NEL = (walker_type==NONCOLLINEAR)?(NAEA+NAEB):NAEA;
    int N_ = (walker_type==NONCOLLINEAR)?2*NMO:NMO;
    if(wfn_type == "occ") {
      //build PsiT_MO
      ComplexType one(1.0,0.0);
      // wfn_type == "occ" implies a single reference now, since integrals can't be UHF
      PsiT_MO.reserve(1);
      PsiT_MO.emplace_back(PsiT_Matrix(tp_ul_ul{N_,N_},tp_ul_ul{0,0},1,Alloc(TGwfn.Node())));
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_combined.size(),NMO},tp_ul_ul{0,0},
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_alpha.size(),NMO},tp_ul_ul{0,0},
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
      PsiT.emplace_back(PsiT_Matrix(tp_ul_ul{active_beta.size(),NMO},tp_ul_ul{0,0},
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

    WALKER_TYPES reference_type = (PsiT.size()==1?CLOSED:COLLINEAR);
    // is PureSD actually faster??? CHECK!!!
    // NOTE: For UHF reference, treat HOps as a 2 determinant wavefunction of a
    //        CLOSED walker type. This way you can treat alpha/beta sectors independently in
    //        HOps through the index corresponding 0/1.
    // never add coulomb to half rotated v2 tensor in PHMSD
    getInitialGuess(dump, name, NMO, NAEA, NAEB, walker_type);
    auto HOps(h.getHamiltonianOperations(false,false,
                                         CLOSED,PsiT,cutvn,cutv2,TGprop,TGwfn,dump));
    TGwfn.node_barrier();
    // setup configuration couplings
    using index_aos = ma::sparse::array_of_sequences<int,int,shared_allocator<int>,ma::sparse::is_root>;
    shared_allocator<int> alloc_{TGwfn.Node()};

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

    return Wavefunction(PHMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(acta2mo),
                        std::move(actb2mo),std::move(abij),std::move(beta_coupled_to_unique_alpha),
                        std::move(alpha_coupled_to_unique_beta),std::move(PsiT),
                        walker_type,NCE,targetNW));
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

ph_excitations<int,ComplexType> WavefunctionFactory::read_ph_wavefunction_hdf(hdf_archive& dump, int& ndets, WALKER_TYPES walker_type,
        boost::mpi3::shared_communicator& comm, int NMO, int NAEA, int NAEB,
        std::vector<PsiT_Matrix>& PsiT, std::string& type)
{

  using Alloc = shared_allocator<ComplexType>;
  assert(walker_type!=UNDEFINED_WALKER_TYPE);
  bool fullMOMat = false;
  bool Cstyle = true;
  int wfn_type = 0;
  int NEL = NAEA;
  bool mixed=false;
  if(walker_type!=CLOSED) NEL+=NAEB;

  /*
   * Expected order of inputs and tags:
   * Reference:
   * Configurations:
   */

  /*
   * type:
   *   - occ: All determinants are specified with occupation numbers
   *
   * wfn_type:
   *   - 0: excitations out of a RHF reference
   *          NOTE: Does not mean perfect pairing, means excitations from a single reference
   *   - 1: excitations out of a UHF reference
   */
  std::vector<ComplexType> ci_coeff;
  getCommonInput(dump, NMO, NAEA, NAEB, ndets, ci_coeff, walker_type, comm.root());
  if(walker_type != COLLINEAR)
    APP_ABORT(" Error: walker_type!=COLLINEAR not yet implemented in read_ph_wavefunction.\n");

  if(!dump.readEntry(type,"type")) {
    app_error()<<" Error in WavefunctionFactory::fromHDF5(): Problems reading type. \n";
    APP_ABORT("");
  }
  if(!dump.readEntry(fullMOMat,"fullmo")) {
    APP_ABORT("Problems reading fullmo in read_ph_wavefunction_hdf.\n");
  }
  if(type == "mixed") mixed = true;

  if(mixed) { // read reference
    int nmo_ = (walker_type==NONCOLLINEAR?2*NMO:NMO);
    if(not comm.root()) nmo_=0; // only root reads matrices
    if(not fullMOMat)
      APP_ABORT("Error: Wavefunction type mixed requires fullMOMat=true.\n");
    PsiT.reserve((wfn_type!=1)?1:2);

    if(!dump.push(std::string("PsiT_")+std::to_string(0),false)) {
      app_error()<<" Error in WavefunctionFactory: Group PsiT not found. \n";
      APP_ABORT("");
    }
    PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix,Alloc>(dump,comm));
    dump.pop();
    if(wfn_type == 1) {
      if(!dump.push(std::string("PsiT_")+std::to_string(1),false)) {
        app_error()<<" Error in WavefunctionFactory: Group PsiT not found. \n";
        APP_ABORT("");
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix,Alloc>(dump,comm));
      dump.pop();
    }
  }

  ComplexType ci;
  // count number of k-particle excitations
  // counts[0] has special meaning, it must be equal to NAEA+NAEB.
  std::vector<size_t> counts_alpha(NAEA+1);
  std::vector<size_t> counts_beta(NAEB+1);
  // ugly but need dynamic memory allocation
  std::vector<std::vector<int>> unique_alpha(NAEA+1);
  std::vector<std::vector<int>> unique_beta(NAEB+1);
  // reference configuration, taken as the first one right now
  std::vector<int> refa;
  std::vector<int> refb;
  // space to read configurations
  std::vector<int> confg;
  // space for excitation string identifying the current configuration
  std::vector<int> exct;
  // record file position to come back
  std::vector<int> Iwork; // work arrays for permutation calculation
  std::streampos start;
  std::vector<int> buff(ndets*(NAEA+NAEB));
  boost::multi::array_ref<int,2> occs(to_address(buff.data()), {ndets, NAEA+NAEB});
  if(comm.root()) {
    if(!dump.readEntry(buff, "occs"))
      APP_ABORT("Error reading occs array.\n");
    confg.reserve(NAEA);
    Iwork.resize(2*NAEA);
    exct.reserve(2*NAEA);
    for(int i=0; i<ndets; i++) {
      ci = ci_coeff[i];
      // alpha
      confg.clear();
      for(int k=0, q=0; k<NAEA; k++) {
        q = occs[i][k];
        if(q < 0 || q >= NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      if(i==0) {
        refa=confg;
      } else {
        int np = get_excitation_number(true,refa,confg,exct,ci,Iwork);
        push_excitation(exct,unique_alpha[np]);
      }
      // beta
      confg.clear();
      for(int k=0, q=0; k<NAEB; k++) {
        q = occs[i][NAEA+k];
        if(q < NMO || q >= 2*NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      if(i==0) {
        refb=confg;
      } else {
        int np = get_excitation_number(true,refb,confg,exct,ci,Iwork);
        push_excitation(exct,unique_beta[np]);
      }
    }
    // now that we have all unique configurations, count
    for(int i=1; i<=NAEA; i++) counts_alpha[i] = unique_alpha[i].size();
    for(int i=1; i<=NAEB; i++) counts_beta[i] = unique_beta[i].size();
  }
  comm.broadcast_n(counts_alpha.begin(),counts_alpha.size());
  comm.broadcast_n(counts_beta.begin(),counts_beta.size());
  // using int for now, but should move to short later when everything works well
  // ph_struct stores the reference configuration on the index [0]
  ph_excitations<int,ComplexType> ph_struct(ndets,NAEA,NAEB,counts_alpha,counts_beta,
                                            shared_allocator<int>(comm));

  if(comm.root()) {
    std::map<int,int> refa2loc;
    for(int i=0; i<NAEA; i++) refa2loc[refa[i]] = i;
    std::map<int,int> refb2loc;
    for(int i=0; i<NAEB; i++) refb2loc[refb[i]] = i;
    // add reference
    ph_struct.add_reference(refa,refb);
    // add unique configurations
    // alpha
    for(int n=1; n<unique_alpha.size(); n++)
      for(std::vector<int>::iterator it=unique_alpha[n].begin(); it<unique_alpha[n].end(); it+=(2*n) )
        ph_struct.add_alpha(n,it);
    // beta
    for(int n=1; n<unique_beta.size(); n++)
      for(std::vector<int>::iterator it=unique_beta[n].begin(); it<unique_beta[n].end(); it+=(2*n) )
        ph_struct.add_beta(n,it);
    // read configurations
    int alpha_index;
    int beta_index;
    int np;
    for(int i=0; i<ndets; i++) {
      ci = ci_coeff[i];
      confg.clear();
      for(int k=0, q=0; k<NAEA; k++) {
        q = occs[i][k];
        if(q < 0 || q >= NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      np = get_excitation_number(true,refa,confg,exct,ci,Iwork);
      alpha_index = ((np==0)?(0):(find_excitation(exct,unique_alpha[np]) +
                                  ph_struct.number_of_unique_smaller_than(np)[0]));
      confg.clear();
      for(int k=0, q=0; k<NAEB; k++) {
        q = occs[i][NAEA+k];
        if(q < NMO || q >= 2*NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      np = get_excitation_number(true,refb,confg,exct,ci,Iwork);
      beta_index = ((np==0)?(0):(find_excitation(exct,unique_beta[np]) +
                                  ph_struct.number_of_unique_smaller_than(np)[1]));
      ph_struct.add_configuration(alpha_index,beta_index,ci);
    }
  }
  comm.barrier();
  return ph_struct;
}

/*
 * Read trial wavefunction information from file.
*/
void WavefunctionFactory::getCommonInput(hdf_archive& dump, int NMO, int NAEA, int NAEB, int& ndets_to_read,
                                         std::vector<ComplexType>& ci, WALKER_TYPES walker_type, bool root)
{
  // check for consistency in parameters
  std::vector<int> dims(5);
  if(!dump.readEntry(dims,"dims")) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Problems reading dims. \n";
    APP_ABORT("");
  }
  if(NMO != dims[0]) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Inconsistent NMO . \n";
    APP_ABORT("");
  }
  if(NAEA != dims[1]) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Inconsistent  NAEA. \n";
    APP_ABORT("");
  }
  if(NAEB != dims[2]) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Inconsistent  NAEB. \n";
    APP_ABORT("");
  }
  if(walker_type != dims[3]) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Inconsistent  walker_type. \n";
    APP_ABORT("");
  }
  if(ndets_to_read < 1) ndets_to_read = dims[4];
  app_log() << " - Number of determinants in trial wavefunction: " << ndets_to_read << "\n";
  if(ndets_to_read > dims[4]) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Inconsistent  ndets_to_read. \n";
    APP_ABORT("");
  }
  ci.resize(ndets_to_read);
  if(!dump.readEntry(ci, "ci_coeffs")) {
    app_error()<<" Error in WavefunctionFactory::getCommonInput(): Problems reading ci_coeffs. \n";
    APP_ABORT("");
  }
  app_log() << " - Coefficient of first determinant: " << ci[0] << "\n";
}

/*
 * Read Initial walker from file.
*/
void WavefunctionFactory::getInitialGuess(hdf_archive& dump, std::string& name, int NMO, int NAEA, int NAEB, WALKER_TYPES walker_type)
{
  auto guess = initial_guess.find(name);
  if( guess == initial_guess.end() ) {
    auto newg = initial_guess.insert(
              std::make_pair(name,boost::multi::array<ComplexType,3>({2,NMO,NAEA})));
    if(!newg.second)
      APP_ABORT(" Error: Problems adding new initial guess. \n");
    using ma::conj;
    std::fill_n((newg.first)->second.origin(),2*NMO*NAEA,ComplexType(0.0,0.0));
    {
      boost::multi::array_ref<ComplexType,2> psi0((newg.first)->second.origin(),{NMO,NAEA});
      if(!dump.readEntry(psi0, "Psi0_alpha")) {
        app_error()<<" Error in WavefunctionFactory: Initial wavefunction Psi0_alpha not found. \n";
        APP_ABORT("");
      }
    }
    if(walker_type==COLLINEAR) {
      boost::multi::array_ref<ComplexType,2> psi0((newg.first)->second.origin()+NMO*NAEA,{NMO,NAEB});
      if(!dump.readEntry(psi0, "Psi0_beta")) {
        app_error()<<" Error in WavefunctionFactory: Initial wavefunction Psi0_beta not found. \n";
        APP_ABORT("");
      }
    }
  } else
    APP_ABORT(" Error: Problems adding new initial guess, already exists. \n");
}


/*
 * Helper function to get HamOps object from file or from scratch.
*/
HamiltonianOperations WavefunctionFactory::getHamOps(bool read, hdf_archive& dump, WALKER_TYPES type, int NMO, int NAEA, int NAEB,
                                                      std::vector<PsiT_Matrix>& PsiT, TaskGroup_& TGprop, TaskGroup_& TGwfn,
                                                      RealType cutvn, RealType cutv2, int ndets_to_read, Hamiltonian& h)
{
  if(read) {
    return loadHamOps(dump,type,NMO,NAEA,NAEB,PsiT,TGprop,TGwfn,cutvn,cutv2);
  } else {
    return h.getHamiltonianOperations(false,ndets_to_read>1,type,PsiT, cutvn,cutv2,TGprop,TGwfn,dump);
  }
}



}

}
