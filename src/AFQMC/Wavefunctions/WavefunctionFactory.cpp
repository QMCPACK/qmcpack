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

#include "io/hdf_archive.h"
#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Wavefunctions/NOMSD.hpp"
#include "AFQMC/HamiltonianOperations/HamOpsIO.hpp"

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
  m_param.add(restart_file,"restart_file","std::string");
  m_param.add(write_trial_density_matrix,"trial_density_matrix","std::string");
  m_param.add(cutv2,"cutoff","double");
  m_param.add(initialDet,"initialDetType","int");
  m_param.add(starting_det,"starting_det","std:string");
  m_param.add(ndets_to_read,"ndet","int");
  m_param.put(cur);

  AFQMCInfo& AFinfo = InfoMap[info];
  auto NCE = h.getNuclearCoulombEnergy();

  int NMO = AFinfo.NMO;
  int NAEA = AFinfo.NAEA;
  int NAEB = AFinfo.NAEB;

  std::vector<ComplexType> ci;
  std::vector<PsiT_Matrix> PsiT;
  std::vector<int> excitations;  
  std::string wfn_type;
  read_wavefunction(filename,ndets_to_read,wfn_type,walker_type,TGwfn.Node(),
                    NMO,NAEA,NAEB,PsiT,ci,excitations);

  using Alloc = boost::mpi3::intranode::allocator<ComplexType>;
  if(type=="msd") {
    if(wfn_type == "occ") {
      int NEL = (walker_type==NONCOLLINEAR)?(NAEA+NAEB):NAEA;  
      int N_ = (walker_type==NONCOLLINEAR)?2*NMO:NMO;
      ComplexType one(1.0,0.0);
      auto it = excitations.begin();  
      if(walker_type==COLLINEAR) 
        PsiT.reserve(2*ndets_to_read);  
      else
        PsiT.reserve(ndets_to_read);  
      for(int i=0; i<ndets_to_read; i++) {
        PsiT.emplace_back(PsiT_Matrix({NEL,N_},{0,0},1,Alloc(TGwfn.Node())));
        if(TGwfn.Node().root())  
          for(int k=0; k<NEL; k++) 
            PsiT.back().emplace_back({k,*it++},one); 
        if(walker_type==COLLINEAR) {
          PsiT.emplace_back(PsiT_Matrix({NAEB,NMO},{0,0},1,Alloc(TGwfn.Node())));
          if(TGwfn.Node().root())  
            for(int k=0; k<NAEB; k++) 
              PsiT.back().emplace_back({k,(*it++)-NMO},one); 
        }
      }  
    } else if(wfn_type == "mixed") {
      // In this case, PsiT is the orbital matrix
      APP_ABORT("Finish WavefunctionFactory \n");
    } else if(wfn_type == "matrix") {
      // nothing to do here  
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
    auto HOps(h.getHamiltonianOperations(wfn_type == "occ",walker_type,PsiT,cutvn,cutv2,TGprop,TGwfn,dump));  
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
                std::make_pair(name,boost::multi_array<ComplexType,3>(extents[2][NMO][NAEA])));
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

    //return Wavefunction{}; 
    return Wavefunction(NOMSD(AFinfo,cur,TGwfn,std::move(HOps),std::move(ci),std::move(PsiT),
                        walker_type,NCE,targetNW)); 
  } else if(type == "fastmsd") {
    app_error()<<" Error: Wavefunction type FastMSD not yet implemented. \n"; 
    APP_ABORT(" Error: Wavefunction type FastMSD not yet implemented. \n");
    return Wavefunction{}; 
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
  if(type=="msd") {

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
                std::make_pair(name,boost::multi_array<ComplexType,3>(extents[2][NMO][NAEA])));
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
  } else if(type == "fastmsd") {
    app_error()<<" Error: Wavefunction type FastMSD not yet implemented. \n";
    APP_ABORT(" Error: Wavefunction type FastMSD not yet implemented. \n");
    return Wavefunction{};
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
