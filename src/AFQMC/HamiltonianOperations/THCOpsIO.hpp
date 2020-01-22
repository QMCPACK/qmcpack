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

#ifndef QMCPLUSPLUS_AFQMC_THCOPSIO_HPP
#define QMCPLUSPLUS_AFQMC_THCOPSIO_HPP

#include<fstream>

#include "type_traits/container_traits_multi.h"
#include "io/hdf_multi.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Matrix/ma_hdf5_readers.hpp"
#include "AFQMC/Numerics/ma_blas_extensions.hpp"

#include "AFQMC/HamiltonianOperations/THCOps.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"

namespace qmcplusplus
{
namespace afqmc
{

// Some code duplication with THCHamiltonian class.
template<typename T>
THCOps<T> loadTHCOps(hdf_archive& dump, WALKER_TYPES type, int NMO, int NAEA, int NAEB, std::vector<PsiT_Matrix>& PsiT, TaskGroup_& TGprop, TaskGroup_& TGwfn, RealType cutvn, RealType cutv2)
{

#if defined(MIXED_PRECISION)
  using SpT = typename to_single_precision<T>::value_type;
  using SpC = typename to_single_precision<ComplexType>::value_type;
#else
  using SpT = T;
  using SpC = ComplexType;
#endif

  using shm_Vmatrix = mpi3_shared_ma_proxy<T>;
  using shm_Cmatrix = mpi3_shared_ma_proxy<ComplexType>;

  if(type==COLLINEAR)
    assert(PsiT.size()%2 == 0);
  int ndet = ((type!=COLLINEAR)?(PsiT.size()):(PsiT.size()/2));
  if(ndet > 1) {
    app_error()<<" Error in loadTHCOps: ndet > 1 not yet implemented in THCOps." <<std::endl;
    APP_ABORT("");
  }

  // fix later for multidet case
  std::vector<int> dims(10);
  ValueType E0;
  std::size_t gnmu,grotnmu,nmu,rotnmu,nmu0,nmuN,rotnmu0,rotnmuN;

  // read from HDF

  if(!dump.push("HamiltonianOperations",false)) {
    app_error()<<" Error in loadTHCOps: Group HamiltonianOperations not found. \n";
    APP_ABORT("");
  }
  if(!dump.push("THCOps",false)) {
    app_error()<<" Error in loadTHCOps: Group THCOps not found. \n";
    APP_ABORT("");
  }
  if(TGwfn.Global().root()) {
    if(!dump.readEntry(dims,"dims")) {
      app_error()<<" Error in loadTHCOps: Problems reading dataset. \n";
      APP_ABORT("");
    }
    assert(dims.size()==7);
    if(dims[0] != NMO) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: NMO. \n";
      APP_ABORT("");
    }
    if(dims[1] != NAEA) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: NAEA. \n";
      APP_ABORT("");
    }
    if(dims[2] != NAEB) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: NAEB. \n";
      APP_ABORT("");
    }
    if(dims[3] != ndet) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: ndet. \n";
      APP_ABORT("");
    }
    if(type == CLOSED && dims[4] != 1) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if(type == COLLINEAR && dims[4] != 2) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    if(type == NONCOLLINEAR && dims[4] != 3) {
      app_error()<<" Error in loadTHCOps: Inconsistent data in file: walker_type. \n";
      APP_ABORT("");
    }
    std::vector<ValueType> et;
    if(!dump.readEntry(et,"E0")) {
      app_error()<<" Error in loadTHCOps: Problems reading dataset. \n";
      APP_ABORT("");
    }
    E0=et[0];
  }
  TGwfn.Global().broadcast_n(dims.data(),7);
  TGwfn.Global().broadcast_value(E0);
  gnmu = size_t(dims[5]);
  grotnmu = size_t(dims[6]);

  // setup partition, in general matrices are partitioned along 'u'
  {
    int node_number = TGwfn.getLocalGroupNumber();
    int nnodes_prt_TG = TGwfn.getNGroupsPerTG();
    std::tie(rotnmu0,rotnmuN) = FairDivideBoundary(std::size_t(node_number),grotnmu,std::size_t(nnodes_prt_TG));
    rotnmu = rotnmuN-rotnmu0;

    node_number = TGprop.getLocalGroupNumber();
    nnodes_prt_TG = TGprop.getNGroupsPerTG();
    std::tie(nmu0,nmuN) = FairDivideBoundary(std::size_t(node_number),gnmu,std::size_t(nnodes_prt_TG));
    nmu = nmuN-nmu0;
  }

  // read 1-body hamiltonian and exchange potential (v0)
  boost::multi::array<ValueType,2> H1({NMO,NMO});
  boost::multi::array<ComplexType,2> v0({NMO,NMO});
  if(TGwfn.Global().root()) {
    if(!dump.readEntry(H1,"H1")) {
      app_error()<<" Error in loadTHCOps: Problems reading dataset. \n";
      APP_ABORT("");
    }
    if(!dump.readEntry(v0,"v0")) {
      app_error()<<" Error in loadTHCOps: Problems reading dataset. \n";
      APP_ABORT("");
    }
  }
  TGwfn.Global().broadcast_n(H1.origin(),H1.num_elements());
  TGwfn.Global().broadcast_n(v0.origin(),v0.num_elements());

  // Until I figure something else, rotPiu and rotcPua are not distributed because a full copy is needed
  size_t nel_ = ((type==CLOSED)?NAEA:(NAEA+NAEB));
  shm_Vmatrix rotMuv(TGwfn.Node(),{rotnmu,grotnmu},{grotnmu,grotnmu},{rotnmu0,0});
  shm_Cmatrix rotPiu(TGwfn.Node(),{size_t(NMO),grotnmu});
  shm_Cmatrix Piu(TGwfn.Node(),{size_t(NMO),nmu},{size_t(NMO),gnmu},{0,nmu0});
  shm_Vmatrix Luv(TGwfn.Node(),{nmu,gnmu},{gnmu,gnmu},{nmu0,0});

  // read Half transformed first
  if(TGwfn.Node().root()) {
    /***************************************/
    auto rpiu_ = rotPiu.get();
    if(!dump.readEntry(rpiu_,"HalfTransformedFullOrbitals")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading HalfTransformedFullOrbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    typename shm_Vmatrix::ma_type muv_(rotMuv.get());
    hyperslab_proxy<typename shm_Vmatrix::ma_type,2> hslab(muv_,
                                                           rotMuv.global_size(),
                                                           rotMuv.shape(),
                                                           rotMuv.global_offset());
    if(!dump.readEntry(hslab,"HalfTransformedMuv")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                <<" Problems reading HalfTransformedMuv. \n";
      APP_ABORT("");
    }
    /***************************************/
    typename shm_Cmatrix::ma_type piu_(Piu.get());
    hyperslab_proxy<typename shm_Cmatrix::ma_type,2> hslab2(piu_,
                                                         Piu.global_size(),
                                                         Piu.shape(),
                                                         Piu.global_offset());
    if(!dump.readEntry(hslab2,"Orbitals")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Orbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    typename shm_Vmatrix::ma_type luv_(Luv.get());
    hyperslab_proxy<typename shm_Vmatrix::ma_type,2> hslab3(luv_,
                                                         Luv.global_size(),
                                                         Luv.shape(),
                                                         Luv.global_offset());
    if(!dump.readEntry(hslab3,"Luv")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Luv. \n";
      APP_ABORT("");
    }
    /***************************************/
  }
  TGwfn.global_barrier();

  // half-rotated Pia
  std::vector<shm_Cmatrix> rotcPua;
  rotcPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    rotcPua.emplace_back(shm_Cmatrix(TGwfn.Node(),{grotnmu,nel_}));
  std::vector<shm_Cmatrix> cPua;
  cPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    cPua.emplace_back(shm_Cmatrix(TGwfn.Node(),{nmu,nel_},{gnmu,nel_},{nmu0,0}));
  if(TGwfn.Node().root()) {
    // simple
    using ma::H;
    if(type==COLLINEAR) {
      boost::multi::array<ComplexType,2> A({NMO,NAEA});
      boost::multi::array<ComplexType,2> B({NMO,NAEB});
      for(int i=0; i<ndet; i++) {
        // cPua = H(Piu) * conj(A)
        ma::Matrix2MA('T',PsiT[2*i],A);
        auto&& cPua_i(cPua[i].get());
        auto&& rotcPua_i(rotcPua[i].get());
        ma::product(H(Piu.get()),A,cPua_i(cPua_i.extension(0),{0,NAEA}));
        ma::product(H(rotPiu.get()),A,rotcPua_i(cPua_i.extension(0),{0,NAEA}));
        ma::Matrix2MA('T',PsiT[2*i+1],B);
        ma::product(H(Piu.get()),B,cPua_i(cPua_i.extension(0),{NAEA,NAEA+NAEB}));
        ma::product(H(rotPiu.get()),B,rotcPua_i(cPua_i.extension(0),{NAEA,NAEA+NAEB}));
      }
    } else {
      boost::multi::array<ComplexType,2> A({PsiT[0].size(1),PsiT[0].size(0)});
      for(int i=0; i<ndet; i++) {
        ma::Matrix2MA('T',PsiT[i],A);
        // cPua = H(Piu) * conj(A)
        ma::product(H(Piu.get()),A,cPua[i].get());
        ma::product(H(rotPiu.get()),A,rotcPua[i].get());
      }
    }
  }
  TGwfn.node_barrier();

  // rotated 1 body hamiltonians
  std::vector<boost::multi::array<ComplexType,1>> hij;
  hij.reserve(ndet);
  int skp=((type==COLLINEAR)?1:0);
  for(int n=0, nd=0; n<ndet; ++n, nd+=(skp+1)) {
    check_wavefunction_consistency(type,&PsiT[nd],&PsiT[nd+skp],NMO,NAEA,NAEB);
    hij.emplace_back(rotateHij(type,&PsiT[nd],&PsiT[nd+skp],H1));
  }

  return THCOps<T>(TGwfn.TG_local(),NMO,NAEA,NAEB,type,std::move(H1),
                                      std::move(hij),std::move(rotMuv),std::move(rotPiu),
                                      std::move(rotcPua),std::move(Luv),
                                      std::move(Piu),std::move(cPua),std::move(v0),E0);
}

// single writer right now
template<class shm_Vmatrix,
         class shm_Cmatrix>
inline void writeTHCOps(hdf_archive& dump, WALKER_TYPES type, int NMO, int NAEA, int NAEB, int ndet,
                              TaskGroup_& TGprop, TaskGroup_& TGwfn,
                              boost::multi::array<ValueType,2> & H1,
                              shm_Cmatrix & rotPiu,
                              shm_Vmatrix & rotMuv,
                              shm_Cmatrix & Piu,
                              shm_Vmatrix & Luv,
                              boost::multi::array<ComplexType,2> & v0,
                              ValueType E0)
{

  if(TGwfn.Global().root()) {
    dump.push("HamiltonianOperations");
    dump.push("THCOps");
    std::vector<int> dims{NMO,NAEA,NAEB,ndet,type,int(Luv.global_size(0)),int(rotMuv.global_size(0))};
    dump.write(dims,"dims");
    std::vector<ValueType> et{E0};
    dump.write(et,"E0");
    dump.write(H1,"H1");
    dump.write(v0,"v0");
    auto rotPiu_(rotPiu.get());
    auto Piu_(Piu.get());
    auto rotMuv_(rotMuv.get());
    auto Luv_(Luv.get());
    dump.write(rotPiu_,"HalfTransformedFullOrbitals");
    ma_hdf5::write_distributed_MA(rotMuv_,rotMuv.global_offset(),rotMuv.global_size(),
                                  dump,"HalfTransformedMuv",TGwfn);
    ma_hdf5::write_distributed_MA(Piu_,Piu.global_offset(),Piu.global_size(),
                                  dump,"Orbitals",TGprop);
    ma_hdf5::write_distributed_MA(Luv_,Luv.global_offset(),Luv.global_size(),
                                  dump,"Luv",TGprop);
    dump.pop();
    dump.pop();
  } else {
    auto Piu_(Piu.get());
    auto rotMuv_(rotMuv.get());
    auto Luv_(Luv.get());
    ma_hdf5::write_distributed_MA(rotMuv_,rotMuv.global_offset(),rotMuv.global_size(),
                                  dump,"HalfTransformedMuv",TGwfn);
    ma_hdf5::write_distributed_MA(Piu_,Piu.global_offset(),Piu.global_size(),
                                  dump,"Orbitals",TGprop);
    ma_hdf5::write_distributed_MA(Luv_,Luv.global_offset(),Luv.global_size(),
                                  dump,"Luv",TGprop);
  }
  TGwfn.Global().barrier();

}

}
}

#endif
