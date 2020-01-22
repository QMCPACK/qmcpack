#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<fstream>
#include<map>
#include<utility>
#include<vector>
#include<numeric>

#include "Configuration.h"

#include "io/hdf_archive.h"
#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"
#include "AFQMC/Matrix/mpi3_shared_ma_proxy.hpp"
#include "AFQMC/HamiltonianOperations/THCOpsIO.hpp"

namespace qmcplusplus
{
namespace afqmc
{

// Right now, cutvn, cutvn2, TGprop and TGwfn are completely ignored.
// Note: addCoulomb only has meaning on the sparse hamiltonians, not in THC
HamiltonianOperations THCHamiltonian::getHamiltonianOperations(bool pureSD, bool addCoulomb,
            WALKER_TYPES type,std::vector<PsiT_Matrix>& PsiT, double cutvn,
            double cutv2,TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& hdf_restart)
{

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if(TGwfn.Global().root()) write_hdf = !hdf_restart.closed();
  //  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if(type==COLLINEAR)
    assert(PsiT.size()%2 == 0);
  int ndet = ((type!=COLLINEAR)?(PsiT.size()):(PsiT.size()/2));
  bool test_Luv = not useHalfRotatedMuv;

  //std::cout<<" test_Luv: " <<std::boolalpha <<test_Luv <<std::endl;

  if(ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  size_t gnmu,grotnmu,nmu,rotnmu,nmu0,nmuN,rotnmu0,rotnmuN;
  hdf_archive dump(TGwfn.Global());
  // right now only Node.root() reads
  if( TG.Node().root() ) {
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      app_error()<<" Error opening integral file in THCHamiltonian. \n";
      APP_ABORT("");
    }
    if(!dump.push("Hamiltonian",false)) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Group not Hamiltonian found. \n";
      APP_ABORT("");
    }
    if(!dump.push("THC",false)) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Group not THC found. \n";
      APP_ABORT("");
    }
  }
  if( TG.Global().root() ) {
    std::vector<int> Idata(3);
    if(!dump.readEntry(Idata,"dims")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading dims. \n";
      APP_ABORT("");
    }
    if(Idata[0] != NMO) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n";
      APP_ABORT(" Error: NMO differs from value in integral file. \n");
    }
    gnmu = size_t(Idata[1]);
    grotnmu = size_t(Idata[2]);
  }
  TG.Global().broadcast_value(gnmu);
  TG.Global().broadcast_value(grotnmu);

  if(test_Luv)
    grotnmu = gnmu;

  // setup partition, in general matrices are partitioned asize_t 'u'
  {
    int node_number = TGwfn.getLocalGroupNumber();
    int nnodes_prt_TG = TGwfn.getNGroupsPerTG();
    std::tie(rotnmu0,rotnmuN) = FairDivideBoundary(size_t(node_number),grotnmu,size_t(nnodes_prt_TG));
    rotnmu = rotnmuN-rotnmu0;

    node_number = TGprop.getLocalGroupNumber();
    nnodes_prt_TG = TGprop.getNGroupsPerTG();
    std::tie(nmu0,nmuN) = FairDivideBoundary(size_t(node_number),gnmu,size_t(nnodes_prt_TG));
    nmu = nmuN-nmu0;
  }

  using shm_Vmatrix = mpi3_shared_ma_proxy<ValueType>;
  using shm_Cmatrix = mpi3_shared_ma_proxy<ComplexType>;

  // INCONSISTENT IN REAL BUILD!!!!!! FIX FIX FIX
  // Until I figure something else, rotPiu and rotcPua are not distributed because a full copy is needed
  size_t nel_ = PsiT[0].size(0) + ((type==CLOSED)?0:(PsiT[1].size(0)));
  shm_Vmatrix rotMuv(TG.Node(),{rotnmu,grotnmu},{grotnmu,grotnmu},{rotnmu0,0});
  shm_Cmatrix rotPiu(TG.Node(),{size_t(NMO),grotnmu});
  std::vector<shm_Cmatrix> rotcPua;
  rotcPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    rotcPua.emplace_back(shm_Cmatrix(TG.Node(),{grotnmu,nel_}));
  shm_Cmatrix Piu(TG.Node(),{size_t(NMO),nmu},{size_t(NMO),gnmu},{0,nmu0});
  shm_Vmatrix Luv(TG.Node(),{nmu,gnmu},{gnmu,gnmu},{nmu0,0});
  // right now only 1 reader. Use hyperslabs and parallel io later
  // read Half transformed first
  if(TG.Node().root()) {
    using ma::conj;
    if(not test_Luv) {
      /***************************************/
      // read full matrix, not distributed for now
      auto piu_ = rotPiu.get();
      if(!dump.readEntry(piu_,"HalfTransformedFullOrbitals")) {
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
    }
  }
  TG.global_barrier();

  if(TG.Node().root()) {
    /***************************************/
    typename shm_Cmatrix::ma_type piu_(Piu.get());
    hyperslab_proxy<typename shm_Cmatrix::ma_type,2> hslab(piu_,
                                                         Piu.global_size(),
                                                         Piu.shape(),
                                                         Piu.global_offset());
    if(!dump.readEntry(hslab,"Orbitals")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Orbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    typename shm_Vmatrix::ma_type luv_(Luv.get());
    hyperslab_proxy<typename shm_Vmatrix::ma_type,2> hslab2(luv_,
                                                         Luv.global_size(),
                                                         Luv.shape(),
                                                         Luv.global_offset());
    if(!dump.readEntry(hslab2,"Luv")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Luv. \n";
      APP_ABORT("");
    }
    /***************************************/
  }
  TG.global_barrier();

  boost::multi::array<ComplexType,2> v0({Piu.size(0),Piu.size(0)});
  if(TGprop.getNGroupsPerTG() > 1)
  {

    // TOO MUCH MEMORY, FIX FIX FIX!!!
    shm_Cmatrix Piu__(TG.Node(),{size_t(NMO),gnmu});
    shm_Vmatrix Luv__(TG.Node(),{gnmu,gnmu});
    if(TG.Node().root()) {
      auto luv_ = Luv__.get();
      if(!dump.readEntry(luv_,"Luv")) {
        app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading Orbitals. \n";
        APP_ABORT("");
      }
      auto piu_ = Piu__.get();
      if(!dump.readEntry(piu_,"Orbitals")) {
        app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading Orbitals. \n";
        APP_ABORT("");
      }
    }
    TG.Node().barrier();

    using ma::H;
    using ma::T;
    using ma::conj;
    size_t c0,cN,nc;
    std::tie(c0,cN) = FairDivideBoundary(size_t(TG.Global().rank()),gnmu,size_t(TG.Global().size()));
    nc = cN-c0;
    boost::multi::array<ComplexType,2> Tuv({gnmu,nc});
    boost::multi::array<ValueType,2> Muv({gnmu,nc});

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work
    ma::product(Luv__.get(),H(Luv__.get().sliced(c0,cN)),Muv);

    if(test_Luv) {
      for(int i=0; i<gnmu; ++i)
        std::copy_n(Muv[i].origin(),nc,rotMuv.get()[i].origin()+c0);
      TG.Node().barrier();
      if(TG.Node().root())
        TG.Cores().all_reduce_in_place_n(rotMuv.origin(),rotMuv.num_elements(),std::plus<>());
    }

    // since generating v0 takes some effort and temporary space,
    // v0(i,l) = -0.5*sum_j <i,j|j,l>
    //         = -0.5 sum_j,u,v ma::conj(Piu(i,u)) ma::conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v ma::conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) ma::conj(Piu(j,v))
    ma::product(H(Piu__.get()),Piu__.get()({0,long(NMO)},{long(c0),long(cN)}),Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for(size_t i=0; i<Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = ma::conj(*itT)*(*itM);
    boost::multi::array<ComplexType,2> T_({Tuv.size(1),size_t(NMO)});
    ma::product(T(Tuv),H(Piu__.get()),T_);
    ma::product(-0.5,T(T_),T(Piu__.get()({0,long(NMO)},{long(c0),long(cN)})),0.0,v0);

    // reduce over Global
    TG.Global().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());

  } else {
    // very simple partitioning until something more sophisticated is in place!!!
    using ma::H;
    using ma::T;
    using ma::conj;
    size_t c0,cN,nc;
    std::tie(c0,cN) = FairDivideBoundary(size_t(TG.Global().rank()),nmu,size_t(TG.Global().size()));
    nc = cN-c0;
    boost::multi::array<ComplexType,2> Tuv({gnmu,nc});
    boost::multi::array<ValueType,2> Muv({gnmu,nc});

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work
    ma::product(Luv.get(),H(Luv.get().sliced(c0,cN)),Muv);

    if(test_Luv) {
      for(int i=0; i<gnmu; ++i)
        std::copy_n(Muv[i].origin(),nc,rotMuv.get()[i].origin()+c0);
      TG.Node().barrier();
      if(TG.Node().root())
        TG.Cores().all_reduce_in_place_n(rotMuv.origin(),rotMuv.num_elements(),std::plus<>());
    }

    // since generating v0 takes some effort and temporary space,
    // v0(i,l) = -0.5*sum_j <i,j|j,l>
    //         = -0.5 sum_j,u,v ma::conj(Piu(i,u)) ma::conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v ma::conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) ma::conj(Piu(j,v))
    ma::product(H(Piu.get()),Piu.get()({0,long(NMO)},{long(c0),long(cN)}),Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for(size_t i=0; i<Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = ma::conj(*itT)*(*itM);
    boost::multi::array<ComplexType,2> T_({Tuv.size(1),size_t(NMO)});
    ma::product(T(Tuv),H(Piu.get()),T_);
    ma::product(-0.5,T(T_),T(Piu.get()({0,long(NMO)},{long(c0),long(cN)})),0.0,v0);

    // reduce over Global
    TG.Global().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());
  }
  TG.global_barrier();

  long naea_ = PsiT[0].size(0);
  long naeb_ = PsiT.back().size(0);

  // half-rotated Pia
  std::vector<shm_Cmatrix> cPua;
  cPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    cPua.emplace_back(shm_Cmatrix(TG.Node(),{nmu,nel_},{gnmu,nel_},{nmu0,0}));
  if(TG.Node().root()) {
    // simple
    using ma::H;
    if(type==COLLINEAR) {
      boost::multi::array<ComplexType,2> A({NMO,naea_});
      boost::multi::array<ComplexType,2> B({NMO,naeb_});
      for(int i=0; i<ndet; i++) {
        // cPua = H(Piu) * ma::conj(A)
        ma::Matrix2MA('T',PsiT[2*i],A);
        ma::product(H(Piu.get()),A,
                    cPua[i].get()({0,long(nmu)},{0,long(naea_)}));
        if(not test_Luv)
          ma::product(H(rotPiu.get()),A,
                      rotcPua[i].get()({0,long(grotnmu)},{0,long(naea_)}));
        ma::Matrix2MA('T',PsiT[2*i+1],B);
        ma::product(H(Piu.get()),B,
                    cPua[i].get()({0,long(nmu)},{naea_,long(nel_)}));
        if(not test_Luv)
          ma::product(H(rotPiu.get()),B,
                      rotcPua[i].get()({0,long(grotnmu)},{naea_,long(nel_)}));
      }
    } else {
      boost::multi::array<ComplexType,2> A({PsiT[0].size(1),PsiT[0].size(0)});
      for(int i=0; i<ndet; i++) {
        ma::Matrix2MA('T',PsiT[i],A);
        // cPua = H(Piu) * ma::conj(A)
        ma::product(H(Piu.get()),A,cPua[i].get());
        if(not test_Luv)
          ma::product(H(rotPiu.get()),A,rotcPua[i].get());
      }
    }
    if(test_Luv) {
      std::copy_n(Piu.origin(),Piu.num_elements(),rotPiu.origin());
      for(int i=0; i<ndet; i++)
        std::copy_n(cPua[i].origin(),cPua[i].num_elements(),rotcPua[i].origin());
    }
  }
  TG.node_barrier();

  ValueType E0 = OneBodyHamiltonian::NuclearCoulombEnergy +
                 OneBodyHamiltonian::FrozenCoreEnergy;

  std::vector<boost::multi::array<ComplexType,1>> hij;
  hij.reserve(ndet);
  int skp=((type==COLLINEAR)?1:0);
  for(int n=0, nd=0; n<ndet; ++n, nd+=(skp+1)) {
    check_wavefunction_consistency(type,&PsiT[nd],&PsiT[nd+skp],NMO,naea_,naeb_);
    hij.emplace_back(rotateHij(type,&PsiT[nd],&PsiT[nd+skp],OneBodyHamiltonian::H1));
  }

  // dense one body hamiltonian
  auto H1 = getH1();

//std::cout<<" nmu: " <<Luv.size(0) <<" " <<rotMuv.size(0) <<std::endl;

  if(write_hdf)
    writeTHCOps(hdf_restart,type,NMO,naea_,naeb_,ndet,TGprop,TGwfn,H1,
                rotPiu,rotMuv,Piu,Luv,v0,E0);

  return HamiltonianOperations(THCOps<ValueType>(TGwfn.TG_local(),NMO,naea_,naeb_,type,std::move(H1),
                                      std::move(hij),std::move(rotMuv),std::move(rotPiu),
                                      std::move(rotcPua),std::move(Luv),
                                      std::move(Piu),std::move(cPua),std::move(v0),E0));
}


}
}
