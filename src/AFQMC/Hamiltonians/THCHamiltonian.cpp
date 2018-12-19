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
HamiltonianOperations THCHamiltonian::getHamiltonianOperations(bool pureSD, 
            WALKER_TYPES type,std::vector<PsiT_Matrix>& PsiT, double cutvn, 
            double cutv2,TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& hdf_restart)
{

  // hack until parallel hdf is in place
  bool write_hdf = false;
  if(TGwfn.Global().root()) write_hdf = (hdf_restart.file_id != hdf_archive::is_closed);
  TGwfn.Global().broadcast_value(write_hdf);

  if(type==COLLINEAR)
    assert(PsiT.size()%2 == 0);
  int ndet = ((type!=COLLINEAR)?(PsiT.size()):(PsiT.size()/2));
  bool test_Luv = not useHalfRotatedMuv; 

  if(ndet > 1)
    APP_ABORT("Error: ndet > 1 not yet implemented in THCHamiltonian::getHamiltonianOperations.\n");

  std::size_t gnmu,grotnmu,nmu,rotnmu,nmu0,nmuN,rotnmu0,rotnmuN; 
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
    if(!dump.read(Idata,"dims")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading dims. \n";
      APP_ABORT("");
    }
    if(Idata[0] != NMO) {
      app_error()<<" ERROR: NMO differs from value in integral file. \n";
      APP_ABORT(" Error: NMO differs from value in integral file. \n");
    }
    gnmu = std::size_t(Idata[1]);  
    grotnmu = std::size_t(Idata[2]);  
  }
  TG.Global().broadcast_value(gnmu);
  TG.Global().broadcast_value(grotnmu);

  if(test_Luv)
    grotnmu = gnmu;

  // setup partition, in general matrices are partitioned along 'u' 
  {
    int node_number = TGwfn.getLocalNodeNumber();
    int nnodes_prt_TG = TGwfn.getNNodesPerTG();
    std::tie(rotnmu0,rotnmuN) = FairDivideBoundary(std::size_t(node_number),grotnmu,std::size_t(nnodes_prt_TG));
    rotnmu = rotnmuN-rotnmu0;

    node_number = TGprop.getLocalNodeNumber();
    nnodes_prt_TG = TGprop.getNNodesPerTG();
    std::tie(nmu0,nmuN) = FairDivideBoundary(std::size_t(node_number),gnmu,std::size_t(nnodes_prt_TG));
    nmu = nmuN-nmu0;
  }

  using shm_Vmatrix = mpi3_shared_ma_proxy<ValueType>;
  using shm_Cmatrix = mpi3_shared_ma_proxy<ComplexType>;

  // INCONSISTENT IN REAL BUILD!!!!!! FIX FIX FIX
  // Until I figure something else, rotPiu and rotcPua are not distributed because a full copy is needed
  size_t nel_ = ((type==CLOSED)?NAEA:(NAEA+NAEB));
  shm_Cmatrix rotMuv(TG.Node(),{rotnmu,grotnmu},{grotnmu,grotnmu},{rotnmu0,0});
  shm_Cmatrix rotPiu(TG.Node(),{size_t(NMO),grotnmu});
  std::vector<shm_Cmatrix> rotcPua;
  rotcPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    rotcPua.emplace_back(shm_Cmatrix(TG.Node(),{grotnmu,nel_}));
  shm_Vmatrix Piu(TG.Node(),{size_t(NMO),nmu},{size_t(NMO),gnmu},{0,nmu0});
  shm_Vmatrix Luv(TG.Node(),{nmu,gnmu},{gnmu,gnmu},{nmu0,0});
  // right now only 1 reader. Use hyperslabs and parallel io later 
  // read Half transformed first
  if(TG.Node().root()) {
    using std::conj;
    if(not test_Luv) {
      /***************************************/
      // read full matrix, not distributed for now
      auto piu_ = rotPiu.get();
      if(!dump.read(piu_,"HalfTransformedFullOrbitals")) {
        app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading HalfTransformedFullOrbitals. \n";
        APP_ABORT("");
      }
      /***************************************/
      typename shm_Cmatrix::ma_type muv_(rotMuv.get());
      hyperslab_proxy<typename shm_Cmatrix::ma_type,2> hslab(muv_,
                                                           rotMuv.global_shape(),
                                                           rotMuv.shape(),
                                                           rotMuv.offset());
      if(!dump.read(hslab,"HalfTransformedMuv")) {
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
    typename shm_Vmatrix::ma_type piu_(Piu.get());
    hyperslab_proxy<typename shm_Vmatrix::ma_type,2> hslab(piu_,
                                                         Piu.global_shape(),
                                                         Piu.shape(),
                                                         Piu.offset());
    if(!dump.read(hslab,"Orbitals")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Orbitals. \n";
      APP_ABORT("");
    }
    /***************************************/
    typename shm_Vmatrix::ma_type luv_(Luv.get());
    hyperslab_proxy<typename shm_Vmatrix::ma_type,2> hslab2(luv_,
                                                         Luv.global_shape(),
                                                         Luv.shape(),
                                                         Luv.offset());
    if(!dump.read(hslab2,"Luv")) {
      app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                 <<" Problems reading Luv. \n";
      APP_ABORT("");
    }
    /***************************************/
  }
  TG.global_barrier();

  boost::multi_array<ComplexType,2> v0(extents[Piu.shape()[0]][Piu.shape()[0]]);
  if(TGprop.getNNodesPerTG() > 1) 
  {

    shm_Vmatrix Piu__(TG.Node(),{size_t(NMO),gnmu});
    shm_Vmatrix Luv__(TG.Node(),{gnmu,gnmu});
    if(TG.Node().root()) {
      auto luv_ = Luv__.get();
      if(!dump.read(luv_,"Luv")) {
        app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading Orbitals. \n";
        APP_ABORT("");
      }
      auto piu_ = Piu__.get();
      if(!dump.read(piu_,"Orbitals")) {
        app_error()<<" Error in THCHamiltonian::getHamiltonianOperations():"
                   <<" Problems reading Orbitals. \n";
        APP_ABORT("");
      }
    }
    TG.Node().barrier();

    using ma::H;
    using ma::T;
    using std::conj;
    size_t c0,cN,nc;
    std::tie(c0,cN) = FairDivideBoundary(size_t(TG.Global().rank()),gnmu,size_t(TG.Global().size()));
    nc = cN-c0;
    boost::multi_array<ValueType,2> Tuv(extents[gnmu][nc]);
    boost::multi_array<ValueType,2> Muv(extents[gnmu][nc]);

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work 
    ma::product(Luv__.get(),H(Luv__.get()[indices[range_t(c0,cN)][range_t()]]),Muv);

    // since generating v0 takes some effort and temporary space, 
    // v0(i,l) = -0.5*sum_j <i,j|j,l> 
    //         = -0.5 sum_j,u,v conj(Piu(i,u)) conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) conj(Piu(j,v))  
    ma::product(H(Piu__.get()),Piu__.get()[indices[range_t()][range_t(c0,cN)]],Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for(size_t i=0; i<Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = conj(*itT)*(*itM);
    boost::multi_array<ValueType,2> T_(extents[Tuv.shape()[1]][size_t(NMO)]);
    ma::product(T(Tuv),H(Piu__.get()),T_);
    ma::product(-0.5,T(T_),T(Piu__.get()[indices[range_t()][range_t(c0,cN)]]),0.0,v0);

    // reduce over Global  
    TG.Global().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());

/*
    using ma::H;
    using ma::T;
    using std::conj;

    // split {u,v} space over all cores.
    int node_number = TGprop.getLocalNodeNumber();
    int nnodes_per_TG = TGprop.getNNodesPerTG();

    std::size_t eqv_node_rank = std::size_t(TGprop.getNodeID()/nnodes_per_TG);
    std::size_t eqv_node_size = std::size_t(TGprop.getTotalNodes()/nnodes_per_TG);
    std::size_t eqv_pos = eqv_node_rank*std::size_t(TGprop.getTotalCores()) + TGprop.getCoreID();
    std::size_t eqv_size = eqv_node_size*std::size_t(TGprop.getTotalCores());;

    // communicator with cores in same position of all nodes in a TG 
    boost::mpi3::communicator core_team(TGprop.TG().split(TGprop.getCoreID()));

    std::fill_n(v0.origin(),v0.num_elements(),ValueType(0.0));

    for(int ni=0; ni<nnodes_per_TG; ++ni) { // loop over nodes in TG
      // split interpolating points of node ni over equivalent nodes and bcast to TG
     
      // find partition of node ni
      std::size_t curr_v0,curr_vN;
      std::tie(curr_v0,curr_vN) = FairDivideBoundary(std::size_t(ni),gnmu,std::size_t(nnodes_per_TG));

      // split range over cores and over equivalent nodes
      std::size_t vi0,viN,nv;
      std::tie(vi0,viN) = FairDivideBoundary(eqv_pos,curr_vN-curr_v0,eqv_size);
      nv = viN-vi0;

      // bcast vi0-vN over core-team and accumulate contribution
      boost::multi_array<ValueType,2> Luv_(extents[nv][gnmu]);      
      boost::multi_array<ValueType,2> Piu_(extents[size_t(NMO)][nv]);      
      if(node_number == ni) {
        std::copy_n((Luv.get())[vi0].origin(),nv*Luv.shape()[1],Luv_.origin());
        for(int i=0; i<NMO; ++i) 
          std::copy_n( (Piu.get())[i].origin()+vi0, nv, Piu_[i].origin() );
      }
      core_team.broadcast_n(Luv_.origin(),Luv_.num_elements(),ni);
      core_team.broadcast_n(Piu_.origin(),Piu_.num_elements(),ni);

      // accumulate contribution
      boost::multi_array<ValueType,2> Muv(extents[nmu][nv]);  
      boost::multi_array<ValueType,2> Tuv(extents[nmu][nv]);  
      ma::product(Luv.get(),H(Luv_),Muv);

      if(test_Luv)
        for(int i=0; i<nmu; ++i)
          std::copy_n(Muv[i].origin(),nv,rotMuv.get()[i].origin()+curr_v0+vi0);

      // since generating v0 takes some effort and temporary space, 
      // v0(i,l) = -0.5*sum_j <i,j|j,l> 
      //         = -0.5 sum_j,u,v conj(Piu(i,u)) conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
      //         = -0.5 sum_u,v conj(Piu(i,u)) W(u,v) Piu(l,u), where
      // W(u,v) = Muv(u,v) * sum_j Piu(j,u) conj(Piu(j,v))  
      ma::product(H(Piu.get()),Piu_,Tuv);
      auto itM = Muv.origin();
      auto itT = Tuv.origin();
      for(size_t i=0; i<Muv.num_elements(); ++i, ++itT, ++itM)
        *(itT) = conj(*itT)*(*itM);   
      boost::multi_array<ValueType,2> T_(extents[Tuv.shape()[1]][size_t(NMO)]);
      ma::product(T(Tuv),H(Piu.get()),T_);
      ma::product(-0.5,T(T_),T(Piu_),1.0,v0); 
    }
    if(test_Luv) // need a new communicator to reduce rotMuv if test_Luv is true
       APP_ABORT("Finish this. \n\n\n"); 
    // reduce over Global  
    TG.Global().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());
*/
  } else {  
    // very simple partitioning until something more sophisticated is in place!!!
    using ma::H;
    using ma::T;
    using std::conj;
    size_t c0,cN,nc;
    std::tie(c0,cN) = FairDivideBoundary(size_t(TG.Global().rank()),nmu,size_t(TG.Global().size()));
    nc = cN-c0;
    boost::multi_array<ValueType,2> Tuv(extents[gnmu][nc]);
    boost::multi_array<ValueType,2> Muv(extents[gnmu][nc]);

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work 
    ma::product(Luv.get(),H(Luv.get()[indices[range_t(c0,cN)][range_t()]]),Muv);

    if(test_Luv)
      for(int i=0; i<gnmu; ++i)
        std::copy_n(Muv[i].origin(),nc,rotMuv.get()[i].origin()+c0);

    // since generating v0 takes some effort and temporary space, 
    // v0(i,l) = -0.5*sum_j <i,j|j,l> 
    //         = -0.5 sum_j,u,v conj(Piu(i,u)) conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) conj(Piu(j,v))  
    ma::product(H(Piu.get()),Piu.get()[indices[range_t()][range_t(c0,cN)]],Tuv);
    auto itM = Muv.origin();
    auto itT = Tuv.origin();
    for(size_t i=0; i<Muv.num_elements(); ++i, ++itT, ++itM)
      *(itT) = conj(*itT)*(*itM);           
    boost::multi_array<ValueType,2> T_(extents[Tuv.shape()[1]][size_t(NMO)]);
    ma::product(T(Tuv),H(Piu.get()),T_);
    ma::product(-0.5,T(T_),T(Piu.get()[indices[range_t()][range_t(c0,cN)]]),0.0,v0);

    // reduce over Global  
    TG.Global().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());

/*
    using ma::H;
    using std::conj;
    int c0,cN;
    std::tie(c0,cN) = FairDivideBoundary(TG.Node().rank(),int(nmu),TG.Node().size()); 
    int k0,kN;
    std::tie(k0,kN) = FairDivideBoundary(TG.Node().rank(),int(Piu.shape()[0]),TG.Node().size()); 
    // disribute globally !!!!! 
    shm_Vmatrix Tuv(TG.Node(),{nmu,nmu});
    shm_Vmatrix Muv(TG.Node(),{nmu,nmu});

    // Muv = Luv * H(Luv)
    // This can benefit from 2D split of work 
    ma::product(Luv.get()[indices[range_t(c0,cN)][range_t()]],H(Luv.get()),
                Muv.get()[indices[range_t(c0,cN)][range_t()]]);

    if(test_Luv)  
      std::copy_n(Muv.get()[c0].origin(),(cN-c0)*Muv.shape()[1],
                  rotMuv.get()[c0].origin());

    // since generating v0 takes some effort and temporary space, 
    // v0(i,l) = -0.5*sum_j <i,j|j,l> 
    //         = -0.5 sum_j,u,v conj(Piu(i,u)) conj(Piu(j,v)) Muv Piu(j,u) Piu(l,v)
    //         = -0.5 sum_u,v conj(Piu(i,u)) W(u,v) Piu(l,u), where
    // W(u,v) = Muv(u,v) * sum_j Piu(j,u) conj(Piu(j,v))  
    ma::product(H(Piu.get()[indices[range_t()][range_t(c0,cN)]]),Piu.get(),
                Tuv.get()[indices[range_t(c0,cN)][range_t()]]);
    auto itT = Tuv.get()[c0].origin();
    auto itM = Muv.get()[c0].origin();
    for(size_t i=0; i<(cN-c0)*nmu; ++i, ++itT, ++itM)
      *itT *= conj(*itM);   
    TG.node_barrier();
    boost::multi_array<ValueType,2> T_(extents[kN-k0][Tuv.shape()[1]]);
    ma::product(Piu.get()[indices[range_t(k0,kN)][range_t()]],Tuv.get(),
                T_);
    ma::product(T_,H(Piu.get()),v0[indices[range_t(k0,kN)][range_t()]]); 
    for(int i=k0; i<kN; i++)
      for(int j=0; j<v0.shape()[1]; j++)
        v0[i][j] = -0.5*conj(v0[i][j]);
    TG.Node().all_reduce_in_place_n(v0.origin(),v0.num_elements(),std::plus<>());
*/
  }
  TG.global_barrier();

  // half-rotated Pia
  std::vector<shm_Cmatrix> cPua;
  cPua.reserve(ndet);
  for(int i=0; i<ndet; i++)
    cPua.emplace_back(shm_Cmatrix(TG.Node(),{nmu,nel_},{gnmu,nel_},{nmu0,0}));  
  if(TG.Node().root()) {  
    // simple
    using ma::H;
    if(type==COLLINEAR) {
      boost::multi_array<ComplexType,2> A(extents[NMO][NAEA]); 
      boost::multi_array<ComplexType,2> B(extents[NMO][NAEB]); 
      for(int i=0; i<ndet; i++) {
        // cPua = H(Piu) * conj(A)
        csr::CSR2MA('T',PsiT[2*i],A);
        ma::product(H(Piu.get()),A,cPua[i].get()[indices[range_t()][range_t(0,NAEA)]]);
        if(not test_Luv)
          ma::product(H(rotPiu.get()),A,rotcPua[i].get()[indices[range_t()][range_t(0,NAEA)]]);
        csr::CSR2MA('T',PsiT[2*i+1],B);
        ma::product(H(Piu.get()),B,cPua[i].get()[indices[range_t()][range_t(NAEA,NAEA+NAEB)]]);
        if(not test_Luv)
          ma::product(H(rotPiu.get()),B,rotcPua[i].get()[indices[range_t()][range_t(NAEA,NAEA+NAEB)]]);
      }
    } else {
      boost::multi_array<ComplexType,2> A(extents[PsiT[0].shape()[1]][PsiT[0].shape()[0]]); 
      for(int i=0; i<ndet; i++) {
        csr::CSR2MA('T',PsiT[i],A);
        // cPua = H(Piu) * conj(A)
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

  std::vector<boost::multi_array<SPComplexType,1>> hij;
  hij.reserve(ndet);
  int skp=((type==COLLINEAR)?1:0);
  for(int n=0, nd=0; n<ndet; ++n, nd+=(skp+1)) {
    check_wavefunction_consistency(type,&PsiT[nd],&PsiT[nd+skp],NMO,NAEA,NAEB);
    hij.emplace_back(rotateHij(type,NMO,NAEA,NAEB,&PsiT[nd],&PsiT[nd+skp],OneBodyHamiltonian::H1));
  }

  // dense one body hamiltonian
  auto H1 = getH1(); 

  if(write_hdf) 
    writeTHCOps(hdf_restart,type,NMO,NAEA,NAEB,ndet,TGprop,TGwfn,H1,
                rotPiu,rotMuv,Piu,Luv,v0,E0);

  return HamiltonianOperations(THCOps<ValueType>(TGwfn.TG_local(),NMO,NAEA,NAEB,type,std::move(H1),
                                      std::move(hij),std::move(rotMuv),std::move(rotPiu),
                                      std::move(rotcPua),std::move(Luv),
                                      std::move(Piu),std::move(cPua),std::move(v0),E0));
}


}
}
