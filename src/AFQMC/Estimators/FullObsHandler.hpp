//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_FULLOBSHANDLER_HPP
#define QMCPLUSPLUS_AFQMC_FULLOBSHANDLER_HPP

#include <vector>
#include <string>
#include <iostream>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Estimators/Observables/Observable.hpp"
#include "AFQMC/config.h"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"


namespace qmcplusplus
{

namespace afqmc
{

/*
 * This class manages a list of "full" observables.
 * Full observables are those that have a walker dependent left-hand side, 
 * which result from back propagation. 
 * This implementation of the class assumes a multi-determinant trial wavefunction,  
 * resulting in the loop over "references" (each determinant in the trial wavefunciton
 * being back-propagated). 
 * Given a walker set and an array of (back propagated) slater matrices, 
 * this routine will calculate and accumulate all requested observables.
 * To make the implementation of the BackPropagated class cleaner, 
 * this class also handles all the hdf5 I/O (given a hdf archive).
 */
class FullObsHandler: public AFQMCInfo
{

  // allocators
  using sharedAllocator = localTG_allocator<ComplexType>;

  using shared_pointer = typename sharedAllocator::pointer;
  using const_shared_pointer = typename sharedAllocator::const_pointer;

  using devCMatrix_ref = boost::multi::array_ref<ComplexType,2,device_ptr<ComplexType>>;

  using sharedCVector = boost::multi::array<ComplexType,1,sharedAllocator>;
  using sharedCVector_ref = boost::multi::array_ref<ComplexType,1,shared_pointer>;
  using sharedCMatrix_ref = boost::multi::array_ref<ComplexType,2,shared_pointer>;
  using sharedC4Tensor_ref = boost::multi::array_ref<ComplexType,4,shared_pointer>;

  using stdCVector = boost::multi::array<ComplexType,1>;
  using stdCMatrix = boost::multi::array<ComplexType,2>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;

  public:

  FullObsHandler(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string name_, xmlNodePtr cur, WALKER_TYPES wlk, 
        Wavefunction& wfn):
                                    AFQMCInfo(info),TG(tg_),name(name_),walker_type(wlk),
                                    wfn0(wfn), writer(false), block_size(1), nave(1),
                                    nspins((walker_type==COLLINEAR)?2:1),
                                    Buff(iextensions<1u>{1},make_localTG_allocator<ComplexType>(TG))
  {

    using std::fill_n;

    xmlNodePtr curRoot = cur;
    if(cur != NULL) {
      ParameterSet m_param;
      m_param.add(nave, "naverages", "int");
      m_param.add(block_size, "block_size", "int");
      m_param.put(cur);
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="OneRDM") {
        properties.emplace_back(Observable(std::move(full1rdm(TG,info,cur,walker_type,nave,block_size))));
      }
      cur = cur->next;
    }

    if(properties.size() == 0)
      APP_ABORT("empty observables list is not allowed.\n");

    Gdims = std::make_tuple(NMO,NMO);
    if(walker_type == NONCOLLINEAR)
      Gdims = std::make_tuple(2*NMO,2*NMO);
    dm_size = nspins * std::get<0>(Gdims) * std::get<1>(Gdims); 

    writer = (TG.getGlobalRank()==0);

    denominator = std::move(stdCVector(iextensions<1u>{nave})); 
    fill_n(denominator.begin(),denominator.num_elements(),ComplexType(0.0,0.0));  

  }

  void print(int iblock, hdf_archive& dump)
  {
    using std::fill_n;
    const int n_zero = 9;

    if( TG.TG_local().root() ) {
      ma::scal(ComplexType(1.0/block_size),denominator);
      TG.TG_heads().reduce_in_place_n(to_address(denominator.origin()),denominator.num_elements(),std::plus<>(),0);
    }

    for(auto& v: properties) v.print(iblock,dump,denominator);
    fill_n(denominator.origin(), denominator.num_elements(), ComplexType(0.0,0.0));
  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR && Refs, CVec && wgt, MatD && DevdetR, bool impsamp)
  {
    if(iav < 0 || iav >= nave)
      APP_ABORT("Runtime Error: iav out of range in full1rdm::accumulate. \n\n\n");

    int nw(wset.size());
    int nrefs(Refs.size(1));
    double LogOverlapFactor(wset.getLogOverlapFactor());
    set_buffer( nw * (dm_size+3) );
    sharedC4Tensor_ref G4D(Buff.origin(), {nw, nspins, std::get<0>(Gdims),std::get<1>(Gdims)});
    sharedCMatrix_ref G2D(Buff.origin(), {nw, dm_size});
    sharedCVector_ref DevOv(G4D.origin()+G4D.num_elements(), {2*nw});

    stdCVector Xw(iextensions<1u>{nw});
    stdCVector Ov(iextensions<1u>{2*nw});
    stdCMatrix detR(DevdetR); 

    using SMType = typename WlkSet::reference::SMType;
    // MAM: The pointer type of GA/GB needs to be device_ptr, it can not be  
    //      one of the shared_memory types. The dispatching in DensityMatrices is done
    //      through the pointer type of the result matrix (GA/GB).
    std::vector<devCMatrix_ref> GA;
    std::vector<devCMatrix_ref> GB;
    std::vector<SMType> RefsA;
    std::vector<SMType> RefsB;
    std::vector<SMType> SMA;
    std::vector<SMType> SMB;
    RefsA.reserve(nw);
    SMA.reserve(nw);
    GA.reserve(nw);
    if(walker_type == COLLINEAR) RefsB.reserve(nw);
    if(walker_type == COLLINEAR) SMB.reserve(nw);
    if(walker_type == COLLINEAR) GB.reserve(nw);


    if(impsamp) 
      denominator[iav] += std::accumulate(wgt.begin(),wgt.end(),ComplexType(0.0));
    else {
      APP_ABORT(" Finish implementation of free projection. \n\n\n");
    }

    for(int iref=0, is=0; iref<nrefs; iref++, is+=nspins) {

      // conjugated here!
      ComplexType CIcoeff( std::conj(wfn0.getReferenceWeight(iref)) );

      //1. Calculate Green functions
      // Refs({wset.size(),nrefs,ref_size}  
      RefsA.clear();  
      RefsB.clear();  
      SMA.clear();  
      SMB.clear();  
      GA.clear();  
      GB.clear();  
      // using SlaterMatrixAux to store References in device memory
      if(walker_type == COLLINEAR) {
        for(int iw=0; iw<nw; iw++) {
          SMA.emplace_back(wset[iw].SlaterMatrixN(Alpha));
          SMB.emplace_back(wset[iw].SlaterMatrixN(Beta));
          GA.emplace_back( devCMatrix_ref(make_device_ptr(G2D[iw].origin()),{NMO,NMO}) );
          GB.emplace_back( devCMatrix_ref(make_device_ptr(G2D[iw].origin())+NMO*NMO,{NMO,NMO}) );
          RefsA.emplace_back(wset[iw].SlaterMatrixAux(Alpha));
          RefsB.emplace_back(wset[iw].SlaterMatrixAux(Beta));
          copy_n(Refs[iw][iref].origin() , RefsA.back().num_elements(), RefsA.back().origin());
          copy_n(Refs[iw][iref].origin()+RefsA.back().num_elements() , 
                 RefsB.back().num_elements() , RefsB.back().origin());
        }
        wfn0.DensityMatrix(RefsA, SMA, GA, DevOv.sliced(0,nw), LogOverlapFactor, false, false);
        wfn0.DensityMatrix(RefsB, SMB, GB, DevOv.sliced(nw,2*nw), LogOverlapFactor, false, false);
      } else {
        for(int iw=0; iw<nw; iw++) {
          SMA.emplace_back(wset[iw].SlaterMatrixN(Alpha));
          GA.emplace_back( devCMatrix_ref(make_device_ptr(G2D[iw].origin()),{NMO,NMO}) );
          RefsA.emplace_back(wset[iw].SlaterMatrixAux(Alpha));
          copy_n(Refs[iw][iref].origin() , RefsA.back().num_elements(), RefsA.back().origin());
        }
        wfn0.DensityMatrix(RefsA, SMA, GA, DevOv.sliced(0,nw), LogOverlapFactor, false, false);
      } 

      //2. calculate and accumulate appropriate weights 
      copy_n( DevOv.origin(), 2*nw, Ov.origin());
      if(walker_type == CLOSED) { 
        for(int iw=0; iw<nw; iw++) 
          Xw[iw] = CIcoeff * Ov[iw] * detR[iw][iref] * detR[iw][iref]; 
      } else if(walker_type == COLLINEAR) {
        for(int iw=0; iw<nw; iw++) { 
          Xw[iw] = CIcoeff * Ov[iw] * Ov[iw+nw] * detR[iw][2*iref] * detR[iw][2*iref+1]; 
        }
      } else if(walker_type == NONCOLLINEAR) {
        for(int iw=0; iw<nw; iw++) 
          Xw[iw] = CIcoeff * Ov[iw] * detR[iw][iref]; 
      } 

      //3. accumulate references 
      for(auto& v: properties) v.accumulate_reference(iav,iref,G4D,wgt,Xw,Ov,impsamp);

    }
    //4. accumulate block (normalize and accumulate sum over references)
    for(auto& v: properties) v.accumulate_block(iav,wgt,impsamp);

  }

  private:

  int nave;

  int block_size;

  std::string name;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  int nspins;
  int dm_size;
  std::tuple<int,int> Gdims;

  Wavefunction& wfn0;

  std::vector<Observable> properties;

  bool writer;

  // denominator (nave, ...)  
  stdCVector denominator;    

  // buffer space
  sharedCVector Buff;

  void set_buffer(size_t N) {
    if(Buff.num_elements() < N)
      Buff = std::move(sharedCVector(iextensions<1u>{N},make_localTG_allocator<ComplexType>(TG)));
    using std::fill_n;
    fill_n(Buff.origin(),N,ComplexType(0.0));
  }

};

}

}

#endif
