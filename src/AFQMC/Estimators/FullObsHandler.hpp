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
                                    wfn0(wfn), writer(false), block_size(1), nave(0),
                                    nspins((walker_type==COLLINEAR)?2:1),
                                    Buff(iextensions<1u>{1},make_localTG_allocator<ComplexType>(TG))
  {

    using std::fill_n;
    std::string obs("");
    if(cur != NULL) {
      ParameterSet m_param;
      m_param.add(nave, "naverages", "int");
      m_param.add(obs, "observables", "std::string");
      m_param.add(block_size, "block_size", "int");
      m_param.put(cur);
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    if(obs == std::string(""))
      APP_ABORT("empty observables list is not allowed.\n");

    // add _XXX_
    obs = std::string("_") + obs + std::string("_");

    if(obs.find("_1rdm_") != std::string::npos) 
      properties.emplace_back(Observable(std::move(full1rdm(TG,info,cur,walker_type))));
/*
    if(obs.find("_2rdm_") != std::string::npos) measure[TwoRDMFull]=true; 
    if(obs.find("_ekt_") != std::string::npos) measure[GFockOpa]=true; 
    if(obs.find("_ekt_") != std::string::npos) measure[GFockOpb]=true; 
    // need to read rotation matrix
    if(obs.find("_1rdmc_") != std::string::npos) measure[OneRDMc]=true; 
    if(obs.find("_2rdmc_") != std::string::npos) measure[TwoRDMc]=true; 
*/

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

    if(writer) {
      ma::scal(ComplexType(1.0/block_size),denominator);
      for(int i=0; i<nave; ++i) {
        dump.push(std::string("BackProp_")+std::to_string(i));
        std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
        stdCVector_ref denom( denominator.origin()+i, {1});
        dump.write(denom, "denominator_"+padded_iblock);
        dump.pop();
      }
    }

    for(auto& v: properties) v.print(iblock,dump);
    fill_n(denominator.origin(), denominator.num_elements(), ComplexType(0.0,0.0));
  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR && Refs, CVec && wgt, MatD && DevdetR, bool impsamp)
  {
    if(iav < 0 || iav >= nave)
      APP_ABORT("Runtime Error: iav out of range in full1rdm::accumulate. \n\n\n");

    int nw(wset.size());
    int nrefs(Refs.size(1));
    set_buffer( nw * (dm_size+3) );
    sharedC4Tensor_ref G4D(Buff.origin(), {nw, nspins, std::get<0>(Gdims),std::get<1>(Gdims)});
    sharedCMatrix_ref G2D(Buff.origin(), {nw, dm_size});
    sharedCVector_ref DevOv(G4D.origin()+G4D.num_elements(), {nw});

    stdCVector Xw(iextensions<1u>{nw});
    stdCVector Ov(iextensions<1u>{nw});
    stdCMatrix detR(DevdetR); 
    std::vector<devCMatrix_ref> RefsA;
    std::vector<devCMatrix_ref> RefsB;
    RefsA.reserve(nw);
    if(walker_type == COLLINEAR) RefsB.reserve(nw);

    for(int iref=0, is=0; iref<nrefs; iref++, is+=nspins) {

      ComplexType CIcoeff(wfn0.getReferenceWeight(iref));

      //1. Calculate Green functions
      // Refs({wset.size(),nrefs,ref_size}  
      RefsA.clear();  
      RefsB.clear();  
      if(walker_type == CLOSED) {
        for(int iw=0; iw<nw; iw++) 
          RefsA.emplace_back(devCMatrix_ref(make_device_ptr(Refs[iw][iref].origin()), 
                                            {NAEA,NMO}));
      } else if(walker_type == COLLINEAR) {
        for(int iw=0; iw<nw; iw++) { 
          RefsA.emplace_back(devCMatrix_ref(make_device_ptr(Refs[iw][iref].origin()), 
                                            {NAEA,NMO}));
          RefsB.emplace_back(devCMatrix_ref(make_device_ptr(Refs[iw][iref].origin()), 
                                            {NAEB,NMO}));
        }
      } else if(walker_type == NONCOLLINEAR) {
        for(int iw=0; iw<nw; iw++) 
          RefsA.emplace_back(devCMatrix_ref(make_device_ptr(Refs[iw][iref].origin()), 
                                            {NAEA+NAEB,2*NMO}));
      }
      wfn0.DensityMatrix(wset, RefsA, RefsB, G2D, DevOv, false, false, true);

      //2. calculate and accumulate appropriate weights 
      copy_n( make_device_ptr(DevOv.origin()), nw, Ov.origin());
      if(walker_type == CLOSED) { 
        for(int iw=0; iw<nw; iw++) 
          Xw[iw] = CIcoeff * Ov[iw] * detR[iw][iref] * detR[iw][iref]; 
      } else if(walker_type == COLLINEAR) {
        for(int iw=0; iw<nw; iw++) 
          Xw[iw] = CIcoeff * Ov[iw] * detR[iw][2*iref] * detR[iw][2*iref+1]; 
      } else if(walker_type == NONCOLLINEAR) {
        for(int iw=0; iw<nw; iw++) 
          Xw[iw] = CIcoeff * Ov[iw] * detR[iw][iref]; 
      } 

      //3. accumulate references 
      for(auto& v: properties) v.accumulate_reference(iav,iref,G4D,wgt,Xw,Ov);

    }
    //4. accumulate block (normalize and accumulate sum over references)
    for(auto& v: properties) v.accumulate_block(iav);

  }

  private:

  int nave;

  int block_size;

  // true: G[nwalk][dm_size], false: G[dm_size][nwalk]
  bool transposedG;

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
