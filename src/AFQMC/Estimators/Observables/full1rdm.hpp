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

#ifndef QMCPLUSPLUS_AFQMC_FULL1RDM_HPP
#define QMCPLUSPLUS_AFQMC_FULL1RDM_HPP

#include "AFQMC/config.h"

namespace qmcplusplus
{

namespace afqmc
{

class full1rdm: public AFQMCInfo
{

  // allocators
  using Allocator = device_allocator<ComplexType>;

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;
  using stdCVector = boost::multi::array<ComplexType,1>;
  using stdCMatrix = boost::multi::array<ComplexType,2>;

  public:

  full1rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk):
                                    AFQMCInfo(info),TG(tg_),walker_type(wlk),
                                    writer(false)
  {


    writer = (TG.getGlobalRank()==0);
    std::pair<int,int> dm_dims;
    int nspin=1;
    if(walker_type == CLOSED) {
      dm_dims = {NMO,NMO};
    } else if(walker_type == COLLINEAR) {
      dm_dims = {2*NMO,NMO};
      nspin=2;
    } else if(walker_type == NONCOLLINEAR) {
      dm_dims = {2*NMO,2*NMO};
    }

    if(writer) {
      RDM1Full = std::move(std4CTensor({nave,nspin,dm_dims.first,dm_dims.second}));
      std::fill_n(RDM1Full.begin(),RDM1Full.num_elements(),ComplexType(0.0,0.0));
    }

  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR && Refs, CVec && wgt, MatD && detR, bool impsamp, bool trasposed)
  {

    if( TG.TG_local.root() )
      if(trasposed) {

      } else {

      }
    }
  }



  void print(int iblock, hdf_archive& dump)
  {
    using std::fill_n;
    const int n_zero = 9;

    if( TG.TG_local.root() ) {  
      for(int i=0; i<counters.size(); i++) {
        if(counters[i] > 0) 
          ma::scal(ComplexType(1.0/counters[i]),DMAverage[i]);
      }
      TG.TG_heads.reduce_in_place_n(DMAverage.origin(),DMAverage.num_averages(),std::plus<>(),0);
      if(writer) { 
        for(int i=0; i<counters.size(); ++i) {
          dump.push(std::string("BackProp_")+std::to_string(i));
          std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          stdCVector_ref DMAverage_( DMAverage[i].origin(), {DMAverage.size(1)});
          dump.write(DMAverage_, "full_one_rdm_"+padded_iblock);
          dump.pop();
        }
      } 
      fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
    }
  }

  private:

  // keep count of how many times each average has been called
  std::vector<int> counters;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  bool writer;

  // OneRDMFull (nave, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  stdCMatrix DMAverage;

  // buffer space
  CVector Buff;

  void set_buffer(size_t N) {
    if(Buff.num_elements() < N)
      Buff = std::move(CVector(iextensions<1u>{N}));
    using std::fill_n;
    fill_n(Buff.origin(),N,ComplexType(0.0));
  }

};

#endif
#endif
