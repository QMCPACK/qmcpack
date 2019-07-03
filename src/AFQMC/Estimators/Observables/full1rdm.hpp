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
#include <vector>
#include <string>
#include <iostream>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"
#include "Utilities/Timer.h"

#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

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

  using CVector_ref = boost::multi::array_ref<ComplexType,1,pointer>;
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;
  using mpi3C4Tensor = boost::multi::array<ComplexType,4,shared_allocator<ComplexType>>;

  public:

  full1rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk):
                AFQMCInfo(info),TG(tg_),walker_type(wlk),writer(false),block_size(1),nave(1),
                denom(iextensions<1u>{0},shared_allocator<ComplexType>{TG.TG_local()}),
                DMWork({0,0,0},shared_allocator<ComplexType>{TG.TG_local()}),
                DMAverage({0,0},shared_allocator<ComplexType>{TG.TG_local()})
  {
    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nave, "naverages", "int");
      m_param.add(block_size, "block_size", "int");
      m_param.put(cur);
    }

    using std::fill_n;
    writer = (TG.getGlobalRank()==0);
    dm_size=NMO*NMO;
    if(walker_type == COLLINEAR) 
      dm_size *= 2; 
    else if(walker_type == NONCOLLINEAR) 
      dm_size *= 2; 

    DMAverage = std::move(mpi3CMatrix({nave,dm_size},shared_allocator<ComplexType>{TG.TG_local()}));
    fill_n(DMAverage.origin(),DMAverage.num_elements(),ComplexType(0.0,0.0));

  }

  template<class MatG, class HostCVec1, class HostCVec2, class HostCVec3>
  void accumulate_reference(int iav, int iref, MatG&& G, HostCVec1&& wgt, 
                                               HostCVec2&& Xw, HostCVec3&& ovlp)
  {
    static_assert(std::decay<MatG>::type::dimensionality == 4, "Wrong dimensionality");
    using std::fill_n;
    using std::copy_n;
    // assumes G[nwalk][spin][M][M]
    int nw(G.size(0));
    assert(G.size(0) == wgt.size(0));
    assert(G[0].num_elements() == dm_size);
    assert(wgt.size(0) == nw);
    assert(Xw.size(0) == nw);
    assert(ovlp.size(0) == nw);
    
    // check structure dimensions
    if(iref == 0) {
      if( denom.size(0) != nw ) {
        denom = std::move(mpi3CVector(iextensions<1u>{nw},
                                      shared_allocator<ComplexType>{TG.TG_local()}));
      }   
      if( DMWork.size(0) != 2 ||
          DMWork.size(1) != nw ||
          DMWork.size(2) != dm_size ) { 
        DMWork = std::move(mpi3CTensor({2,nw,dm_size},
                                      shared_allocator<ComplexType>{TG.TG_local()}));
      }
      fill_n(denom.origin(),nw,ComplexType(0.0,0.0));  
      fill_n(DMWork.origin(),nw,ComplexType(0.0,0.0));  
    } else {
      if( denom.size(0) != nw ||
          DMWork.size(0) != 2 || 
          DMWork.size(1) != nw ||
          DMWork.size(2) != dm_size ||
          DMAverage.size(0) != nave ||  
          DMAverage.size(1) != dm_size ) 
        APP_ABORT(" Error: Invalid state in accumulate_reference. \n\n\n");
    }

    int i0,iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),int(G.num_elements()),
                                         TG.TG_local().size());
    copy_n( make_device_ptr(G.origin())+i0, iN-i0, to_address(DMWork[0].origin())+i0);
    TG.TG_local().barrier();

    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),dm_size,TG.TG_local().size());
    
    for(int iw=0; iw<nw; iw++) {
      if(TG.TG_local().root()) denom[iw] += Xw[iw];
      ma::axpy( Xw[iw]*wgt[iw], DMWork[0][iw].sliced(i0,iN), DMWork[1][iw].sliced(i0,iN) );
    }
  }

  void accumulate_block(int iav) 
  {
    int nw(denom.size(0));
    int i0,iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),dm_size,TG.TG_local().size());
    if(TG.TG_local().root()) 
      for(int iw=0; iw<nw; iw++) 
        denom[iw] = ComplexType(1.0,0.0)/denom[iw]; 
    TG.TG_local().barrier();

    // DMAverage[iav][ij] = sum_iw DMWork[1][iw][ij] * denom[iw] = T( DMWork[1] ) * denom
    ma::product( ma::T( DMWork[1]( {0, nw}, {i0,iN}) ),  denom, DMAverage[iav].sliced(i0,iN)); 
    TG.TG_local().barrier();
  }

  void print(int iblock, hdf_archive& dump)
  {
    using std::fill_n;
    const int n_zero = 9;

    if( TG.TG_local().root() ) {  
      ma::scal(ComplexType(1.0/block_size),DMAverage);
      TG.TG_heads().reduce_in_place_n(to_address(DMAverage.origin()),DMAverage.num_elements(),std::plus<>(),0);
      if(writer) { 
        for(int i=0; i<nave; ++i) {
          dump.push(std::string("BackProp_")+std::to_string(i));
          std::string padded_iblock = 
                std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          stdCVector_ref DMAverage_( to_address(DMAverage[i].origin()), {dm_size});
          dump.write(DMAverage_, "full_one_rdm_"+padded_iblock);
          dump.pop();
        }
      } 
    }
    TG.TG_local().barrier();
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
  }

  private:

  int nave;

  int block_size;  

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  int dm_size;

  bool writer;

  // DMAverage (nave, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMAverage;
  // DMWork (k, nwalk, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  // k:0 communication buffer
  //   1 accumulate over references 
  mpi3CTensor DMWork;

  mpi3CVector denom; 

  // buffer space
  CVector Buff;

  void set_buffer(size_t N) {
    if(Buff.num_elements() < N)
      Buff = std::move(CVector(iextensions<1u>{N}));
    using std::fill_n;
    fill_n(Buff.origin(),N,ComplexType(0.0));
  }

};

} // namespace afqmc
} // namespace qmcplusplus

#endif
