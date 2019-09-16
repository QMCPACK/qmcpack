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

/* 
 * Observable class that calculates the walker averaged "full" 1 RDM.
 * In this context, "full" means that no contraction over the RDM is
 * being performed. The resulting RDM will be [spin][x*NMO][x*NMO],
 * where x:2 for NONCOLLINEAR and 1 for everything else.
 */
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
  using stdCMatrix_ref = boost::multi::array_ref<ComplexType,2>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;
  using mpi3C4Tensor = boost::multi::array<ComplexType,4,shared_allocator<ComplexType>>;

  public:

  full1rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk, 
           int nave_=1, int bsize=1):
                AFQMCInfo(info),TG(tg_),walker_type(wlk),writer(false),
                block_size(bsize),nave(nave_),counter(0),
                hdf_walker_output(""),nskip_walker_output(0),ncnt_walker_output(0),
                denom(iextensions<1u>{0},shared_allocator<ComplexType>{TG.TG_local()}),
                DMWork({0,0},shared_allocator<ComplexType>{TG.TG_local()}),
                DMAverage({0,0},shared_allocator<ComplexType>{TG.TG_local()})
  {

    app_log()<<"  --  Adding Back Propagated Full 1RDM (OneRDM) estimator. -- \n";

    if(cur != NULL) {
      ParameterSet m_param;
      m_param.add(hdf_walker_output, "walker_output", "std::string");
      m_param.add(nskip_walker_output, "nskip_output", "int");
      m_param.put(cur);
    }

    if(hdf_walker_output != std::string("")) { 
      hdf_walker_output = "G"+std::to_string(TG.TG_heads().rank())+"_"+hdf_walker_output;
      hdf_archive dump;  
      if(not dump.create(hdf_walker_output)) {
        app_log()<<"Problems creating walker output hdf5 file: " << hdf_walker_output <<std::endl;
        APP_ABORT("Problems creating walker output hdf5 file.\n");
      }
      dump.push("FullOneRDM");
      dump.push("Metadata");
      dump.write(NMO, "NMO");
      dump.write(NAEA, "NUP");
      dump.write(NAEB, "NDOWN");
      int wlk_t_copy = walker_type; // the actual data type of enum is implementation-defined. convert to int for file
      dump.write(wlk_t_copy, "WalkerType");
      dump.pop();
      dump.pop();
      dump.close();
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

  template<class MatG, class MatG_host, class HostCVec1, class HostCVec2, class HostCVec3>
  void accumulate_reference(int iav, int iref, MatG&& G, MatG_host&& G_host, HostCVec1&& wgt, 
                                               HostCVec2&& Xw, HostCVec3&& ovlp, bool impsamp)
  {
    static_assert(std::decay<MatG>::type::dimensionality == 4, "Wrong dimensionality");
    static_assert(std::decay<MatG_host>::type::dimensionality == 4, "Wrong dimensionality");
    using std::fill_n;
    using std::copy_n;
    // assumes G[nwalk][spin][M][M]
    int nw(G.size(0));
    assert(G.size(0) == wgt.size(0));
    assert(G[0].num_elements() == dm_size);
    assert(wgt.size(0) == nw);
    assert(Xw.size(0) == nw);
    assert(ovlp.size(0) >= nw);
    assert(G.num_elements() == G_host.num_elements());
    assert(G.extensions() == G_host.extensions());
    
    // check structure dimensions
    if(iref == 0) {
      if( denom.size(0) != nw ) {
        denom = std::move(mpi3CVector(iextensions<1u>{nw},
                                      shared_allocator<ComplexType>{TG.TG_local()}));
      }   
      if( DMWork.size(0) != nw ||
          DMWork.size(1) != dm_size ) { 
        DMWork = std::move(mpi3CMatrix({nw,dm_size},
                                      shared_allocator<ComplexType>{TG.TG_local()}));
      }
      fill_n(denom.origin(),denom.num_elements(),ComplexType(0.0,0.0));  
      fill_n(DMWork.origin(),DMWork.num_elements(),ComplexType(0.0,0.0));  
    } else {
      if( denom.size(0) != nw ||
          DMWork.size(0) != nw ||
          DMWork.size(1) != dm_size ||
          DMAverage.size(0) != nave ||  
          DMAverage.size(1) != dm_size ) 
        APP_ABORT(" Error: Invalid state in accumulate_reference. \n\n\n");
    }

    int i0, iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),dm_size,TG.TG_local().size());

    stdCMatrix_ref G2D( to_address(G_host.origin()), {nw, dm_size});
    
    for(int iw=0; iw<nw; iw++) {
      if(TG.TG_local().root()) denom[iw] += Xw[iw];
      ma::axpy( Xw[iw], G2D[iw].sliced(i0,iN), DMWork[iw].sliced(i0,iN) );
    }
    TG.TG_local().barrier();

  }

  template<class HostCVec>
  void accumulate_block(int iav, HostCVec&& wgt, bool impsamp) 
  {
    int nw(denom.size(0));
    int i0,iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),dm_size,TG.TG_local().size());

    if(iav == 0) ncnt_walker_output++; 
    if(hdf_walker_output != std::string("") && ncnt_walker_output%(nskip_walker_output+1)==0) { 
      const int n_zero = 9;
      for(int iw=0; iw<nw; iw++) 
        ma::scal(ComplexType(1.0,0.0)/denom[iw], DMWork[iw].sliced(i0,iN));
      TG.TG_local().barrier();

      if(TG.TG_local().root()) {
        hdf_archive dump;
        if(iav==0) counter++;
        if(not dump.open(hdf_walker_output)) {
          app_log()<<"Problems opening walker output hdf5 file: " 
                   << hdf_walker_output <<std::endl;
          APP_ABORT("Problems opening walker output hdf5 file.\n");
        }
        dump.push("FullOneRDM");
        dump.push(std::string("Group")+std::to_string(TG.TG_heads().rank()));
        dump.push(std::string("Average_")+std::to_string(iav));
        std::string padded_num = std::string(n_zero-std::to_string(counter).length(),'0')+
                                    std::to_string(counter); 
        dump.write(wgt,"weights_"+padded_num);
        stdCMatrix_ref DM(to_address(DMWork.origin()),{nw,dm_size}); 
        dump.write(DM,"one_rdm_"+padded_num);
        dump.pop();
        dump.pop();
        dump.pop();
        dump.close();

        // adjust denom
        for(int iw=0; iw<nw; iw++) 
          denom[iw] = wgt[iw]; 
      }
      TG.TG_local().barrier();
    } else {
      if(TG.TG_local().root()) 
        for(int iw=0; iw<nw; iw++) 
          denom[iw] = wgt[iw]/denom[iw];
    }
    TG.TG_local().barrier();

    // DMAverage[iav][ij] += sum_iw DMWork[iw][ij] * denom[iw] = T( DMWork ) * denom
    ma::product( ComplexType(1.0,0.0), ma::T( DMWork( {0, nw}, {i0,iN}) ),  denom, 
                 ComplexType(1.0,0.0), DMAverage[iav].sliced(i0,iN)); 
    TG.TG_local().barrier();
  }

  template< class HostCVec>
  void print(int iblock, hdf_archive& dump, HostCVec&& Wsum)
  {
    using std::fill_n;
    const int n_zero = 9;

    if( TG.TG_local().root() ) {  
      ma::scal(ComplexType(1.0/block_size),DMAverage);
      TG.TG_heads().reduce_in_place_n(to_address(DMAverage.origin()),DMAverage.num_elements(),std::plus<>(),0);
      if(writer) { 
        dump.push(std::string("FullOneRDM"));
        for(int i=0; i<nave; ++i) {
          dump.push(std::string("Average_")+std::to_string(i));
          std::string padded_iblock = 
                std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          stdCVector_ref DMAverage_( to_address(DMAverage[i].origin()), {dm_size});
          dump.write(DMAverage_, "one_rdm_"+padded_iblock);
          dump.write(Wsum[i], "denominator_"+padded_iblock);
          dump.pop();
        }
        dump.pop();
      } 
    }
    TG.TG_local().barrier();
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
  }

  private:

  int block_size;  

  int nave;

  int counter;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  int dm_size;

  bool writer;

  std::string hdf_walker_output;  

  int nskip_walker_output;

  int ncnt_walker_output;

  // DMAverage (nave, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMAverage;
  // DMWork (nwalk, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMWork;

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
