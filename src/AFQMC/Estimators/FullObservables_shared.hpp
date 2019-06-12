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

#ifndef QMCPLUSPLUS_AFQMC_FULLOBSERVABLE_SHARED_HPP
#define QMCPLUSPLUS_AFQMC_FULLOBSERVABLE_SHARED_HPP

#include <vector>
#include <string>
#include <iostream>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Estimators/Observables_config.h"
#include "AFQMC/config.h"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"


namespace qmcplusplus
{

namespace afqmc
{

class FullObservables_shared: public AFQMCInfo
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
  using stdCTensor = boost::multi::array<ComplexType,3>;
  using std4CTensor = boost::multi::array<ComplexType,4>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;

  public:

  FullObservables_shared(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string name_, xmlNodePtr cur, WALKER_TYPES wlk, 
        Wavefunction& wfn):
                                    AFQMCInfo(info),TG(tg_),name(name_),walker_type(wlk),
                                    wfn0(wfn), writer(false)
  {
    std::fill(measure.begin(),measure.end(),false);

    std::string obs("");
    if(cur != NULL) {
      ParameterSet m_param;
      m_param.add(nave, "naverages", "int");
      m_param.add(obs, "observables", "std::string");
      m_param.put(cur);
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    if(obs == std::string(""))
      APP_ABORT("empty observables list is not allowed.\n");

    if(obs.find("_1rdm_") != std::string::npos) measure[OneRDMFull]=true; 
    if(obs.find("_2rdm_") != std::string::npos) measure[TwoRDMFull]=true; 
    if(obs.find("_ekt_") != std::string::npos) measure[GFockOpa]=true; 
    if(obs.find("_ekt_") != std::string::npos) measure[GFockOpb]=true; 
    // need to read rotation matrix
    if(obs.find("_1rdmc_") != std::string::npos) measure[OneRDMc]=true; 
    if(obs.find("_2rdmc_") != std::string::npos) measure[TwoRDMc]=true; 

    if(std::find(measure.begin(),measure.end(),true) == measure.end())
      APP_ABORT("empty observables list is not allowed.\n");

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

    denominator = = std::move(stdCMatrix({nave,1})); 
    std::fill_n(denominator.begin(),denominator.num_elements(),ComplexType(0.0,0.0));  

    // Full 1-RDM
    if(measure[OneRDMFull]) {
      if(writer) { 
        RDM1Full = std::move(std4CTensor({nave,nspin,dm_dims.first,dm_dims.second})); 
        std::fill_n(RDM1Full.begin(),RDM1Full.num_elements(),ComplexType(0.0,0.0));  
      }
    }

    // ...

  }

  void print(hdf_archive& dump)
  {
    if(!writer) return;

    

    ma::scal(ComplexType(1.0/block_size),denom_average);

    if(measure[OneRDMFull]) {
      ma::axpy(ComplexType(1.0),buff,DMAverage);

    } 
      if(iblock%block_size == 0) {
//        for(int i = 0; i < DMAverage.size(); i++)
//          DMAverage[i] /= block_size;
//        denom_average[0] /= block_size;
        ma::scal(ComplexType(1.0/block_size),DMAverage);
        dump.push("BackPropagated");
        for(int i=0; i<nave; ++i) {
          dump.push(std::string("NumBackProp_")+std::to_string(nback_prop_steps[i]));
          std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          stdCVector_ref DMAverage_( DMAverage[i].origin(), {DMAverage.size(1)});
          stdCVector_ref denom_average_( denom_average[i].origin(), {denom_average.size(1)});
          stdCVector_ref wOvlp_( wOvlp[i].origin(), {wOvlp.size(1)*wOvlp.size(2)});
          stdCVector_ref wDMsum_( wDMsum[i].origin(), {wDMsum.size(1)*wDMsum.size(2)});
          stdCVector_ref wdetR_( wdetR[i].origin(), {wdetR.size(1)*wdetR.size(2)});
          dump.write(DMAverage_, "one_rdm_"+padded_iblock);
          dump.write(denom_average_, "one_rdm_denom_"+padded_iblock);
          dump.write(wOvlp_, "one_rdm_walker_overlaps_"+padded_iblock);
          dump.write(wDMsum_, "one_rdm_walker_dm_sums_"+padded_iblock);
          dump.write(wdetR_, "one_rdm_detR_"+padded_iblock);
          dump.pop();
        }
        dump.pop();
        using std::fill_n;
        fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
        fill_n(denom_average.origin(), denom_average.num_elements(), ComplexType(0.0,0.0));
        fill_n(wdetR.origin(),wdetR.num_elements(),ComplexType(0.0));
        fill_n(wOvlp.origin(),wOvlp.num_elements(),ComplexType(0.0));
        fill_n(wDMsum.origin(),wDMsum.num_elements(),ComplexType(0.0));
      }
  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR && Refs, CVec && wgt, MatD && detR, bool impsamp)
  {

  }



  private:

  int nave;

  std::string name;

  std::array<bool,10> measure;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  Wavefunction& wfn0;

  bool writer;

  // denominator (nave, ...)  
  stdCMatrix denominator;    

  /**************************************************
   *                 Accumulators                   * 
   **************************************************/

  // OneRDMFull (nave, spin, x*NMO, x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  std4CTensor RDM1Full; 



  /**************************************************/

  // buffer space
  CVector Buff;

  void set_buffer(size_t N) {
    if(Buff.num_elements() < N)
      Buff = std::move(CVector(iextensions<1u>{N}));
    using std::fill_n;
    fill_n(Buff.origin(),N,ComplexType(0.0));
  }

};

}

}

#endif
