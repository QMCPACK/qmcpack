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

  FullObsHandler(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string name_, xmlNodePtr cur, WALKER_TYPES wlk, 
        Wavefunction& wfn):
                                    AFQMCInfo(info),TG(tg_),name(name_),walker_type(wlk),
                                    wfn0(wfn), writer(false)
  {

    int nave=0;
    std::string obs("");
    if(cur != NULL) {
      ParameterSet m_param;
      m_param.add(nave, "naverages", "int");
      m_param.add(obs, "observables", "std::string");
      m_param.put(cur);
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    counters.resize(nave);

    if(obs == std::string(""))
      APP_ABORT("empty observables list is not allowed.\n");

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

    writer = (TG.getGlobalRank()==0);

    denominator = = std::move(stdCVector(iextensions<1u>{nave})); 
    std::fill_n(denominator.begin(),denominator.num_elements(),ComplexType(0.0,0.0));  

  }

  void print(int iblock, hdf_archive& dump)
  {
    const int n_zero = 9;

    if(writer) {
      for(int i=0; i<counters.size(); i++) {
        if(counters[i] > 0)
          denominator[i] /= counters[i];
      }
      for(int i=0; i<counters.size(); ++i) {
        dump.push(std::string("BackProp_")+std::to_string(i));
        std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
        stdCVector_ref denom( denominator[i].origin(), {1});
        dump.write(denom, "denominator_"+padded_iblock);
        dump.pop();
      }
    }

    for(auto& v: properties) v.print(iblock,dump);

    using std::fill_n;
    fill_n(denominator.origin(), denominator.num_elements(), ComplexType(0.0,0.0));
    fill_n(counters.begin(),counters.end(),0);
  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR && Refs, CVec && wgt, MatD && detR, bool impsamp)
  {
    if(iav < 0 || iav >= counters.size())
      APP_ABORT("Runtime Error: iav out of range in full1rdm::accumulate. \n\n\n");
    counters[iav]++;

    //1. Calculate Green functions

    //2. accumulate observables

  }

  private:

  std::vector<int> counters;

  std::string name;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  Wavefunction& wfn0;

  std::vector<Observable> properties;

  bool writer;

  // denominator (nave, ...)  
  stdCVector denominator;    

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
