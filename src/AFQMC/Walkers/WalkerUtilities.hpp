//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERUTILITIES_HPP
#define QMCPLUSPLUS_AFQMC_WALKERUTILITIES_HPP
    
#include<cassert>
#include <memory>
#include <mpi.h>
#include<AFQMC/config.0.h>
#include <Utilities/FairDivide.h>

#include "mpi3/communicator.hpp"

namespace qmcplusplus
{

namespace afqmc
{

template<class WlkBucket, 
         class DVec = std::vector<ComplexType>
         >
inline void BasicWalkerData(WlkBucket& wlk, DVec& curData, communicator& comm)
{
  assert(curData.size() >= 7);
  std::fill(curData.begin(),curData.begin()+7,0);
  int nW = wlk.size();
  ComplexType enume = 0, edeno = 0;
  std::vector<double> data(8,0.);
  RealType sumo=0.0;
  int wi=0;
  for(typename WlkBucket::iterator it=wlk.begin(); it!=wlk.end(); ++it) { 
    auto w0 = *it;
    ComplexType weight = w0.weight();
    ComplexType ovlp = w0.overlap(); 
    ComplexType eloc = w0.pseudo_energy();
    data[6]++;   // all walkers
    if( std::abs(weight) <= 1e-6 || (!std::isfinite( std::abs(ovlp) )) || (!std::isfinite( (weight*eloc).real() )) || (!std::isfinite( (weight*eloc).imag() )) ) {
        w0.weight() = ComplexType(0.0,0.0);
        w0.overlap() = ComplexType(1.0,0.0);
        w0.pseudo_energy() = ComplexType(0.0,0.0);
        continue;
    }
    data[0] += (weight*eloc).real();
    data[1] += (weight*eloc).imag();
    data[2] += weight.real();
    data[3] += weight.imag();
    data[4] += std::abs(weight);
    data[5] += std::abs(ovlp);
    data[7]++;   // healthy walkers
  }
  comm.all_reduce_in_place_n(data.begin(),data.size(),std::plus<>()); 
  curData[0] = ComplexType(data[4]/static_cast<RealType>(wlk.get_global_target_population()),0.0);
  curData[1] = ComplexType(data[0]/data[6],data[1]/data[6]);
  curData[2] = ComplexType(data[2]/data[6],data[3]/data[6]);
  curData[3] = data[4];
  curData[4] = data[5]/data[6];
  curData[5] = data[6];
  curData[6] = data[7];
}

template<class WlkBucket,
         class IVec = std::vector<int>
         >
inline void CountWalkers(WlkBucket& wlk, IVec& WCnt, communicator& comm)
{
  WCnt.resize(comm.size(),0);
  int nw = wlk.size();
  comm.all_gather_value(nw,WCnt.begin());
}

template<class WlkBucket,
         typename = typename std::enable_if<(WlkBucket::fixed_population)>::type        
         >
inline void getGlobalListOfWalkerWeights(WlkBucket& wlk, std::vector<std::pair<double,int>>& buffer, communicator& comm)
{
  using Type = std::pair<double,int>;
  int target = wlk.get_TG_target_population();
  if(buffer.size() < target*comm.size())
    APP_ABORT(" Error in getGlobalListOfWalkerWeights(): Array dimensions.\n"); 
  if(wlk.size() > target) 
    APP_ABORT(" Error in getGlobalListOfWalkerWeights(): size > target.\n"); 
  std::vector<Type> blocal(target);
  std::vector<Type>::iterator itv = blocal.begin();
  for(typename WlkBucket::iterator it=wlk.begin(); it!=wlk.end(); ++it, ++itv) 
    *itv = {std::abs(it->weight()),1};
  //comm.gather_n(blocal.data(),blocal.size(),buffer.data(),0);
  MPI_Allgather(blocal.data(), blocal.size()*sizeof(Type), MPI_CHAR, 
                  buffer.data(), blocal.size()*sizeof(Type), MPI_CHAR,
                  &comm);

}

}  // namespace afqmc

} // namespace qmcplusplus

#endif
