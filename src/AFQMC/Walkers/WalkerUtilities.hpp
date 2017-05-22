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
    
    
#include<cassert>
#include <memory>
#include <mpi.h>
#include<AFQMC/config.0.h>
#include <Utilities/UtilityFunctions.h>

namespace qmcplusplus
{

namespace afqmc
{

template<class WlkBucket, 
         class DVec = std::vector<ComplexType>
         >
// eventually generalize MPI_Comm to a MPI wrapper
inline void BasicWalkerData(WlkBucket& wlk, DVec& curData, MPI_Comm comm)
{
  assert(curData.size() >= 7);
  std::fill(curData.begin(),curData.begin()+7,0);
  int nW = wlk.numWalkers(true);
  ComplexType enume = 0, edeno = 0;
  std::vector<double> data(16,0);
  ComplexType w,oa,ob,eloc;
  RealType sumo=0.0;
  for(int i=0; i<nW; i++) {
    ComplexType* dum = wlk.getWalker(i,w,eloc,oa,ob);
    if( !wlk.isAlive(i) ) continue; 
    data[6]++;   // all walkers
    //if( std::abs(w) <= 1e-6 || std::abs(oa*ob)<1e-8 || (!std::isfinite( std::abs(oa*ob) )) || (!std::isfinite( (w*eloc).real() )) ) continue;
    if( std::abs(w) <= 1e-6 || (!std::isfinite( std::abs(oa*ob) )) || (!std::isfinite( (w*eloc).real() )) ) continue;
    data[0] += (w*eloc).real();
    data[1] += (w*eloc).imag();
    data[2] += w.real();
    data[3] += w.imag();
    data[4] += std::abs(w);
    data[5] += std::abs(oa*ob);
    data[7]++;   // healthy walkers
  }
  
  MPI_Allreduce(data.data(),data.data()+8,8,MPI_DOUBLE,MPI_SUM,comm);
  curData[0] = ComplexType(data[12]/static_cast<RealType>(wlk.get_global_target_population()),0.0);
  curData[1] = ComplexType(data[8]/data[14],data[9]/data[14]);
  curData[2] = ComplexType(data[10]/data[14],data[11]/data[14]);
  curData[3] = data[12];
  curData[4] = data[13]/data[14];
  curData[5] = data[14];
  curData[6] = data[15];

}

template<class WlkBucket,
         class IVec = std::vector<int>
         >
inline void CountWalkers(WlkBucket& wlk, IVec& WCnt, MPI_Comm comm)
{
  int npr, rk;
  MPI_Comm_size(comm,&npr);
  MPI_Comm_rank(comm,&rk);
  WCnt.resize(npr);

  std::fill(WCnt.begin(), WCnt.end(), 0);
  int nw = wlk.numWalkers(false);
  MPI_Allgather(&nw, 1, MPI_INT, WCnt.data(), 1, MPI_INT, comm);
}



}  // namespace afqmc

} // namespace qmcplusplus
