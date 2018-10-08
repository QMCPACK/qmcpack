//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/DMC/DMCOMP.h"
#include "Message/OpenMP.h"

#ifdef QMC_CUDA
#include "QMCDrivers/DMC/DMC_CUDA.h"
#endif

//#define PETA_DMC_TEST
namespace qmcplusplus
{

template<>
QMCDriver<Batching::SINGLE>* DMCFactory<Batching::SINGLE>::create(MCWalkerConfiguration& w, TrialWaveFunction<Batching::SINGLE>& psi
                              , QMCHamiltonian& h, HamiltonianPool<Batching::SINGLE>& hpool,WaveFunctionPool& ppool)
{
  app_log() << "Creating DMCOMP for the qmc driver" << std::endl;
  QMCDriver<Batching::SINGLE>*  qmc = new DMCOMP(w,psi,h,ppool);
  qmc->setUpdateMode(PbyPUpdate);
  return qmc;
}

#ifdef QMC_CUDA
template<>
QMCDriver<Batching::BATCHED>* DMCFactory<Batching::BATCHED>::create(MCWalkerConfiguration& w, TrialWaveFunction<Batching::BATCHED>& psi
                              , QMCHamiltonian& h, HamiltonianPool<Batching::BATCHED>& hpool,WaveFunctionPool& ppool)
{
  return new DMCcuda (w, psi, h,ppool);
}
#endif
  
}
