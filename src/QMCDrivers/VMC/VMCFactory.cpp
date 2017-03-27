//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#if defined(QMC_BUILD_COMPLETE)
//REMOVE Broken warping
//#if !defined(QMC_COMPLEX)
//#include "QMCDrivers/VMC/VMCMultipleWarp.h"
//#include "QMCDrivers/VMC/VMCPbyPMultiWarp.h"
//#endif
//#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#endif
#include "Message/OpenMP.h"

#ifdef QMC_CUDA
#include "QMCDrivers/VMC/VMC_CUDA.h"
#endif

namespace qmcplusplus
{

QMCDriver* VMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                              QMCHamiltonian& h, ParticleSetPool& ptclpool, HamiltonianPool& hpool, WaveFunctionPool& ppool)
{
  int np=omp_get_max_threads();
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  QMCDriver* qmc=0;
#ifdef QMC_CUDA
  if (VMCMode & 16)
    qmc = new VMCcuda(w,psi,h,ppool);
  else
#endif
    if(VMCMode == 0 || VMCMode == 1) //(0,0,0) (0,0,1)
    {
      qmc = new VMCSingleOMP(w,psi,h,hpool,ppool);
    }
  //else if(VMCMode == 2) //(0,1,0)
  //{
  //  qmc = new VMCMultiple(w,psi,h);
  //}
  //else if(VMCMode == 3) //(0,1,1)
  //{
  //  qmc = new VMCPbyPMultiple(w,psi,h);
  //}
    else if(VMCMode ==2 || VMCMode ==3)
    {
      qmc = new CSVMC(w,psi,h,hpool,ppool);
    }
//#if !defined(QMC_COMPLEX)
//    else if(VMCMode == 6) //(1,1,0)
//    {
//      qmc = new VMCMultipleWarp(w,psi,h, ptclpool);
//    }
//    else if(VMCMode == 7) //(1,1,1)
//    {
//      qmc = new VMCPbyPMultiWarp(w,psi,h, ptclpool);
//    }
//#endif
//     else if(VMCMode == 8) //(only possible for WFMC run)
//     {
//       qmc = new WFMCSingleOMP(w,psi,h,hpool,ppool);
//     }
  qmc->setUpdateMode(VMCMode&1);
  return qmc;
}
}
