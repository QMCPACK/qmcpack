//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author: jnkim $
 * $Revision: 1.3 $   $Date: 2006/04/05 00:49:59 $
 * $Id: DMCFactory.cpp,v 1.3 2006/04/05 00:49:59 jnkim Exp $
 ***************************************************************************/
