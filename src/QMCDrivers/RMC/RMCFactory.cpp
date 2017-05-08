//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/RMC/RMCFactory.h"
#include "QMCDrivers/RMC/RMCSingleOMP.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{


  QMCDriver *RMCFactory::create (MCWalkerConfiguration & w,
				 TrialWaveFunction & psi, QMCHamiltonian & h,
				 ParticleSetPool & ptclpool,
				 HamiltonianPool & hpool,
				 WaveFunctionPool & ppool)
  {
    int np = omp_get_max_threads ();
    //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
    QMCDriver *qmc = 0;
#ifdef QMC_CUDA
    APP_ABORT("RMCFactory::create. RMC is not supported on GPU.\n");
#endif

    if (RMCMode == 0 || RMCMode == 1)	//(0,0,0) (0,0,1) pbyp and all electron
      {
	qmc = new RMCSingleOMP (w, psi, h, hpool, ppool);
      }
#if defined(QMC_BUILD_COMPLETE)
//else if(RMCMode == 2) //(0,1,0)
//{
//  qmc = new RMCMultiple(w,psi,h);
//}
//else if(RMCMode == 3) //(0,1,1)
//{
//  qmc = new RMCPbyPMultiple(w,psi,h);
//}
// else if(RMCMode ==2 || RMCMode ==3)
// {
//   qmc = new CSRMC(w,psi,h);
// }
// #if !defined(QMC_COMPLEX)
// else if(RMCMode == 6) //(1,1,0)
// {
//   qmc = new RMCMultipleWarp(w,psi,h, ptclpool);
// }
// else if(RMCMode == 7) //(1,1,1)
// {
//   qmc = new RMCPbyPMultiWarp(w,psi,h, ptclpool);
// }
// #endif
#endif
    qmc->setUpdateMode (RMCMode & 1);
    return qmc;
  }
}
