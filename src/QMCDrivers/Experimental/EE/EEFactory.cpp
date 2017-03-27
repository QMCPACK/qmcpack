//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/EE/EEFactory.h"
#include "QMCDrivers/EE/VMCRenyiOMP.h"
#include "QMCDrivers/EE/VMCMultiRenyiOMP.h"
#include "Message/OpenMP.h"


namespace qmcplusplus
{

QMCDriver* EEFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                             QMCHamiltonian& h, ParticleSetPool& ptclpool, HamiltonianPool& hpool, WaveFunctionPool& ppool)
{
  int np=omp_get_max_threads();
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  QMCDriver* qmc=0;
  if(VMCMode == 0 || VMCMode == 1) //(0,0,0) (0,0,1)
  {
    qmc = new VMCRenyiOMP(w,psi,h,hpool,ppool);
  }
  else
    if(VMCMode == 2 || VMCMode == 3)  //(0,1,0) (0,1,1)
    {
      qmc = new VMCMultiRenyiOMP(w,psi,h,hpool,ppool);
    }
  qmc->setUpdateMode(VMCMode&1);
  return qmc;
}
}
