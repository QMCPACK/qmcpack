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
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author: jnkim $
 * $Revision: 1.3 $   $Date: 2006/04/05 00:49:59 $
 * $Id: DMCFactory.cpp,v 1.3 2006/04/05 00:49:59 jnkim Exp $
 ***************************************************************************/
