//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
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

QMCDriver* DMCFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi
                              , QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
{
#ifdef QMC_CUDA
  if (GPU)
    return new DMCcuda (w, psi, h,ppool);
#endif
  app_log() << "Creating DMCMP for the qmc driver" << std::endl;
  QMCDriver*  qmc = new DMCOMP(w,psi,h,hpool,ppool);
  qmc->setUpdateMode(PbyPUpdate);
  return qmc;
}
}
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
