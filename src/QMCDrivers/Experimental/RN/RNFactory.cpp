//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/RNFactory.h"
#include "QMCDrivers/DMC/RNDMCOMP.h"
#include "Message/OpenMP.h"
//#define PETA_DMC_TEST
namespace qmcplusplus
{

QMCDriver* RNFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide
                             , QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
{
  app_log() << "Creating RNDMCMP for the qmc driver" << std::endl;
  QMCDriver*  qmc = new RNDMCOMP(w,psi,guide,h,hpool,ppool);
  qmc->setUpdateMode(PbyPUpdate);
  qmc->put(myNode);
  return qmc;
}
}
/***************************************************************************
 * $RCSfile: DMCFactory.cpp,v $   $Author: jmcminis $
 * $Revision: 4215 $   $Date: 2009-09-21 22:38:50 -0500 (Mon, 21 Sep 2009) $
 * $Id: DMCFactory.cpp 4215 2009-09-22 03:38:50Z jmcminis $
 ***************************************************************************/
