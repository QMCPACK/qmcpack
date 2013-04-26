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
#include "QMCDrivers/DMC/RNFactory.h"
#include "QMCDrivers/DMC/RNDMCOMP.h"
#include "Message/OpenMP.h"
//#define PETA_DMC_TEST
namespace qmcplusplus
{

QMCDriver* RNFactory::create(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide
                             , QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
{
  app_log() << "Creating RNDMCMP for the qmc driver" << endl;
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
