//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
#ifndef QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#define QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 * @brief Implements the VMC algorithm
 */
class CSVMCUpdatePbyP: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdatePbyP(MCWalkerConfiguration& w,
                   vector<TrialWaveFunction*>& psi,
                   vector<QMCHamiltonian*>& h,
                  RandomGenerator_t& rg);

  ~CSVMCUpdatePbyP();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: CSVMCUpdatePbyP.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
