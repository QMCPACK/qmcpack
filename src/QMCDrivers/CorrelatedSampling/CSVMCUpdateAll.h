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
/**@file CSVMCUpdateAll.h
 * @brief Definition of CSVMCUpdateAll
 */
#ifndef QMCPLUSPLUS_CS_VMC_UPDATEALL_H
#define QMCPLUSPLUS_CS_VMC_UPDATEALL_H

#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class CSVMCUpdateAll: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdateAll(MCWalkerConfiguration& w,  vector<TrialWaveFunction*>& psi, vector<QMCHamiltonian*>& h,
                 RandomGenerator_t& rg);

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
};


class CSVMCUpdateAllWithDrift: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdateAllWithDrift(MCWalkerConfiguration& w,  vector<TrialWaveFunction*>& psi, vector<QMCHamiltonian*>& h,
                 RandomGenerator_t& rg);

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
