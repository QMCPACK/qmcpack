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
/**@file CSVMC.h
 * @brief Definition of CSVMC
 */
#ifndef QMCPLUSPLUS_CS_VMCMULTIPLE_H
#define QMCPLUSPLUS_CS_VMCMULTIPLE_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class CSEnergyEstimator;
class CSUpdateBase;

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class CSVMC: public QMCDriver
{
public:
  /// Constructor.
  CSVMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur);

private:
  ///blocks over which normalization factors are accumulated
  int equilBlocks;
  /// Copy Constructor (disabled)
  CSVMC(const CSVMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  CSVMC& operator=(const CSVMC&)
  {
    return *this;
  }

  CSEnergyEstimator *multiEstimator;
  CSUpdateBase* Mover;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: CSVMC.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
