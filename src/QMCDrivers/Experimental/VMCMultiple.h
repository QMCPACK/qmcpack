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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file VMCMultiple.h
 * @brief Definition of VMCMultiple
 */
#ifndef QMCPLUSPLUS_VMCMULTIPLE_H_PSI_H
#define QMCPLUSPLUS_VMCMULTIPLE_H_PSI_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class MultipleEnergyEstimator;

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class VMCMultiple: public QMCDriver
{
public:
  /// Constructor.
  VMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  void advanceWalkerByWalker();
  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCMultiple(const VMCMultiple& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCMultiple& operator=(const VMCMultiple&)
  {
    return *this;
  }

  MultipleEnergyEstimator *multiEstimator;
  ///temporary storage
  int nPsi;
  ///number of blocks to compute the normalization factor
  int equilBlocks;
  vector<RealType> logpsi;
  vector<RealType> sumratio;
  vector<RealType> invsumratio;
  vector<RealType> Norm;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
