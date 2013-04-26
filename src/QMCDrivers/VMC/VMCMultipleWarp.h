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
/**@file VMCMultipleWarp.h
 * @brief Definition of VMCMultipleWarp
 */
#ifndef QMCPLUSPLUS_QMC_VMCWARP_H_PSI_H
#define QMCPLUSPLUS_QMC_VMCWARP_H_PSI_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/SpaceWarp.h"
namespace qmcplusplus
{

class MultipleEnergyEstimator;
class ParticleSetPool;

/** @ingroup QMCDrivers WalkerByWalker MultiplePsi
 * @brief Implements the VMC algorithm using umbrella sampling.
 *
 * Energy difference method with multiple H/Psi.
 * Consult S. Chiesa's note.
 */
class VMCMultipleWarp: public QMCDriver
{
public:
  /// Constructor.
  VMCMultipleWarp(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, ParticleSetPool& ptclPool);

  void advanceWalkerByWalker();
  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCMultipleWarp(const VMCMultipleWarp& a, ParticleSetPool& ptclPool):
    QMCDriver(a), PtclPool(ptclPool) { }
  /// Copy operator (disabled).
  VMCMultipleWarp& operator=(const VMCMultipleWarp&)
  {
    return *this;
  }

  ParticleSetPool& PtclPool;
  MultipleEnergyEstimator *multiEstimator;
  string refSetName;

  SpaceWarp PtclWarp;

  ///temporary storage
  int nPsi,equilBlocks;
  int nptcl;
  int JACOBIAN;
  vector<RealType> logpsi;
  vector<RealType> sumratio;
  vector<RealType> invsumratio;
  vector<RealType> Jacobian;
  vector<RealType> Norm;
  vector<ParticleSet*> WW;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultipleWarp.h 1593 2007-01-04 23:23:27Z jnkim $
 ***************************************************************************/
