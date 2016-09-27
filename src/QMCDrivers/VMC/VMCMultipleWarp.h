//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  std::string refSetName;

  SpaceWarp PtclWarp;

  ///temporary storage
  int nPsi,equilBlocks;
  int nptcl;
  int JACOBIAN;
  std::vector<RealType> logpsi;
  std::vector<RealType> sumratio;
  std::vector<RealType> invsumratio;
  std::vector<RealType> Jacobian;
  std::vector<RealType> Norm;
  std::vector<ParticleSet*> WW;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultipleWarp.h 1593 2007-01-04 23:23:27Z jnkim $
 ***************************************************************************/
