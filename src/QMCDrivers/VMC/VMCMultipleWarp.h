//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
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
