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
    
    


#ifndef QMCPLUSPLUS_DMC_CUDA_H
#define QMCPLUSPLUS_DMC_CUDA_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the DMC algorithm using particle-by-particle move.
 */
class DMCcuda: public QMCDriver
{
public:
  /// Constructor.
  DMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool);
  bool run();
  bool runWithNonlocal();
  bool put(xmlNodePtr cur);
  void resetUpdateEngine();

private:
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  std::string ScaleWeight;
  /// tau/mass
  RealType m_tauovermass;
  ///steps before branching
  int BranchInterval;
  ///number of warmup steps
  int myWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///update engine
  QMCUpdateBase* Mover;
  /// Copy Constructor (disabled)
  DMCcuda(const DMCcuda& a): QMCDriver(a),
    ResizeTimer("DMCcuda::resize"),
    DriftDiffuseTimer("DMCcuda::Drift/Diffuse"),
    BranchTimer("DMCcuda::Branch"),
    HTimer("DMCcuda::Hamiltonian")

  { }

  bool checkBounds (const PosType &newpos);
  void checkBounds (std::vector<PosType> &newpos, std::vector<bool> &valid);

  /// Copy operator (disabled).
  DMCcuda& operator=(const DMCcuda&)
  {
    return *this;
  }
  ///hide initialization from the main function
  void resetRun();
  NonLocalTOperator NLop;

  NewTimer ResizeTimer, DriftDiffuseTimer, BranchTimer, HTimer;
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCcuda.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: DMCcuda.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
