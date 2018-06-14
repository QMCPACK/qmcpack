//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
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
  bool put(xmlNodePtr cur);
  void resetUpdateEngine();

private:
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
  DMCcuda(const DMCcuda& a) = delete;

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
  ///use T-moves
  int UseTMove;

  NewTimer ResizeTimer, DriftDiffuseTimer, BranchTimer, HTimer;
};
}

#endif
