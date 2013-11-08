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
  ///input string to determine to use nonlocal move
  string NonLocalMove;
  string ScaleWeight;
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
  void checkBounds (vector<PosType> &newpos, vector<bool> &valid);

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
