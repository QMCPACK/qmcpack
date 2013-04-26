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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_WFMCSINGLE_OMP_H
#define QMCPLUSPLUS_WFMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class WFMCSingleOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  WFMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                HamiltonianPool& hpool,WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
private:
  ///number of warmup steps
  int myWarmupSteps;
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///Index for Energy Index
  int Eindex;
  ///option to enable/disable drift equation for VMC
  string UseDrift,reweight;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  WFMCSingleOMP(const WFMCSingleOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  WFMCSingleOMP& operator=(const WFMCSingleOMP&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: VMCSingleOMP.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCSingleOMP.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/
