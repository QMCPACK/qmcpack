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
#ifndef QMCPLUSPLUS_VMCSINGLE_OMP_H
#define QMCPLUSPLUS_VMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCSingleOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
               HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
  //inline vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  int prevSteps;
  int prevStepsBetweenSamples;

  ///Ways to set rn constant
  RealType logoffset,logepsilon;
  ///option to enable/disable drift equation or RN for VMC
  string UseDrift;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMCSingleOMP(const VMCSingleOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  VMCSingleOMP& operator=(const VMCSingleOMP&)
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
