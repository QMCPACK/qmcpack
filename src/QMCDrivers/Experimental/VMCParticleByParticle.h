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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCParticleByParticle: public QMCDriver
{
public:
  /// Constructor.
  VMCParticleByParticle(MCWalkerConfiguration& w,
                        TrialWaveFunction& psi,
                        QMCHamiltonian& h);
  bool run();
  bool put(xmlNodePtr cur);

private:
  IndexType nAcceptTot;
  IndexType nRejectTot;
  IndexType nAllRejected;
  IndexType nSubSteps;

  string UseDrift;

  ParticleSet::ParticleGradient_t G, dG;
  ParticleSet::ParticleLaplacian_t L, dL;

  void runBlockWithDrift();
  void runBlockWithoutDrift();
  //void advanceWalkerByWalker();
  /// Copy Constructor (disabled)
  VMCParticleByParticle(const VMCParticleByParticle& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCParticleByParticle& operator=(const VMCParticleByParticle&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCParticleByParticle.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
