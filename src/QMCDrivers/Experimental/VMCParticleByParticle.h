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

  std::string UseDrift;

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
