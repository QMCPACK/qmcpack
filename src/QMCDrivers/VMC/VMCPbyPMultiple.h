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
    
    


#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_MULTIPLE_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_MULTIPLE_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class MultipleEnergyEstimator;

/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 * @brief Implements the VMC algorithm
 */
class VMCPbyPMultiple: public QMCDriver
{
public:
  /// Constructor.
  VMCPbyPMultiple(MCWalkerConfiguration& w,
                  TrialWaveFunction& psi,
                  QMCHamiltonian& h);
  ~VMCPbyPMultiple();

  bool run();
  bool put(xmlNodePtr cur);

private:
  ///number of configurations
  int nPsi;
  ///number of blocks to equilibrate
  int equilBlocks;
  ///turn on/off drift for moves
  bool useDrift;
  ///option
  std::string useDriftOpt;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
  typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
  ParticleGradient_t dG;
  std::vector<ParticleGradient_t*> G;
  std::vector<ParticleLaplacian_t*> dL;
  std::vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
  MultipleEnergyEstimator *multiEstimator;
  ///resize the containers
  void resize(int ncopy, int nptcls);
  /// Copy Constructor (disabled)
  VMCPbyPMultiple(const VMCPbyPMultiple& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  VMCPbyPMultiple& operator=(const VMCPbyPMultiple&)
  {
    return *this;
  }
};

}

#endif
