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
    
    


#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_WARP_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_WARP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/SpaceWarp.h"
namespace qmcplusplus
{

class MultipleEnergyEstimator;
class ParticleSetPool;

/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 * @brief Implements the VMC algorithm
 */
class VMCPbyPMultiWarp: public QMCDriver
{
public:
  /// Constructor.
  VMCPbyPMultiWarp(MCWalkerConfiguration& w,
                   TrialWaveFunction& psi,
                   QMCHamiltonian& h,
                   ParticleSetPool& ptclPool);
  ~VMCPbyPMultiWarp();

  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  VMCPbyPMultiWarp(const VMCPbyPMultiWarp& a, ParticleSetPool& ptclPool): QMCDriver(a), PtclPool(ptclPool) { }
  /// Copy operator (disabled).
  VMCPbyPMultiWarp& operator=(const VMCPbyPMultiWarp&)
  {
    return *this;
  }

  ParticleSetPool& PtclPool;

  int nPsi,nptcl,JACOBIAN,equilBlocks;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
  typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
  ParticleGradient_t dG;
  std::vector<ParticleGradient_t*> G;
  std::vector<ParticleLaplacian_t*> dL;
  std::vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
  std::vector<ParticleSet*> WW;
  MultipleEnergyEstimator *multiEstimator;
  std::string refSetName;

  SpaceWarp PtclWarp;

  inline void resize(int ncopy, int nptcls)
  {
    int m=ncopy*(ncopy-1)/2;
    ratio.resize(ncopy);
    logpsi2.resize(ncopy);
    UmbrellaWeight.resize(ncopy);
    invsumratio.resize(ncopy);
    sumratio.resize(ncopy);
    ratioij.resize(m);
    for(int i=0; i<ncopy; i++)
    {
      G.push_back(new ParticleGradient_t(nptcls));
      dL.push_back(new ParticleLaplacian_t(nptcls));
    }
  }
};

}

#endif
