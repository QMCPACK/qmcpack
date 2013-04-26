//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
  vector<ParticleGradient_t*> G;
  vector<ParticleLaplacian_t*> dL;
  vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
  vector<ParticleSet*> WW;
  MultipleEnergyEstimator *multiEstimator;
  string refSetName;

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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCPbyPMultiWarp.h 1593 2007-01-04 23:23:27Z jnkim $
 ***************************************************************************/
