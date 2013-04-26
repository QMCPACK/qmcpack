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
  string useDriftOpt;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
  typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
  ParticleGradient_t dG;
  vector<ParticleGradient_t*> G;
  vector<ParticleLaplacian_t*> dL;
  vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCPbyPMultiple.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
