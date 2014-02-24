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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file CSUpdateBase.h
 * @brief Definition of CSVUpdateBase
 */
#ifndef QMCPLUSPLUS_CORRELATEDSAMPLINGBASE_H
#define QMCPLUSPLUS_CORRELATEDSAMPLINGBASE_H

#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

class CSEnergyEstimator;
class CSUpdateBase: public QMCUpdateBase
{
public:

  CSUpdateBase(MCWalkerConfiguration& w, vector<TrialWaveFunction*>& psi, vector<QMCHamiltonian*>& h,
               RandomGenerator_t& rg);

  virtual ~CSUpdateBase();

  int nPsi;
  bool useDrift;
  vector<RealType> logpsi;
  vector<RealType> sumratio;
  vector<RealType> invsumratio;
  vector<RealType> avgNorm;
  vector<RealType> avgWeight;
  vector<RealType> logNorm;
  vector<RealType> cumNorm;
  vector<RealType> instRij;
  Matrix<RealType> ratioIJ;

  string useDriftOption;

  ///multiple estimator
  CSEnergyEstimator* multiEstimator;
  ///a list of TrialWaveFunctions for multiple method
  vector<TrialWaveFunction*> Psi1;
  ///a list of QMCHamiltonians for multiple method
  vector<QMCHamiltonian*> H1;

  vector<ParticleSet::ParticleGradient_t*> G1;
  vector<ParticleSet::ParticleLaplacian_t*> L1;

  void resizeWorkSpace(int nw, int nptcls);
  void updateNorms();
  void updateAvgWeights();
  void initCSWalkers(WalkerIter_t it, WalkerIter_t it_end, bool resetNorms);
  void initCSWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end, bool resetNorms);
  void updateCSWalkers(WalkerIter_t it, WalkerIter_t it_end);
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
