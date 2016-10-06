//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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

  CSUpdateBase(MCWalkerConfiguration& w, std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& h,
               RandomGenerator_t& rg);

  virtual ~CSUpdateBase();

  int nPsi;
  bool useDrift;
  std::vector<RealType> logpsi;
  std::vector<RealType> sumratio;
  std::vector<RealType> invsumratio;
  std::vector<RealType> avgNorm;
  std::vector<RealType> avgWeight;
  std::vector<RealType> logNorm;
  std::vector<RealType> cumNorm;
  std::vector<RealType> instRij;
  Matrix<RealType> ratioIJ;
  std::string useDriftOption;

  ///multiple estimator
  CSEnergyEstimator* multiEstimator;
  ///a list of TrialWaveFunctions for multiple method
  std::vector<TrialWaveFunction*> Psi1;
  ///a list of QMCHamiltonians for multiple method
  std::vector<QMCHamiltonian*> H1;

  std::vector<ParticleSet::ParticleGradient_t*> G1;
  std::vector<ParticleSet::ParticleLaplacian_t*> L1;

  //a list of single particle gradients for multiple wavefunctions;
  std::vector<GradType> g1_old;
  std::vector<GradType> g1_new;

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
