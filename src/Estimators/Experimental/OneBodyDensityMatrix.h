//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
#ifndef QMCPLUSPLUS_nofr_ESTIMATOR_H
#define QMCPLUSPLUS_nofr_ESTIMATOR_H
#include "Estimators/CompositeEstimators.h"
//#include "Estimators/VectorEstimatorImpl.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus
{
class TrialWaveFunction;
class nofrEstimator: public CompositeEstimatorBase
{
  //typedef VectorEstimatorImpl<RealType> VectorEstimatorType;
  ///true if source == target
  bool Symmetric;
  /** number of centers */
  int Centers;
  /** number of distinct pair types */
  int NumPairTypes;
  /** number bins for gofr */
  int NumBins;
  /** maximum distance */
  RealType Dmax;
  /** bin size */
  RealType Delta;
  /** one of bin size */
  RealType DeltaInv;
  ///save the source particleset
  ParticleSet* sourcePtcl;
  ///save the target particleset
  ParticleSet* targetPtcl;
  /** distance table */
  const DistanceTableData*  myTable;
  /** local copy of pair index */
  std::vector<int> PairID;
  /** normalization factor for each bin*/
  Vector<RealType> normFactor;
  /** instantaneous gofr */
  Matrix<RealType> gofrInst;

public:
  void PutInBox (PosType &v);
  /** constructor
   * @param source particleset
   */
  nofrEstimator(TrialWaveFunction &psi,MCWalkerConfiguration& w);
  TrialWaveFunction &Psi;
  MCWalkerConfiguration& W;
  /** constructor
   * @param source particleset
   * @param target particleset
   */
  //    nofrEstimator(ParticleSet& source, ParticleSet& target);

  /** virtal destrctor */
  ~nofrEstimator();

  /** create gofr group */
  hid_t createGroup(hid_t gid);
  //void writeHeaders(hid_t gid);
  void resetTargetParticleSet(ParticleSet& p);
  /** prepare data collect */
  void startAccumulate();
  /** accumulate the observables */
  void accumulate(ParticleSet& p);
  /** reweight of the current cummulative  values */
  void stopAccumulate(RealType wgtnorm);
  void setBound(RealType dr);
  int NumSamples;
  CompositeEstimatorBase* clone();

private:

};
}

#endif
