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
    
    
    
    
    
    
    
    



#ifndef QMCPLUSPLUS_PAIRCORRELATION_ESTIMATOR_H
#define QMCPLUSPLUS_PAIRCORRELATION_ESTIMATOR_H
#include "Estimators/CompositeEstimators.h"
//#include "Estimators/VectorEstimatorImpl.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus
{

class GofREstimator: public CompositeEstimatorBase
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

  /** constructor
   * @param source particleset
   */
  GofREstimator(ParticleSet& source);

  /** constructor
   * @param source particleset
   * @param target particleset
   */
  GofREstimator(ParticleSet& source, ParticleSet& target);

  /** virtal destrctor */
  ~GofREstimator();

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

  CompositeEstimatorBase* clone();

private:
  GofREstimator(const GofREstimator& pc): sourcePtcl(pc.sourcePtcl) {}
};
}

#endif
