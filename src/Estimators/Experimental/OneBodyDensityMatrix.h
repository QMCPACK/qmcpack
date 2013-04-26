//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
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
  vector<int> PairID;
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $
 ***************************************************************************/
