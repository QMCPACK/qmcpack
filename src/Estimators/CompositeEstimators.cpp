//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Estimators/CompositeEstimators.h"
#include "Estimators/PairCorrEstimator.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  CompositeEstimatorSet::CompositeEstimatorSet(ParticleSet& p): targetPtcl(p)
  {
    for(int i=0; i<p.DistTables.size(); i++)
    {
      if(p.DistTables[i]->origin().tag() == p.tag())
        Estimators.push_back(new PairCorrEstimator(targetPtcl));
      else
        Estimators.push_back(new PairCorrEstimator(p.DistTables[i]->origin(),targetPtcl));
    }
  }

  CompositeEstimatorSet::~CompositeEstimatorSet()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
  }

  void CompositeEstimatorSet::resetTargetParticleSet(ParticleSet& p)
  {
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->resetTargetParticleSet(p);
  }

  /** initiate sum over walkers
   */
  void CompositeEstimatorSet::startAccumulate()
  {
    curWeight=0.0;
    curStep++;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->startAccumulate();
  }

  /** accumulate data over walkers
   * @param P ParticleSet for the current walker
   * @param wgt weight
   */
  void CompositeEstimatorSet::accumulate(ParticleSet& P, RealType wgt) 
  {
    curWeight += wgt;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->accumulate(P,wgt);
  }

  /** After a sweep over walkers is completed, need to reweight the estimators
   * @param wgtnorm
   *
   * CompositeEstimatorSet calculates the sum over the walkers and pass it
   * to each estimator for proper reweighted sums.
   */
  void CompositeEstimatorSet::stopAccumulate(RealType wgtnorm)
  {
    wgtnorm=1.0/curWeight;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->stopAccumulate(wgtnorm);
  }


  /** reset the internal data of all the estimators for new averages
  */
  void CompositeEstimatorSet::startBlock(int steps) 
  {
    totSteps=steps;
    curStep=0;
    totWeight=0.0;
    curWeight=0.0;
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->startBlock(steps);
  }

  void CompositeEstimatorSet::stopBlock(RealType wgtnorm)
  {
    cout << " Block: steps =" << totSteps << " walker =" << curWeight << endl;
    //need to correct by the number of steps per block
    wgtnorm=1.0/static_cast<RealType>(totSteps);
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->stopBlock(wgtnorm);
  }

  //bool CompositeEstimatorSet::put(xmlNodePtr cur) 
  //{
  //  if(Estimators.empty())
  //  {
  //    //we can count the number of tables
  //    Estimators.push_back(new PairCorrEstimator(targetPtcl));
  //  }
  //  return true;
  //}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1800 $   $Date: 2007-02-22 17:33:31 -0600 (Thu, 22 Feb 2007) $
 * $Id: CompositeEstimatorSet.cpp 1800 2007-02-22 23:33:31Z jnkim $ 
 ***************************************************************************/
