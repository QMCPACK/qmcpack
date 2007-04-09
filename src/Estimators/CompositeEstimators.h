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
#ifndef QMCPLUSPLUS_COMPOSITE_ESTIMATORBASE_H
#define QMCPLUSPLUS_COMPOSITE_ESTIMATORBASE_H
#include "Configuration.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {

  /** Abstract class for an estimator of an operator.
   */
  struct CompositeEstimatorBase: public QMCTraits {

    bool CollectSum;

    /** default constructor */
    CompositeEstimatorBase(): CollectSum(false){}

    /** virtal destrctor */
    virtual ~CompositeEstimatorBase() {}

    /** reassign the target particleset */
    virtual void resetTargetParticleSet(ParticleSet& p) = 0;

    /** start accumulate */
    virtual void startAccumulate()=0;

    /** accumulate the observables */
    virtual void accumulate(ParticleSet& p, RealType wgt)=0;

    /** stop accumulate for an ensemble and reweight the data */
    virtual void stopAccumulate(RealType wgtnorm)=0;

    /** start a block */
    virtual void startBlock(int steps)=0;

    /** stop a block */
    virtual void stopBlock(RealType wgtnorm)=0;
  };

  /**Class to manage a set of ScalarEstimators */
  struct CompositeEstimatorSet: public CompositeEstimatorBase
  {

    typedef CompositeEstimatorBase EstimatorType;

    CompositeEstimatorSet(ParticleSet& p);
    ~CompositeEstimatorSet();

    ///implement virutal functions
    void resetTargetParticleSet(ParticleSet& p);

    void startAccumulate();
    void accumulate(ParticleSet& P, RealType wgt);
    void stopAccumulate(RealType wgtnorm);

    void report(int iter);
    void reset();

    void startBlock(int steps);
    void stopBlock(RealType wgtnorm);
    ///number of steps per block
    int totSteps;
    ///current step
    int curStep;
    ///total weight during a block
    RealType totWeight;
    ///weight during a step
    RealType curWeight;
    ///target particles
    ParticleSet& targetPtcl;
    ///estimators
    vector<EstimatorType*> Estimators;
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/
