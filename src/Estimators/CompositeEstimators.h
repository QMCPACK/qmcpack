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
#include "OhmmsData/HDFAttribIO.h"
#include "Particle/MCWalkerConfiguration.h"

/* accumulate weighted squares
 * @param first starting iterator of input data
 * @param last ending iterator for input data
 * @param target starting iterator of the accumulation
 * @param w weight
 *
 * target[i] += w*soure[i]*source[i];
 */
template<typename IT1, typename IT2, typename T>
inline void accumulate2(IT1 first, IT1 last, IT2 target, T w)
{
  while(first != last)
  {
    *target += w*(*first)*(*first); 
    ++target; ++first;
  }
}


namespace qmcplusplus {

  /** Abstract class for an estimator of an operator.
   */
  struct CompositeEstimatorBase: public QMCTraits {

    ///hdf5 handle of the object
    hid_t GroupID;
    ///name of the object
    string Title;
    /** default constructor */
    CompositeEstimatorBase(): GroupID(-1){}

    /** virtal destrctor */
    virtual ~CompositeEstimatorBase() {}

    /** reassign the target particleset */
    virtual void resetTargetParticleSet(ParticleSet& p) = 0;

    /** initialize the estimator */
    virtual void open(hid_t hroot)=0;

    /** finalize the estimator */
    virtual void close()=0;

    /** start accumulate */
    virtual void startAccumulate()=0;

    /** accumulate the observables */
    virtual void accumulate(ParticleSet& p, RealType wgt)=0;

    /** stop accumulate for an ensemble and reweight the data */
    virtual void stopAccumulate(RealType wgtnorm)=0;

    /** start a block */
    virtual void startBlock(int steps)=0;

    /** stop a block 
     * @param wgtnorm for average
     * @param errnorm for error normalization 1.0/(samples-1)
     */
    virtual void stopBlock(RealType wgtnorm, RealType errnorm)=0;
  };

  /**Class to manage a set of ScalarEstimators */
  struct CompositeEstimatorSet: public CompositeEstimatorBase
  {

    typedef CompositeEstimatorBase EstimatorType;

    CompositeEstimatorSet(ParticleSet& p);
    ~CompositeEstimatorSet();

    void resetTargetParticleSet(ParticleSet& p);
    void open(hid_t hroot);
    void close();
    void startAccumulate();
    void accumulate(ParticleSet& P, RealType wgt);
    void stopAccumulate(RealType wgtnorm);

    void report(int iter);
    void reset();

    void startBlock(int steps);
    void stopBlock(RealType wgtnorm, RealType errnorm);
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
