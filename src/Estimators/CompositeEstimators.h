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
    virtual void accumulate(ParticleSet& p)=0;

    /** stop accumulate for an ensemble and reweight the data */
    virtual void stopAccumulate(RealType wgtnorm)=0;

    /** start a block */
    virtual void startBlock(int steps)=0;

    /** stop a block 
     * @param wgtnorm for average
     * @param errnorm for error normalization 1.0/(samples-1)
     */
    virtual void stopBlock(RealType wgtnorm, RealType errnorm)=0;

    /* accumulate weighted squares
     * @param first starting iterator of input data
     * @param last ending iterator for input data
     * @param v starting iterator for the sum
     * @param v2 starting iterator for the squred sum
     * @param w weight
     *
     * v[i] += w*soure[i];
     * v2[i] += w*soure[i]*source[i];
     */
    template<typename IT1, typename IT2, typename T>
      inline void collect(IT1 first, IT1 last, IT2 v, IT2 v2, T w)
      {
        while(first != last)
        {
          *v2++ += w*(*first)*(*first); 
          *v++  += w*(*first++);
        }
      }

  };

  /** Class to manage a set of CompositeEstimatorBase
   */
  struct CompositeEstimatorSet: public QMCTraits
  {

    ///typedef estimator type is CompositeEstimatorBase
    typedef CompositeEstimatorBase EstimatorType;
    ///true if the move was particle by particle
    bool PbyP;
    ///number of steps per block
    int totSteps;
    ///current step
    int curStep;
    ///total weight during a block
    RealType totWeight;
    ///weight during a step
    RealType curWeight;
    ///hdf5 handle of the object
    hid_t GroupID;
    ///name of the object
    string Title;
    ///estimators
    vector<EstimatorType*> Estimators;
    ///name map
    map<string,int> EstimatorMap;

    ///constructor
    //CompositeEstimatorSet(ParticleSet& p);
    CompositeEstimatorSet();
    ///destructor
    ~CompositeEstimatorSet();

    /** return true if aname does not exisit
     */
    bool missing(const string& aname)
    {
      return EstimatorMap.find(aname) == EstimatorMap.end();
    }
    /** add estimator
     * @param est a new estimator
     * @param aname the name of the new estimator
     */
    void add(EstimatorType* est, const string& aname) ;
    ///reset the target particle set
    void resetTargetParticleSet(ParticleSet& p);

    ///open a h5group to record the estimators.
    void open(hid_t hroot);
    ///close GroupID;
    void close();

    /** start a block to record
     * @param steps number of steps for a block
     */
    void startBlock(int steps);
    /** accumulate the measurements */
    void accumulate(MCWalkerConfiguration& W);
    /** stop recording the block */
    void stopBlock();
    void reset();
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/
