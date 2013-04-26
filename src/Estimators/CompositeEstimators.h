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
#include "Estimators/VectorEstimatorImpl.h"

namespace qmcplusplus
{

/** Abstract class for an estimator of an operator.
 */
struct CompositeEstimatorBase: public QMCTraits
{

  typedef VectorEstimatorImpl<RealType> VectorEstimatorType;

  ///hdf5 handle of the object
  hid_t GroupID;
  ///name of the object
  string Title;
  ///name of associated data
  vector<string>                            nList;
  ///VectorEstimatorType
  vector<VectorEstimatorType*>              dList;
  ///h5 engine
  vector<HDFAttribIO<VectorEstimatorType>*> oList;

  /** default constructor */
  CompositeEstimatorBase(): GroupID(-1) {}

  /** virtal destrctor */
  virtual ~CompositeEstimatorBase();

  /** virtual function to enable deep copy for threaded applications*/
  virtual CompositeEstimatorBase* clone() = 0;
  /** reassign the target particleset */
  virtual void resetTargetParticleSet(ParticleSet& p) = 0;

  /** start accumulate */
  virtual void startAccumulate()=0;

  /** accumulate the observables */
  virtual void accumulate(ParticleSet& p)=0;

  /** stop accumulate for an ensemble and reweight the data */
  virtual void stopAccumulate(RealType wgtnorm)=0;

  /** create a group for a st of estimators */
  virtual hid_t createGroup(hid_t gid)=0;
  //virtual void writeHeaders(hid_t gid)=0;
  /** initialize the estimator IO */
  void open(hid_t hroot);

  /** finalize the estimator */
  void close();

  /** start a block */
  void startBlock(int steps);

  /** stop a block
   * @param wgtnorm for average
   * @param errnorm for error normalization 1.0/(samples-1)
   */
  void stopBlock(RealType wgtnorm);

  ///record block
  void recordBlock();

  /** collect data from eth
   * @param eth CompositeEstimatorBase on the data collected by a thread
   */
  void collectBlock(CompositeEstimatorBase* eth);

  /** return size of the data handled by this estimator
   */
  inline int size() const
  {
    int n=0;
    for(int i=0; i<dList.size(); i++)
      n += dList[i]->size();
    return n;
  }

  template<typename ForwardIterator>
  ForwardIterator putMessage(ForwardIterator cur) const
  {
    for(int i=0; i<dList.size(); i++)
      cur=dList[i]->putMessage(cur);
    return cur;
  }

  template<typename ForwardIterator>
  ForwardIterator getMessage(ForwardIterator cur)
  {
    for(int i=0; i<dList.size(); i++)
      cur=dList[i]->getMessage(cur);
    return cur;
  }

  void print(ostream& os)
  {
    for(int i=0; i<dList.size(); i++)
    {
      os << setw(3) << i;
      dList[i]->print(os);
    }
  }
};

/** Class to manage a set of CompositeEstimatorBase
 */
struct CompositeEstimatorSet: public QMCTraits
{

  ///typedef estimator type is CompositeEstimatorBase
  typedef CompositeEstimatorBase EstimatorType;
  ///number of steps for block average
  int NumSteps;
  /////total number of steps
  //int totSteps;
  /////current step in a block
  //int curSteps;
  ///1.0/NumSteps
  RealType OneOverNumSteps;
  /////total weight during a block
  //RealType totWeight;
  /////weight during a step
  //RealType curWeight;
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
  ///copy constructor
  CompositeEstimatorSet(const CompositeEstimatorSet& ce);
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
   * @param wgtnorm normalization factor
   */
  void startBlock(int steps);

  /** return size of the data handled by this estimator
   */
  inline int size() const
  {
    int n=0;
    for(int i=0; i<Estimators.size(); i++)
      n+= Estimators[i]->size();
    return n;
  }
  /** accumulate measurements */
  void accumulate(MCWalkerConfiguration& W, RealType wgtnorm);
  /** accumulate measurements
   * @param W particleset to evaluate quantities
   * @param it first walker
   * @param it_end last walker
   */
  void accumulate(MCWalkerConfiguration& W,
                  MCWalkerConfiguration::iterator it, MCWalkerConfiguration::iterator it_end, RealType wgtnorm);
  /** stop recording the block
   * @param wgtnorm normalization factor
   */
  void stopBlock();
  /** collect blocks from other estimators
   * @param eth estimator to be added
   * @param wgtnorm normalization factor
   *
   * For threaded applications.
   */
  void collectBlock(CompositeEstimatorSet* eth);
  void recordBlock();
  void reset();

  template<typename ForwardIterator>
  ForwardIterator putMessage(ForwardIterator cur) const
  {
    for(int i=0; i<Estimators.size(); i++)
      cur=Estimators[i]->putMessage(cur);
    return cur;
  }

  template<typename ForwardIterator>
  ForwardIterator getMessage(ForwardIterator cur)
  {
    for(int i=0; i<Estimators.size(); i++)
      cur=Estimators[i]->getMessage(cur);
    return cur;
  }

  void print(ostream& os)
  {
    for(int i=0; i<Estimators.size(); i++)
      Estimators[i]->print(os);
  }
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $
 ***************************************************************************/
