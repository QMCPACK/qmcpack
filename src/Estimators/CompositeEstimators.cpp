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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Estimators/CompositeEstimators.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  CompositeEstimatorSet::CompositeEstimatorSet():GroupID(-1),
  totSteps(0), curSteps(0)
  {
  }

  CompositeEstimatorSet::~CompositeEstimatorSet()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
  }

  ///not checking the map again, assuming that missing function
  void CompositeEstimatorSet::add(EstimatorType* est,const string& aname)
  {
    //map<strin,int>::iterator it(EstimatorMap.find(aname));
    //if(it == EstimatorMap.end())
    //{
      EstimatorMap[aname]=Estimators.size();
      Estimators.push_back(est);
    //}
  }

  void CompositeEstimatorSet::open(hid_t hroot)
  {
    GroupID=hroot;
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->open(hroot);
  }

  void CompositeEstimatorSet::close()
  {
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->close();
  }

  void CompositeEstimatorSet::resetTargetParticleSet(ParticleSet& p)
  {
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->resetTargetParticleSet(p);
  }

  /** accumulate data over walkers
   * @param W MCWalkerConfiguration
   * @param wgt weight
   */
  void CompositeEstimatorSet::accumulate(MCWalkerConfiguration& W) 
  {
    //cleanup
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->startAccumulate();

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
    if(PbyP)
    {
      while(wit != wit_end)
      {
        Walker_t& thisWalker(**wit);
        Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
        w_buffer.rewind();
        W.copyFromBuffer(w_buffer);
        for(int i=0; i< Estimators.size(); i++) Estimators[i]->accumulate(W);
        ++wit;
      }
    }
    else
    {
      while(wit != wit_end)
      {
        Walker_t& thisWalker(**wit);
        W.R=thisWalker.R;
        W.update();
        for(int i=0; i< Estimators.size(); i++) Estimators[i]->accumulate(W);
        ++wit;
      }
    }

    ///weight it by the total number of walkers per group
    RealType wgtnorm=1.0/static_cast<RealType>(W.getGlobalNumWalkers());
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->stopAccumulate(wgtnorm);

    curSteps++;
  }

  /** reset the internal data of all the estimators for new averages
   */
  void CompositeEstimatorSet::startBlock(int steps) 
  {
    totSteps=steps;
    curSteps=0;
    totWeight=0.0;
    curWeight=0.0;
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->startBlock(steps);
  }

  void CompositeEstimatorSet::stopBlock()
  {
    totSteps += curSteps;
    RealType wgtnorm=1.0/static_cast<RealType>(curSteps);
    RealType errnorm=1.0/static_cast<RealType>(totSteps-1);
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->stopBlock(wgtnorm,errnorm);
    curSteps=0;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1800 $   $Date: 2007-02-22 17:33:31 -0600 (Thu, 22 Feb 2007) $
 * $Id: CompositeEstimatorSet.cpp 1800 2007-02-22 23:33:31Z jnkim $ 
 ***************************************************************************/
