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
//#include "Estimators/PairCorrEstimator.h"
#include "Estimators/GofREstimator.h"
#include "Estimators/SkEstimator.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  CompositeEstimatorSet::CompositeEstimatorSet(ParticleSet& p): targetPtcl(p), GroupID(-1)
  {
    //Disable gofr for th moment
    //for(int i=0; i<p.DistTables.size(); i++)
    //{
    //  if(p.DistTables[i]->origin().tag() == p.tag())
    //    Estimators.push_back(new GofREstimator(targetPtcl));
    //  else
    //    Estimators.push_back(new GofREstimator(p.DistTables[i]->origin(),targetPtcl));
    //}
    if(p.Lattice.SuperCellEnum) 
      Estimators.push_back(new SkEstimator(p));
  }

  CompositeEstimatorSet::~CompositeEstimatorSet()
  { 
    delete_iter(Estimators.begin(), Estimators.end());
  }

  void CompositeEstimatorSet::open(hid_t hroot)
  {
    if(hroot<0)
    {
      if(GroupID>=0) H5Fclose(GroupID);
      GroupID = hroot 
        = H5Fcreate("stuff.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    }
    for(int i=0; i< Estimators.size(); i++) Estimators[i]->open(hroot);
  }

  void CompositeEstimatorSet::close()
  {
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->close();
    if(GroupID>=0)
    {//responsible to close it
      H5Fclose(GroupID);
      GroupID=-1;
    }
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

  void CompositeEstimatorSet::stopBlock()
  {
    //need to correct by the number of steps per block
    RealType wgtnorm=1.0/static_cast<RealType>(totSteps);
    RealType errnorm=1.0/static_cast<RealType>(totSteps-1);
    for(int i=0; i< Estimators.size(); i++) 
      Estimators[i]->stopBlock(wgtnorm,errnorm);
  }

  ///** initiate sum over walkers
  // */
  //void CompositeEstimatorSet::startAccumulate()
  //{
  //  curWeight=0.0;
  //  curStep++;
  //  for(int i=0; i< Estimators.size(); i++) 
  //    Estimators[i]->startAccumulate();
  //}
  //
  ///** accumulate data over walkers
  // * @param P ParticleSet for the current walker
  // * @param wgt weight
  // */
  //void CompositeEstimatorSet::accumulate(ParticleSet& P, RealType wgt) 
  //{
  //  curWeight += wgt;
  //  for(int i=0; i< Estimators.size(); i++) 
  //    Estimators[i]->accumulate(P,wgt);
  //}
  //
  ///** After a sweep over walkers is completed, need to reweight the estimators
  // * @param wgtnorm
  // *
  // * CompositeEstimatorSet calculates the sum over the walkers and pass it
  // * to each estimator for proper reweighted sums.
  // */
  //void CompositeEstimatorSet::stopAccumulate(RealType wgtnorm)
  //{
  //  wgtnorm=(wgtnorm>0)? wgtnorm :1.0/curWeight;
  //  for(int i=0; i< Estimators.size(); i++) 
  //    Estimators[i]->stopAccumulate(wgtnorm);
  //}
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
