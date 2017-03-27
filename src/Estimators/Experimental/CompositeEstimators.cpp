//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    



#include "Estimators/CompositeEstimators.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

CompositeEstimatorBase::~CompositeEstimatorBase()
{
  close();
  delete_iter(dList.begin(), dList.end());
}

/** start a new block
 * @param nsteps number of steps per block
 *
 * Clear all the internal data to start blocking.
 */
void CompositeEstimatorBase::startBlock(int steps)
{
  for(int p=0; p<dList.size(); p++)
    dList[p]->init();
}

/** stop a block
 * @param wgtnorm normalization factor
 */
void CompositeEstimatorBase::stopBlock(RealType wgtnorm)
{
  for(int p=0; p<dList.size(); p++)
    dList[p]->takeBlockAverage(wgtnorm);
}

/** record a block data
 */
void CompositeEstimatorBase::recordBlock()
{
  for(int p=0; p<oList.size(); p++)
    oList[p]->write(GroupID,0);
}

void CompositeEstimatorBase::open(hid_t hroot)
{
  if(GroupID<0)
  {
    int n=nList.size();
    oList.resize(n,0);
    GroupID=createGroup(hroot);
    for(int p=0; p<n; p++)
    {
      hid_t gid = H5Gcreate(GroupID,nList[p].c_str(),0);
      oList[p]= new HDFAttribIO<VectorEstimatorType>(*dList[p]);
      oList[p]->reserve(gid);
    }
  }
}

void CompositeEstimatorBase::close()
{
  if(GroupID>-1)
  {
    //for(int p=0; p<hList.size(); p++) H5Gclose(hList[p]);
    delete_iter(oList.begin(),oList.end());
    oList.clear();
    //hList.clear();
    H5Gclose(GroupID);//close all the group
    GroupID=-1;
  }
}

//  ///add measurements
//  void CompositeEstimatorBase::addMeasurements(int n)
//  {
//    delete_iter(dList.begin(), dList.end());
//    dList.resize(n,0);
//  }

///add block
void CompositeEstimatorBase::collectBlock(CompositeEstimatorBase* eth)
{
  for(int i=0; i<dList.size(); i++)
  {
    dList[i]->d_sum += eth->dList[i]->d_sum;
    //dList[i]->d_sum2 += eth->dList[i]->d_sum2;
  }
}

///////////////////////////////////////////////////////////////
// definitions of CompositeEstimatorSet
///////////////////////////////////////////////////////////////
CompositeEstimatorSet::CompositeEstimatorSet():
  GroupID(-1), NumSteps(1),OneOverNumSteps(1.0)
{
}

CompositeEstimatorSet::CompositeEstimatorSet(const CompositeEstimatorSet& ce):
  GroupID(-1), NumSteps(1),OneOverNumSteps(1.0)
{
  std::map<std::string,int>::const_iterator it(ce.EstimatorMap.begin());
  std::map<std::string,int>::const_iterator it_end(ce.EstimatorMap.end());
  while(it != it_end)
  {
    add(ce.Estimators[(*it).second]->clone(),(*it).first);
    ++it;
  }
}

CompositeEstimatorSet::~CompositeEstimatorSet()
{
  delete_iter(Estimators.begin(), Estimators.end());
  //close the file
  if(GroupID>=-1)
    H5Fclose(GroupID);
}

///not checking the map again, assuming that missing function
void CompositeEstimatorSet::add(EstimatorType* est,const std::string& aname)
{
  //map<strin,int>::iterator it(EstimatorMap.find(aname));
  //if(it == EstimatorMap.end())
  //{
  est->Title=aname; //assign the title
  EstimatorMap[aname]=Estimators.size();
  Estimators.push_back(est);
  //}
}

void CompositeEstimatorSet::open(hid_t hroot)
{
  GroupID=hroot;
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->open(hroot);
}

void CompositeEstimatorSet::close()
{
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->close();
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
void CompositeEstimatorSet::accumulate(MCWalkerConfiguration& W, RealType wgtnorm)
{
  std::accumulate(W,W.begin(),W.end(),wgtnorm);
}

void CompositeEstimatorSet::accumulate(MCWalkerConfiguration& W,
                                       MCWalkerConfiguration::iterator wit,
                                       MCWalkerConfiguration::iterator wit_end,
                                       RealType wgtnorm)
{
  //initialize temporary data
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->startAccumulate();
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  if(W.updatePbyP())
  {
    while(wit != wit_end)
    {
      Walker_t& thisWalker(**wit);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      for(int i=0; i< Estimators.size(); i++)
        Estimators[i]->accumulate(W);
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
      for(int i=0; i< Estimators.size(); i++)
        Estimators[i]->accumulate(W);
      ++wit;
    }
  }
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->stopAccumulate(wgtnorm);
  //curSteps++;
}

/** reset the internal data of all the estimators for new averages
 */
void CompositeEstimatorSet::startBlock(int steps)
{
  NumSteps=steps;
  OneOverNumSteps=1.0/static_cast<RealType>(steps);
  //  curSteps=0;
  //  totSteps=0;
  //  totWeight=0.0;
  //  curWeight=0.0;
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->startBlock(steps);
}

void CompositeEstimatorSet::stopBlock()
{
  //  totSteps += curSteps;
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->stopBlock(OneOverNumSteps);
  //  curSteps=0;
}

void CompositeEstimatorSet::collectBlock(CompositeEstimatorSet* eth)
{
  //  curSteps+=eth->curSteps;
  for(int i=0; i< Estimators.size(); i++)
  {
    Estimators[i]->collectBlock(eth->Estimators[i]);
  }
}

void CompositeEstimatorSet::recordBlock()
{
  for(int i=0; i< Estimators.size(); i++)
    Estimators[i]->recordBlock();
  H5Fflush(GroupID,H5F_SCOPE_LOCAL);
}
}

