//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Jeongnim Kim
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
#include <Estimators/CollectablesEstimator.h>
#include <QMCHamiltonians/observable_helper.h>

namespace qmcplusplus
{

CollectablesEstimator::CollectablesEstimator(QMCHamiltonian& h)
  : refH(h)
{
  scalars.resize(h.sizeOfCollectables());
  scalars_saved.resize(h.sizeOfCollectables());
}

void CollectablesEstimator::registerObservables(vector<observable_helper*>& h5desc
    , hid_t gid)
{
  int loc=h5desc.size();
  refH.registerCollectables(h5desc,gid);
  for(int i=loc; i<h5desc.size(); ++i)
    h5desc[i]->lower_bound += FirstIndex;
}

CollectablesEstimator* CollectablesEstimator::clone()
{
  return new CollectablesEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 *
 * Do not add Collectables to record
 */
void CollectablesEstimator::add2Record(RecordListType& record)
{
  FirstIndex=record.size();
  LastIndex=FirstIndex+scalars.size();
  //FirstIndex = record.size();
  //for(int i=0; i<refH.sizeOfCollectables(); ++i)
  //{
  //  ostringstream o;
  //  o<<"a"<<i;
  //  int dummy=record.add(o.str());
  //}
  //LastIndex = record.size();
  clear();
}

//void CollectablesEstimator::accumulate(const MCWalkerConfiguration& W
//    , WalkerIterator first, WalkerIterator last , RealType wgt)
//{
//  for(int i=0; i<refH.sizeOfCollectables(); ++i)
//    scalars[i](W.Collectables[i],wgt);
//}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3020 $   $Date: 2008-08-18 16:49:48 -0500 (Mon, 18 Aug 2008) $
 * $Id: LocalEnergyEstimatorHDF.cpp 3020 2008-08-18 21:49:48Z jnkim $
 ***************************************************************************/
