//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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

void CollectablesEstimator::registerObservables(std::vector<observable_helper*>& h5desc
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
  //  std::ostringstream o;
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
