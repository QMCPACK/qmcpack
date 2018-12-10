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
    
    



#include "Estimators/LocalEnergyEstimator.h"
#include <Utilities/IteratorUtility.h>

namespace qmcplusplus
{

LocalEnergyEstimator::LocalEnergyEstimator(QMCHamiltonian& h, bool use_hdf5)
  :UseHDF5(use_hdf5),refH(h)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  scalars.resize(SizeOfHamiltonians+LE_MAX);
  scalars_saved.resize(SizeOfHamiltonians+LE_MAX);
}

ScalarEstimatorBase* LocalEnergyEstimator::clone()
{
  return new LocalEnergyEstimator(*this);
}

void  LocalEnergyEstimator::registerObservables(std::vector<observable_helper*>& h5desc, hid_t gid)
{
  if(!UseHDF5)
    return;
  int loc=h5desc.size();
  //add LocalEnergy and LocalPotential
  h5desc.push_back(new observable_helper("LocalEnergy"));
  h5desc.push_back(new observable_helper("LocalEnergy_sq"));
  h5desc.push_back(new observable_helper("LocalPotential"));
  std::vector<int> onedim(1,1);
  h5desc[loc]->set_dimensions(onedim,FirstIndex);
  h5desc[loc++]->open(gid);
  h5desc[loc]->set_dimensions(onedim,FirstIndex+1);
  h5desc[loc++]->open(gid);
  h5desc[loc]->set_dimensions(onedim,FirstIndex+2);
  h5desc[loc++]->open(gid);
  //hamiltonian adds more
  refH.registerObservables(h5desc,gid);
  int correction=FirstIndex+3;
  for(int i=loc; i<h5desc.size(); ++i)
    h5desc[i]->lower_bound += correction;
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void LocalEnergyEstimator::add2Record(RecordListType& record)
{
  FirstIndex = record.size();
  int dumy=record.add("LocalEnergy");
  dumy=record.add("LocalEnergy_sq");
  dumy=record.add("LocalPotential");
  for(int i=0; i<SizeOfHamiltonians; ++i)
    record.add(refH.getObservableName(i));
  LastIndex=record.size();
  clear();
}

}
