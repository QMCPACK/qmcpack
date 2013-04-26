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
#include <Estimators/LocalEnergyEstimatorHDF.h>
#include <Utilities/IteratorUtility.h>

namespace qmcplusplus
{

LocalEnergyEstimatorHDF::LocalEnergyEstimatorHDF(QMCHamiltonian& h)
  :refH(h)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  scalars.resize(SizeOfHamiltonians+LE_MAX);
  scalars_saved.resize(SizeOfHamiltonians+LE_MAX);
}

void LocalEnergyEstimatorHDF::registerObservables(vector<observable_helper*>& h5desc
    ,hid_t gid)
{
  int loc=h5desc.size();
  //add LocalEnergy and LocalPotential
  h5desc.push_back(new observable_helper("LocalEnergy"));
  h5desc.push_back(new observable_helper("LocalPotential"));
  std::vector<int> onedim(1,1);
  h5desc[loc]->set_dimensions(onedim,FirstIndex);
  h5desc[loc++]->open(gid);
  h5desc[loc]->set_dimensions(onedim,FirstIndex+1);
  h5desc[loc++]->open(gid);
  //hamiltonian adds more
  refH.registerObservables(h5desc,gid);
  int correction=FirstIndex+2;
  for(int i=loc; i<h5desc.size(); ++i)
    h5desc[i]->lower_bound += correction;
}

ScalarEstimatorBase* LocalEnergyEstimatorHDF::clone()
{
  return new LocalEnergyEstimatorHDF(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void LocalEnergyEstimatorHDF::add2Record(RecordListType& record)
{
  FirstIndex = record.size();
  int dumy=record.add("LocalEnergy");
  dumy=record.add("LocalPotential");
  for(int i=0; i<SizeOfHamiltonians; ++i)
    record.add(refH.getObservableName(i));
  LastIndex=record.size();
  clear();
}

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3020 $   $Date: 2008-08-18 16:49:48 -0500 (Mon, 18 Aug 2008) $
 * $Id: LocalEnergyEstimatorHDF.cpp 3020 2008-08-18 21:49:48Z jnkim $
 ***************************************************************************/
