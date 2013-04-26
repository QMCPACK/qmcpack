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
#include "Estimators/FWEstimator.h"

namespace qmcplusplus
{

ForwardWalkingEstimator::ForwardWalkingEstimator(QMCHamiltonian& h)
  :refH(h)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  scalars.resize(2.0*(SizeOfHamiltonians+LE_MAX));
  scalars_saved.resize(2.0*(SizeOfHamiltonians+LE_MAX));
  //elocal_name is removed since refH is used
  //elocal_name.push_back("LocalEnergy");
  //elocal_name.push_back("LocalPotential");
  //for(int i=0; i<SizeOfHamiltonians; ++i)
  //  elocal_name.push_back(h.getObservableName(i));
}

ScalarEstimatorBase* ForwardWalkingEstimator::clone()
{
  return new ForwardWalkingEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void ForwardWalkingEstimator::add2Record(RecordListType& record)
{
  FirstIndex = record.size();
  record.add("LocalEnergy");
  record.add("LocalEnergy_sq");
  record.add("LocalPotential");
  record.add("LocalPotential_sq");
  for(int i=0; i<SizeOfHamiltonians; ++i)
  {
    record.add(refH.getObservableName(i));
    record.add(refH.getObservableName(i)+"_sq");
  }
  LastIndex=record.size();
  clear();
}

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3503 $   $Date: 2009-02-02 11:24:37 -0600 (Mon, 02 Feb 2009) $
 * $Id: ForwardWalkingEstimator.cpp 3503 2009-02-02 17:24:37Z jnkim $
 ***************************************************************************/
