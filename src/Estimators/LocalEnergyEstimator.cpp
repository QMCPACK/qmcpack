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
#include "Estimators/LocalEnergyEstimator.h"

namespace qmcplusplus {

  LocalEnergyEstimator::LocalEnergyEstimator(QMCHamiltonian& h)
    :refH(h)
  { 
    SizeOfHamiltonians = h.sizeOfObservables();
    FirstHamiltonian = h.startIndex();
    scalars.resize(SizeOfHamiltonians+LE_MAX);
    scalars_saved.resize(SizeOfHamiltonians+LE_MAX);
    //elocal_name is removed since refH is used
    //elocal_name.push_back("LocalEnergy");
    //elocal_name.push_back("LocalPotential");
    //for(int i=0; i<SizeOfHamiltonians; ++i)
    //  elocal_name.push_back(h.getObservableName(i));
  }

  ScalarEstimatorBase* LocalEnergyEstimator::clone()
  {
    return new LocalEnergyEstimator(*this);
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
    for(int i=0; i<SizeOfHamiltonians; ++i) record.add(refH.getObservableName(i));
    LastIndex=record.size();
    clear();
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
