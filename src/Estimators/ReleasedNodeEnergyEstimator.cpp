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
#include "Estimators/ReleasedNodeEnergyEstimator.h"

namespace qmcplusplus {

  ReleasedNodeEnergyEstimator::ReleasedNodeEnergyEstimator(QMCHamiltonian& h)
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

  ScalarEstimatorBase* ReleasedNodeEnergyEstimator::clone()
  {
    return new ReleasedNodeEnergyEstimator(*this);
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   * @param record storage of scalar records (name,value)
   */
  void ReleasedNodeEnergyEstimator::add2Record(RecordListType& record) 
  {
    FirstIndex = record.size();
    int dumy=record.add("LocalEnergy");
    dumy=record.add("LocalEnergy_sq");
    dumy=record.add("LocalPotential");
    dumy=record.add("FermionEnergy");
    dumy=record.add("AverageSign");
    for(int i=0; i<SizeOfHamiltonians; ++i) record.add(refH.getObservableName(i));
    LastIndex=record.size();
    clear();
  }

}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4163 $   $Date: 2009-08-31 05:47:46 -0500 (Mon, 31 Aug 2009) $
 * $Id: ReleasedNodeEnergyEstimator.cpp 4163 2009-08-31 10:47:46Z jmcminis $ 
 ***************************************************************************/
