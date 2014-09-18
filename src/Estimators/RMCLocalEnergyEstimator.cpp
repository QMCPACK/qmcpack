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
#include "Estimators/RMCLocalEnergyEstimator.h"

namespace qmcplusplus
{

RMCLocalEnergyEstimator::RMCLocalEnergyEstimator(QMCHamiltonian& h, int nobs)
  :refH(h), NObs(nobs)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  RMCSpecificTerms=8;
  scalars.resize(2*SizeOfHamiltonians+RMCSpecificTerms);
  scalars_saved.resize(2*SizeOfHamiltonians+RMCSpecificTerms);
}

ScalarEstimatorBase* RMCLocalEnergyEstimator::clone()
{
  return new RMCLocalEnergyEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void RMCLocalEnergyEstimator::add2Record(RecordListType& record)
{
  FirstIndex = record.size();
  record.add("LocalEnergy");
  record.add("LocalEnergy_sq");
  record.add("LocalEnergy_p");
  record.add("LocalEnergy_sq_p");
  record.add("LocalEnergy_sq_cross");
  record.add("LocalPotential");
  record.add("LocalPotential_pure");
  
  record.add("OldestBead");
  //for(int j=0; j <= NObs; j++)
  for(int i=0; i < SizeOfHamiltonians; i++)
  {
    ostringstream ss;
    ss << refH.getObservableName(i)<<"_m";
    record.add(ss.str());
  }
  for(int i=0; i < SizeOfHamiltonians; i++)
  {
    ostringstream ss;
    ss << refH.getObservableName(i)<<"_p";
    // app_log()<<"Registering observable "<<ss.str()<<endl;
    record.add(ss.str());
  }
  LastIndex=record.size();
  clear();
}

}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4165 $   $Date: 2009-08-31 05:47:46 -0500 (Mon, 31 Aug 2009) $
 * $Id: RMCLocalEnergyEstimator.cpp 4165 2009-08-31 10:47:46Z jmcminis $
 ***************************************************************************/
