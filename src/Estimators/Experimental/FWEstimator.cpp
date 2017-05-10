//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    



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
