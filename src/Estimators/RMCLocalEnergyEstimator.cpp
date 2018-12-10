//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
    std::ostringstream ss;
    ss << refH.getObservableName(i)<<"_m";
    record.add(ss.str());
  }
  for(int i=0; i < SizeOfHamiltonians; i++)
  {
    std::ostringstream ss;
    ss << refH.getObservableName(i)<<"_p";
    // app_log()<<"Registering observable "<<ss.str()<< std::endl;
    record.add(ss.str());
  }
  LastIndex=record.size();
  clear();
}

}
