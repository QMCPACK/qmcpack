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
    
    



#include "Estimators/WFMCOnlyEstimator.h"

namespace qmcplusplus
{

WFMCOnlyEstimator::WFMCOnlyEstimator(QMCHamiltonian& h)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  elocal_name.push_back("LocalEnergy");
  elocal_name.push_back("LocalPotential");
  for(int i=0; i<SizeOfHamiltonians; ++i)
  {
    if (h.getObservableName(i)=="Psi")
      PsiIndex=i;
    elocal_name.push_back(h.getObservableName(i));
  }
  elocal_name.push_back("DMC_Psi_W");
  elocal_name.push_back("DMC_W");
  elocal_name.push_back("URW_LE");
  scalars.resize(elocal_name.size());
  scalars_saved.resize(elocal_name.size());
}

ScalarEstimatorBase* WFMCOnlyEstimator::clone()
{
  return new WFMCOnlyEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void WFMCOnlyEstimator::add2Record(RecordListType& record)
{
  FirstIndex = record.add(elocal_name[0].c_str());
  for(int i=1; i<elocal_name.size(); i++)
    record.add(elocal_name[i].c_str());
  LastIndex=FirstIndex + elocal_name.size();
  clear();
}

void WFMCOnlyEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}
