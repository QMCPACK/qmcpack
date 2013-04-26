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

void WFMCOnlyEstimator::registerObservables(vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3449 $   $Date: 2009-01-11 08:47:20 -0600 (Sun, 11 Jan 2009) $
 * $Id: WFMCOnlyEstimator.cpp 3449 2009-01-11 14:47:20Z jnkim $
 ***************************************************************************/
