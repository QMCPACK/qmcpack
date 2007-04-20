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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Estimators/LocalEnergyEstimator.h"

namespace qmcplusplus {

  LocalEnergyEstimator::LocalEnergyEstimator(QMCHamiltonian& h)
  { 
    int hterms(h.size());
    SizeOfHamiltonians = hterms;
    FirstHamiltonian = h.startIndex();
    d_data.resize((SizeOfHamiltonians+LE_MAX)*2);
    elocal_name.resize((SizeOfHamiltonians+LE_MAX)*2);

    int target=0;
    elocal_name[target++] = "LocalEnergy";
    elocal_name[target++] = "LocalEnergy2";
    elocal_name[target++] = "LocalPotential";
    elocal_name[target++] = "LocalPotential2";
    for(int i=0; i < SizeOfHamiltonians; i++) 
    {
      elocal_name[target++] = h.getName(i);
      elocal_name[target++] = h.getName(i)+"2";
    }
  }

  LocalEnergyEstimator::LocalEnergyEstimator(const LocalEnergyEstimator& est):
    ScalarEstimatorBase(est),
  SizeOfHamiltonians(est.SizeOfHamiltonians),
  FirstHamiltonian(est.FirstHamiltonian),
  elocal_name(est.elocal_name)
  {
  }

  ScalarEstimatorBase* LocalEnergyEstimator::clone()
  {
    return new LocalEnergyEstimator(*this);
  }

  /** calculate the averages and reset to zero
   *\param record a container class for storing scalar records (name,value)
   *\param wgtinv the inverse weight
   */
  void LocalEnergyEstimator::report(RecordListType& record, RealType wgtinv, BufferType& msg) {
    msg.get(d_data.begin(),d_data.end());
    report(record,wgtinv);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
