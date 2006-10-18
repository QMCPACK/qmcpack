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

  LocalEnergyEstimator::LocalEnergyEstimator(QMCHamiltonian& h):Href(h) { 
    int hterms(h.size());
    SizeOfHamiltonians = hterms;
    FirstHamiltonian = h.startIndex();
    elocal.resize(SizeOfHamiltonians+LE_MAX);
    elocal_name.resize(SizeOfHamiltonians+LE_MAX);
    elocal_name[ENERGY_INDEX] = "LocalEnergy";
    elocal_name[ENERGY_SQ_INDEX] = "Variance";
    elocal_name[POTENTIAL_INDEX] = "LocalPotential";
    int ii(LE_MAX);
    for(int i=0; i < SizeOfHamiltonians; i++) elocal_name[ii++] = h.getName(i);
  }


  /** calculate the averages and reset to zero
   *\param record a container class for storing scalar records (name,value)
   *\param wgtinv the inverse weight
   */
  void LocalEnergyEstimator::report(RecordListType& record, RealType wgtinv, BufferType& msg) {
    register int ir=LocalEnergyIndex;
    b_average =  elocal[ENERGY_INDEX]*wgtinv;
    b_variance = elocal[ENERGY_SQ_INDEX]*wgtinv-b_average*b_average;
    record[ir++] = b_average;
    record[ir++] = b_variance;
    record[ir++] = elocal[POTENTIAL_INDEX]*wgtinv;
    for(int i=0, ii=LE_MAX; i<SizeOfHamiltonians; i++,ii++) {
      record[ir++] = elocal[ii]*wgtinv;
    }
    reset();
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
