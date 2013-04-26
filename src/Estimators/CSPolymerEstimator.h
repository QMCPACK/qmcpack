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
#ifndef QMCPLUSPLUS_CORRELATED_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_CORRELATED_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"

namespace qmcplusplus
{

struct CSPolymerEstimator: public PolymerEstimator
{

  CSPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
  //probably not needed
  CSPolymerEstimator(const CSPolymerEstimator& mest);

  /*@{*/
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt);
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: CSPolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
