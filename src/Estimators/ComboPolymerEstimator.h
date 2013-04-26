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
#ifndef QMCPLUSPLUS_COMBO_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_COMBO_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"
#include "Particle/MCWalkerConfiguration.h"
#include "ReptationEstimators/ReptileEstimator.h"


namespace qmcplusplus
{

struct ComboPolymerEstimator: public PolymerEstimator
{

  ComboPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
  ComboPolymerEstimator(const ComboPolymerEstimator& mest);


  void put(xmlNodePtr cur, MCWalkerConfiguration& refWalker,int Rlength);

  /*@{*/
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt);
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/
  void evaluateDiff();

private:

  std::vector<string> scalars_name;
  std::vector<int> scalars_index;
  std::vector<ReptileEstimator*> RepEstimators;
  int FirstHamiltonian;
  int SizeOfHamiltonians;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: MJPolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
