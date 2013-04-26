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
#ifndef QMCPLUSPLUS_HFDHE2_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_HFDHE2_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"

namespace qmcplusplus
{

struct HFDHE2PolymerEstimator: public PolymerEstimator
{

  HFDHE2PolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
  HFDHE2PolymerEstimator(const HFDHE2PolymerEstimator& mest);
  void setpNorm(RealType pn)
  {
    pNorm = pn;
  }
  void settruncLength(int pn)
  {
    truncLength = pn;
    app_log()<<"  Truncation length set to: "<<truncLength<<endl;
  }
  void setrLen(int pn)
  {
    ObsCont.resize(pn,0.0);
    ObsContAvg.resize(pn,0.0);
    ObsCont2.resize(pn,0.0);
    ObsContAvg2.resize(pn,0.0);
  }

  /*@{*/
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt);
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/


private:
  ///vector to contain the names of all the constituents of the local energy
  std::vector<string> elocal_name;
  std::vector<RealType> ObsCont,ObsContAvg;
  std::vector<RealType> ObsCont2,ObsContAvg2;
  int FirstHamiltonian;
  int SizeOfHamiltonians;
  RealType KEconst, pNorm, ObsEvals ;
  int HDFHE2index, Pindex, truncLength;
//       QMCHamiltonian* Hpointer;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: HFDHE2PolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
