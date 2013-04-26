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
#ifndef QMCPLUSPLUS_HF_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_HF_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"
#include "Particle/MCWalkerConfiguration.h"


namespace qmcplusplus
{

struct HFPolymerEstimator: public PolymerEstimator
{

  HFPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
  HFPolymerEstimator(const HFPolymerEstimator& mest);
  void setpNorm(RealType pn)
  {
    pNorm = pn;
  }

  void resize_HF(int pn, int nH)
  {
    nbds=pn;
    ObsSumL.resize(pn+1,nH);
    ObsSumR.resize(pn+1,nH);
  }
  void add_HF_Observables(int nchains, int strd);

//     inline void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker) { }

  /*@{*/
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt);
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/


  void evaluateDiff();

private:
  ///vector to contain the names of all the constituents of the local energy
  std::vector<string> elocal_name;
  int FirstHamiltonian,SizeOfHamiltonians;
  int stride,nobs,nbds;
  RealType KEconst,pNorm;
  Matrix<RealType> ObsSumL,ObsSumR;
  QMCHamiltonian& H;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: HFPolymerEstimator.h 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
