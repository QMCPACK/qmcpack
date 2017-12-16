//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/


  void evaluateDiff();

private:
  ///vector to contain the names of all the constituents of the local energy
  std::vector<std::string> elocal_name;
  int FirstHamiltonian,SizeOfHamiltonians;
  int stride,nobs,nbds;
  RealType KEconst,pNorm;
  Matrix<RealType> ObsSumL,ObsSumR;
  QMCHamiltonian& H;
};

}
#endif
