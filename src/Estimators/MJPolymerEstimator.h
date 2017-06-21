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
    
    



#ifndef QMCPLUSPLUS_FULL_POLYMERESTIMATOR_H
#define QMCPLUSPLUS_FULL_POLYMERESTIMATOR_H

#include "Estimators/PolymerEstimator.h"
#include "Particle/MCWalkerConfiguration.h"


namespace qmcplusplus
{

struct MJPolymerEstimator: public PolymerEstimator
{

  MJPolymerEstimator(QMCHamiltonian& h, int hcopy=1, MultiChain* polymer=0);
  MJPolymerEstimator(const MJPolymerEstimator& mest);
  void setpNorm(RealType pn)
  {
    pNorm = pn;
  }

  void setrLen(int pn)
  {
    ObsCont.resize(pn,0.0);
    ObsContAvg.resize(pn,0.0);
    ObsCont2.resize(pn,0.0);
    ObsContAvg2.resize(pn,0.0);
    ObsEvals=0;
    truncLength.resize(3);
    truncLength[0]= static_cast<int>(0.25*pn);
    truncLength[1]= static_cast<int>(0.5*pn);
    truncLength[2]= static_cast<int>(0.75*pn);
  }

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
  int FirstHamiltonian;
  int SizeOfHamiltonians;
  std::vector<int> truncLength;
  int Findex;
  RealType KEconst,pNorm,ObsEvals;
  std::vector<RealType> ObsCont,ObsContAvg;
  std::vector<RealType> ObsCont2,ObsContAvg2;
//       QMCHamiltonian* Hpointer;
};

}
#endif
