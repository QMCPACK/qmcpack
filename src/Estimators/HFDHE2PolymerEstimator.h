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
    app_log()<<"  Truncation length set to: "<<truncLength<< std::endl;
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
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/


private:
  ///vector to contain the names of all the constituents of the local energy
  std::vector<std::string> elocal_name;
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
