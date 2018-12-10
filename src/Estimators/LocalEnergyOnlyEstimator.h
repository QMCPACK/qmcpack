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
    
    



#ifndef QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGY_ONLY_ESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
namespace qmcplusplus
{

/** Estimator for local energy only
 */
struct LocalEnergyOnlyEstimator: public ScalarEstimatorBase
{

  inline LocalEnergyOnlyEstimator()
  {
    scalars.resize(2);
    scalars_saved.resize(2);
  }

  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first != last; ++first)
    {
      scalars[0]((*first)->Properties(LOCALENERGY),wgt);
      scalars[1]((*first)->Properties(LOCALPOTENTIAL),wgt);
    }
  }

  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
  {}
  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   * @param record storage of scalar records (name,value)
   */
  inline void add2Record(RecordListType& record)
  {
    FirstIndex = record.add("LocalEnergy");
    int s1=record.add("LocalPotential");
    LastIndex = FirstIndex+2;
    // int s2=record.add("KineticEnergy");
    //LastIndex = FirstIndex+3;
    clear();
  }
  ScalarEstimatorBase* clone()
  {
    return new LocalEnergyOnlyEstimator();
  }

};
}
#endif
