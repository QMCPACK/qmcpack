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
struct LocalEnergyOnlyEstimator : public ScalarEstimatorBase
{
  using WP = WalkerProperties::Indexes;

  inline LocalEnergyOnlyEstimator()
  {
    scalars.resize(2);
    scalars_saved.resize(2);
  }

  std::string getName() const override { return "LocalEnergyOnlyEstimator"; }
  
  inline void accumulate(const MCWalkerConfiguration& W,
                         WalkerIterator first,
                         WalkerIterator last,
                         RealType wgt) override
  {
    for (; first != last; ++first)
    {
      scalars[0]((*first)->Properties(WP::LOCALENERGY), wgt);
      scalars[1]((*first)->Properties(WP::LOCALPOTENTIAL), wgt);
    }
  }

  inline void accumulate(const RefVector<MCPWalker>& walkers) override
  {
    for (MCPWalker& walker : walkers)
    {
      scalars[0](walker.Properties(WP::LOCALENERGY), 1.0);
      scalars[1](walker.Properties(WP::LOCALPOTENTIAL), 1.0);
    }
  }

  void registerObservables(std::vector<ObservableHelper>& h5dec, hdf_archive& file) override {}

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   * @param record storage of scalar records (name,value)
   */
  inline void add2Record(RecordListType& record) override
  {
    FirstIndex = record.add("LocalEnergy");
    int s1     = record.add("LocalPotential");
    LastIndex  = FirstIndex + 2;
    // int s2=record.add("KineticEnergy");
    //LastIndex = FirstIndex+3;
    clear();
  }

  LocalEnergyOnlyEstimator* clone() override { return new LocalEnergyOnlyEstimator(); }

  const std::string type_str = "LocalEnergyOnlyEstimatorNotSupportedInBatchedVersion";
  const std::string& getSubTypeStr() const override { return type_str; }
};
} // namespace qmcplusplus
#endif
