//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FAKEESTIMATOR_H
#define QMCPLUSPLUS_FAKEESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"
namespace qmcplusplus
{
class FakeEstimator : public ScalarEstimatorBase
{
public:
  std::string getName() const override { return "FakeEstimator"; }

  void accumulate(const MCWalkerConfiguration& W, WalkerIterator first, WalkerIterator last, RealType wgt) override {}

  void accumulate(const RefVector<MCPWalker>& walkers) override {}

  void add2Record(RecordNamedProperty<RealType>& record) override {}

  void registerObservables(std::vector<ObservableHelper>& h5dec, hdf_archive& file) override {}

  FakeEstimator* clone() override { return new FakeEstimator; }
  
  std::string type_{"fake"};
  const std::string& getSubTypeStr() const override { return type_; }
};

} // namespace qmcplusplus
#endif
