//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
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
  void accumulate(const MCWalkerConfiguration& W, WalkerIterator first, WalkerIterator last, RealType wgt) override {}

  void accumulate(const RefVector<MCPWalker>& walkers) override {}

  void add2Record(RecordNamedProperty<RealType>& record) override {}

  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid) override {}

  ScalarEstimatorBase* clone() override { return new FakeEstimator; }
};

} // namespace qmcplusplus
#endif
