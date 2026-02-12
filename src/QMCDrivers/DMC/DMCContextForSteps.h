//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <NonLocalTOperator.h>

namespace qmcplusplus
{
class DMCBatched::DMCContextForSteps : public ContextForSteps
{
public:
  DMCContextForSteps(RandomBase<FullPrecRealType>& random_gen, NonLocalTOperator&& non_local_ops)
      : ContextForSteps(random_gen), non_local_ops(non_local_ops)
  {}

  ///non local operator
  NonLocalTOperator non_local_ops;
};
} // namespace qmcplusplus
