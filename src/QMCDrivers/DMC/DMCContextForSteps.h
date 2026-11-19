//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <NonLocalTOperator.h>
#include "L2Diffusion.h"

namespace qmcplusplus
{
class DMCBatched::DMCContextForSteps : public ContextForSteps
{
public:
  DMCContextForSteps(RandomBase<FullPrecRealType>& random_gen, NonLocalTOperator&& non_local_ops, bool use_l2_diffusion)
      : ContextForSteps(random_gen),
        non_local_ops(non_local_ops),
        l2_workspace(use_l2_diffusion ? std::make_unique<L2Diffusion::Workspace>() : nullptr)
  {}

  ///non local operator
  NonLocalTOperator non_local_ops;
  /// Mutable L2 scratch data is crowd-local because crowd tasks run concurrently.
  std::unique_ptr<L2Diffusion::Workspace> l2_workspace;
};
} // namespace qmcplusplus
