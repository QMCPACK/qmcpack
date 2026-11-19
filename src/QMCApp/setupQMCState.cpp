//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "setupQMCState.h"
#include <sstream>
#include "CommandLineOptions.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "Utilities/qmc_common.h"

namespace qmcplusplus
{

void setupQMCState(const CommandLineOptions& options, const Communicate& comm, QMCState& qmc_state)
{
  qmc_state.dryrun  = options.get_dryrun();
  qmc_state.io_node = (comm.rank() == 0 && !options.get_noprint());

  const int max_vnode_ranks = comm.getVirtualNodeRankCount();
  const int vnode_ranks     = options.get_vnode_ranks().value_or(max_vnode_ranks);
  if (max_vnode_ranks % vnode_ranks != 0)
  {
    std::ostringstream msg;
    msg << "setupQMCState(). --vnode-ranks (" << vnode_ranks
        << ") can only be a factor of the uniform node rank count (" << max_vnode_ranks << ").";
    throw UniformCommunicateError(msg.str());
  }

  const auto& inputs = options.get_inputs();
  if (inputs.size() > 1)
  {
    const int total_blocks = comm.size() / vnode_ranks;
    if (total_blocks < inputs.size())
    {
      std::ostringstream msg;
      msg << "setupQMCState(). Current " << comm.size() << " MPI ranks (" << total_blocks
          << " blocks of size " << vnode_ranks << ") cannot accommodate all the "
          << inputs.size() << " individual calculations in the ensemble. "
          << "Increase the number of MPI ranks or reduce the number of calculations.";
      throw UniformCommunicateError(msg.str());
    }
  }
}

} // namespace qmcplusplus
