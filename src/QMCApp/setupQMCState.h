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

#ifndef QMCPLUSPLUS_SETUP_QMC_STATE_H
#define QMCPLUSPLUS_SETUP_QMC_STATE_H

class Communicate;

namespace qmcplusplus
{

class CommandLineOptions;
struct QMCState;

/** Setup QMCState based on parsed command line options and communication topology.
 * @param options parsed command line options
 * @param comm communicator (e.g. world communicator)
 * @param qmc_state global QMCState object to configure
 */
void setupQMCState(const CommandLineOptions& options, const Communicate& comm, QMCState& qmc_state);

} // namespace qmcplusplus

#endif
