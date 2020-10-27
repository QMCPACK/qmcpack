//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ACCELERATORS_HPP
#define QMCPLUSPLUS_ACCELERATORS_HPP

#include "Message/Communicate.h"

namespace qmcplusplus
{
/* assign one default accelerator to each MPI rank
 * @param NodeComm an MPI communicator including all the ranks within a shared memory node
 * @return the number of accelerators within the node
 *
 * on a node with multiple accelerators, each MPI rank will be assgined a default accelerator.
 * If the number of MPI ranks is more than the number of accelerators,
 * MPI ranks will distributed among the accelerators.
 * This routine also initialze the device runtime library.
 */
int assignAccelerators(Communicate& NodeComm);
} // namespace qmcplusplus
#endif
