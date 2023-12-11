//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef BUFFER_MANAGERS_H
#define BUFFER_MANAGERS_H

#include "AFQMC/Memory/host_buffer_manager.hpp"
#include "AFQMC/Memory/device_buffer_manager.hpp"
#include "AFQMC/Memory/localTG_buffer_manager.hpp"

namespace qmcplusplus
{
namespace afqmc
{
void setup_memory_managers(mpi3::shared_communicator& local, size_t size);
void setup_memory_managers(mpi3::shared_communicator& node, size_t size, int nc);
void update_memory_managers();
void release_memory_managers();

} // namespace afqmc
} // namespace qmcplusplus

#endif
