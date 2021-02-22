//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MemoryUsage.h"
#include <cstring>
#include <string>
#include "Host/sysutil.h"
#include "OMPTarget/OMPallocator.hpp"
#ifdef ENABLE_CUDA
#include "CUDA/CUDAallocator.hpp"
#include "CUDA/CUDAruntime.h"
#endif

namespace qmcplusplus
{

void print_mem(const char* title, std::ostream& log)
{
  char msg[256];
  std::string line_separator;
  for (int i = 0; i < strlen(title) + 30; i++)
    line_separator += "#";
  log << line_separator << std::endl;
  sprintf(msg, "### Memory usage report : %s ###\n", title);
  log << msg;
  log << line_separator << std::endl;
  sprintf(msg, "Available memory on node 0, free + buffers : %7zu MiB\n", freemem() >> 20);
  log << msg;
  sprintf(msg, "Memory footprint by rank 0 on node 0       : %7zu MiB\n", memusage() >> 10);
  log << msg;
#ifdef ENABLE_CUDA
  sprintf(msg, "Device memory allocated via CUDA allocator : %7zu MiB\n", getCUDAdeviceMemAllocated() >> 20);
  log << msg;
  sprintf(msg, "Free memory available on default device    : %7zu MiB\n", getCUDAdeviceFreeMem() >> 20);
  log << msg;
#endif
#ifdef ENABLE_OFFLOAD
  sprintf(msg, "Device memory allocated via OpenMP offload : %7zu MiB\n", getOMPdeviceMemAllocated() >> 20);
  log << msg;
#endif
  log << line_separator << std::endl;
}

}
