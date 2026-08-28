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


#include <cstddef>
#include <atomic>

namespace qmcplusplus
{

class MemoryUsageAccount
{
  std::atomic<size_t> mem_allocated{0};

public:
  void creditUsage(size_t n) { mem_allocated += n; }
  void debitUsage(size_t n) { mem_allocated -= n; }
  size_t getBalance() const { return mem_allocated; }
};
} // namespace qmcplusplus
