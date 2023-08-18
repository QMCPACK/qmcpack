//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_LATHREADPROTECTOR_H
#define QMCPLUSPLUS_LATHREADPROTECTOR_H

#include <memory>
#include <Concurrency/UtilityFunctions.hpp>

namespace qmcplusplus
{
/** For linear algebra only. A service class to restore active avaiable threads upon destruction as the thread count recorded during construction
 * It protects any side effects from linear algebra library calls changing the number of active avaiable threads.
 * Known trouble maker: OpenBLAS https://github.com/xianyi/OpenBLAS/issues/3940
 */
class OMPThreadCountProtectorLA
{
  using Protector = Concurrency::ThreadCountProtector<Executor::OPENMP>;
  std::unique_ptr<Protector> handle_;

public:
  OMPThreadCountProtectorLA();
  ~OMPThreadCountProtectorLA();
};
} // namespace qmcplusplus
#endif
