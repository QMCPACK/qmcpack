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
/** service class for explicitly managing the threading of BLAS/LAPACK calls from OpenMP parallel region
 */
class LAThreadProtector
{
  std::unique_ptr<Concurrency::ThreadCountProtector<>> handle_;
public:
  LAThreadProtector();
  ~LAThreadProtector();
};
} // namespace qmcplusplus
#endif
