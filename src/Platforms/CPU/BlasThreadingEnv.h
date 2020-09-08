//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_BLAS_THREADING_ENV_H
#define QMCPLUSPLUS_BLAS_THREADING_ENV_H

namespace qmcplusplus
{
/** service class for explicitly managing the threading of BLAS/LAPACK calls from OpenMP parallel region
 *
 * intended to use only locally around heavy calls.
 */
class BlasThreadingEnv
{
  int old_state_;

public:
  /// Constructor, obtains the number of threads at the next level
  BlasThreadingEnv(int num_threads);
  ~BlasThreadingEnv();

  static bool NestedThreadingSupported();
};
} // namespace qmcplusplus
#endif
