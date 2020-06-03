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

#include "config.h"
#include "Message/OpenMP.h"
#ifdef HAVE_MKL
#include <mkl_service.h>
#endif

namespace qmcplusplus
{
/** service class for explicitly managing the threading of BLAS/LAPACK calls from OpenMP parallel region
 *
 * intended to use only locally around heavy calls.
 */
class BlasThreadingEnv
{
  int old_state;

public:
  /// Constructor, obtains the number of threads at the next level
  BlasThreadingEnv(int num_threads)
  {
#ifdef HAVE_MKL
    old_state = mkl_set_num_threads_local(num_threads);
#endif
  }

  ~BlasThreadingEnv()
  {
#ifdef HAVE_MKL
    mkl_set_num_threads_local(old_state);
#endif
  }

  static bool NestedThreadingSupported()
  {
#ifdef HAVE_MKL
    return true;
#else
    return false;
#endif
  }
};

} // namespace qmcplusplus
#endif
