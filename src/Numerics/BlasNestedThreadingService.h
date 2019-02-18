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

#ifndef QMCPLUSPLUS_BLAS_SERVICE_H
#define QMCPLUSPLUS_BLAS_SERVICE_H

namespace qmcplusplus
{

/** service class for managing threaded BLAS/LAPACK calls from OpenMP parallel region
 *
 * intended to use only locally around heavy calls.
 * omp_get_nested() protects every call.
 */
class BlasNestedThreadingService
{
  int num_threads;
  int old_state;

  public:
  /// Constructor, obtains the number of threads at the next level
  BlasNestedThreadingService();

  /// get num_threads
  int getNumThreads() { return num_threads; }

  /// Detect the threading environment and preset BLAS libraray
  void presetBLASNumThreads();

  /// Detect the threading environment and unset BLAS libraray
  void unsetBLASNumThreads();

  /// Return true, if nested threading is supported by the BLAS/LAPACK
  bool NestedThreadingSupported();
};

}
#endif
