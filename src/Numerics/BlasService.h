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

#include "config.h"

namespace qmcplusplus
{

class BlasService
{
  int num_threads;
  int old_state;

  public:
  /// Constructor, obtains the number of threads at the next level
  BlasService();

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
