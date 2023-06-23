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


#include "BlasThreadingEnv.h"
#include "config.h"
#ifdef HAVE_MKL
#include <mkl_service.h>
#endif

namespace qmcplusplus
{
BlasThreadingEnv::BlasThreadingEnv(int num_threads)
{
#ifdef HAVE_MKL
  old_state_ = mkl_set_num_threads_local(num_threads);
#else
  old_state_ = 0;
#endif
}

BlasThreadingEnv::~BlasThreadingEnv()
{
#ifdef HAVE_MKL
  mkl_set_num_threads_local(old_state_);
#endif
}

bool BlasThreadingEnv::NestedThreadingSupported()
{
#ifdef HAVE_MKL
  return true;
#else
  return false;
#endif
}

} // namespace qmcplusplus
