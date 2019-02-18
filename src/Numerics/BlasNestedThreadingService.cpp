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


#include "Numerics/BlasNestedThreadingService.h"
#include <omp.h>
#include "config.h"
#ifdef HAVE_MKL
#include<mkl_service.h>
#endif

namespace qmcplusplus
{

BlasNestedThreadingService::BlasNestedThreadingService() : num_threads(1), old_state(0)
{
  if(omp_get_nested())
  {
    #pragma omp parallel
    {
      #pragma omp master
      num_threads = omp_get_num_threads();
    }
  }
}

void BlasNestedThreadingService::presetBLASNumThreads()
{
  if(omp_get_nested() && num_threads>1)
  {
#ifdef HAVE_MKL
    old_state = mkl_set_num_threads_local(num_threads);
#endif
  }
}

void BlasNestedThreadingService::unsetBLASNumThreads()
{
  if(omp_get_nested() && num_threads>1)
  {
#ifdef HAVE_MKL
    mkl_set_num_threads_local(old_state);
#endif
  }
}

bool BlasNestedThreadingService::NestedThreadingSupported()
{
#ifdef HAVE_MKL
  return true;
#else
  return false;
#endif
}

}
