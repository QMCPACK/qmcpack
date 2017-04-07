//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




#ifndef OHMMS_OPENMP_H
#define OHMMS_OPENMP_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(ENABLE_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num()
{
  return 0;
}
inline omp_int_t omp_get_max_threads()
{
  return 1;
}
inline omp_int_t omp_get_num_threads()
{
  return 1;
}
#endif
#endif // OHMMS_COMMUNICATE_H
