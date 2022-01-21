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

#ifdef _OPENMP
#include <omp.h>
#else
using omp_int_t = int;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
inline omp_int_t omp_get_level() { return 0; }
inline omp_int_t omp_get_ancestor_thread_num(int level) { return 0; }
inline omp_int_t omp_get_max_active_levels() { return 1; }
inline void omp_set_num_threads(int num_threads) {}
#endif

/// get the number of threads at the next parallel level
inline int getNextLevelNumThreads()
{
  int num_threads = 1;
#pragma omp parallel
  {
#pragma omp master
    num_threads = omp_get_num_threads();
  }
  return num_threads;
}

#endif // OHMMS_COMMUNICATE_H
