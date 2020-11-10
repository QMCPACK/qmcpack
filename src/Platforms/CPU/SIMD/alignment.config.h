//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef ALIGNMENT_CONFIG_H
#define ALIGNMENT_CONFIG_H

#if defined(__INTEL_COMPILER)
  #if defined(__AVX512F__)
    #define QMC_CLINE 64
    #define ASSUME_ALIGNED(x) __assume_aligned(x,64)
  #else
    #define QMC_CLINE 32
    #define ASSUME_ALIGNED(x) __assume_aligned(x,32)
  #endif
#elif defined(__GNUC__) && !defined(__ibmxl__)
  #if defined(__AVX512F__)
    #define QMC_CLINE 64
    #define ASSUME_ALIGNED(x) (x) = (__typeof__(x)) __builtin_assume_aligned(x,64)
  #else
    #define QMC_CLINE 32
    #define ASSUME_ALIGNED(x) (x) = (__typeof__(x)) __builtin_assume_aligned(x,32)
  #endif
#else
  #define QMC_CLINE 32
  #define ASSUME_ALIGNED(x)
#endif

#endif
