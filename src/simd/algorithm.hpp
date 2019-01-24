//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file algorithm.hpp
 *
 * SIMD version of functions in algorithm
 */
#ifndef QMCPLUSPLUS_SIMD_ALGORITHM_HPP
#define QMCPLUSPLUS_SIMD_ALGORITHM_HPP

namespace qmcplusplus {

  namespace simd 
  {

    /** simd version of copy_n( InputIt first, Size count, OutputIt result)
     * @param first starting address of the input
     * @param count number of elements to copy
     * @param result starting address of the output
     * One hopes this optimization is worth it, caller is responsible for safety
     * first must of been aligned_alloc and this is still
     * technically a heap buffer read violation if first did not use all of the aligned block.
     * What is in those bytes is undefined.  So if your result isn't aligned_alloc and 
     * doesn't get a matching element count from somewhere you are going to have garbage in result.
     */
    template<typename T1, typename T2>
    __attribute__((no_sanitize_address)) inline void copy_n(const T1* restrict first, size_t count, T2* restrict result) 
      {
//#pragma omp simd

	//If anything is going to emit SIMD instructions
	//std::memcpy should
	std::memcpy(result, first, count * sizeof(T1));
      }

    template<typename T1, typename T2>
    inline T2 accumulate_n(const T1* restrict in, size_t n, T2 res)
      {
#pragma omp simd reduction(+:res)
        for(int i=0; i<n; ++i)
          res += in[i];
        return res;
      }

    ///inner product
    template<typename T1, typename T2, typename T3>
      inline T3 inner_product_n(const T1* restrict a, const T2* restrict b, int n, T3 res)
      {
        for(int i=0; i<n; ++i) res += a[i]*b[i];
        return res;
      }

  } //simd namepsace
}
#endif
