//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file Mallocator.hpp
 */

#include "spline2/einspline_allocator.h"

namespace qmcplusplus
{

template <class T>
struct Mallocator {
  typedef T value_type;
  Mallocator() = default;
  template <class U> Mallocator(const Mallocator<U>&) {}
  T* allocate(std::size_t n) { return static_cast<T*>(einspline_alloc(n*sizeof(T), QMC_CLINE)); }
  void deallocate(T* p, std::size_t) { einspline_free(p); }
};

}
