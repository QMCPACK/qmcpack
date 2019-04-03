//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OMPstd.hpp
 */

#ifndef QMCPLUSPLUS_OPENMP_STD_H
#define QMCPLUSPLUS_OPENMP_STD_H

namespace OMPstd
{
template<typename T>
inline void fill_n(T* x, size_t count, const T& value)
{
  PRAGMA_OFFLOAD("omp for")
  for (size_t id = 0; id < count; id++)
    x[id] = value;
}
} // namespace OMPstd
#endif
