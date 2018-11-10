//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file MemorySpace.hpp
 */
#ifndef QMCPLUSPLUS_MEMORYSPACE_H
#define QMCPLUSPLUS_MEMORYSPACE_H

namespace qmcplusplus
{
  struct MemorySpace
  {
    enum { HOST = 0, CUDA };
  };

  template<typename T, unsigned MemType = MemorySpace::HOST>
  class MemoryInstance
  {
    T* const inst;

    public:
    inline MemoryInstance(T* const in) : inst(in) { };
    inline T* getPointer() const { return inst; };
  };
}

#endif
