//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef POINTER_POOL_H
#define POINTER_POOL_H

#include <vector>

template<typename CONT>
class PointerPool
{
public:
  typedef CONT buffer_type;
  typedef typename CONT::pointer pointer;

  // Local data routines
  pointer getPointer (int index, buffer_type &buffer)
  {
    return &(buffer[offsets[index]]);
  }

  void allocate (buffer_type &buffer)
  {
    buffer.resize(totalSize);
  }

  // Shared-data routines

  // Reserves size elements and returns the offset to the member
  // in the buffer
  size_t reserve (size_t size)
  {
    if (size % 32)
    {
      std::cerr << "Unaligned reservation in PointerPool.  size = "
           << size << std::endl;
      size += 32 - (size % 32);
    }
    size_t off = totalSize;
    offsets.push_back(off);
    totalSize += size;
    return off;
  }

  void reset()
  {
    offsets.resize(0);
    totalSize = 0;
  }

  size_t getTotalSize()
  {
    return totalSize;
  }

  PointerPool() : totalSize(0)
  { }

protected:
  size_t totalSize;
  std::vector<size_t> offsets;

};

#endif
