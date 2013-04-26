//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

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
      cerr << "Unaligned reservation in PointerPool.  size = "
           << size << endl;
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
  vector<size_t> offsets;

};

#endif
