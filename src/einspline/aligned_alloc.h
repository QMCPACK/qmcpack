//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef ALIGNED_ALLOC_H
#define ALIGNED_ALLOC_H

#include <stdlib.h>
#include "config.h"

#ifdef HAVE_POSIX_MEMALIGN
inline void *
aligned_alloc (size_t size, size_t alignment)
{
  void *ptr;
  posix_memalign (&ptr, alignment, size);
  return ptr;
}

inline void
aligned_free (void *ptr)
{
  free (ptr);
}

#else

inline void *
aligned_alloc (size_t size, size_t alignment)
{
  size += (alignment-1)+sizeof(void*);
  void *ptr = malloc (size);
  if (ptr == NULL)
    return NULL;
  else
  {
    void *shifted = ptr + sizeof(void*);
    size_t offset = alignment - (size_t)shifted%(size_t)alignment;
    void *aligned = shifted + offset;
    *((void**)aligned-1) = ptr;
    return aligned;
  }
}

inline void
aligned_free (void *aligned)
{
  void *ptr = *((void**)aligned-1);
  free (ptr);
}
#endif


#endif
