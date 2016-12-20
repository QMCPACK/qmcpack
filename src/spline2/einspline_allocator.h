/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//                                                                         //
//  This program is free software; you can redistribute it and/or modify   //
//  it under the terms of the GNU General Public License as published by   //
//  the Free Software Foundation; either version 2 of the License, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  This program is distributed in the hope that it will be useful,        //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//  GNU General Public License for more details.                           //
//                                                                         //
//  You should have received a copy of the GNU General Public License      //
//  along with this program; if not, write to the Free Software            //
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,                     //
//  Boston, MA  02110-1301  USA                                            //
/////////////////////////////////////////////////////////////////////////////
/** @file einspline_allocator.h
 *
 * Rename aligned_alloc/aligned_free as einspline_alloc/einspline_free to
 * avoid naming conflicts with the standards
 */
#ifndef EINSPLINE_ALIGNED_ALLOC_H
#define EINSPLINE_ALIGNED_ALLOC_H

#include <stdlib.h>
#include "config.h"
#include <omp.h>

#if defined(__INTEL_COMPILER)
inline void *
einspline_alloc (size_t size, size_t alignment)
{
  return _mm_malloc(size,alignment);
}

inline void
einspline_free (void *ptr)
{
  _mm_free(ptr);
}
#else

#ifdef HAVE_POSIX_MEMALIGN

#if (defined(__IBMCPP__)) && ( __IBMCPP__ <= 1210 )
#else
int posix_memalign(void **memptr, size_t alignment, size_t size);
#endif

inline void *
einspline_alloc (size_t size, size_t alignment)
{
  void *ptr;
  posix_memalign (&ptr, alignment, size);
  return ptr;
}

inline void
einspline_free (void *ptr)
{
  free (ptr);
}

#else

inline void *
einspline_alloc (size_t size, size_t alignment)
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
einspline_free (void *aligned)
{
  void *ptr = *((void**)aligned-1);
  free (ptr);
}
#endif


#endif
#endif
