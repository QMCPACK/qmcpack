//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_VECTOR_VIEWER_H
#define QMCPLUSPLUS_VECTOR_VIEWER_H

namespace qmcplusplus
{

/** a class to map a memory sequence to a vector
 * @tparam T datatype
 */
template<typename T>
struct VectorViewer
{
#if (__cplusplus >= 201103L)
  using value_type=T;
#else
  typedef T value_type;
#endif

  T* data_ptr;
  size_t data_size;

  inline VectorViewer(T* a, size_t n):data_ptr(a),data_size(n) {}

  inline T* data()
  {
    return data_ptr;
  }

  inline size_t size() const
  {
    return data_size;
  }

  inline T& operator[](size_t i)
  {
    return data_ptr[i];
  }

  inline T operator[](size_t i) const
  {
    return data_ptr[i];
  }
};

}

#endif
