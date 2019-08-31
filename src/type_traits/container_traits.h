//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Lab
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONTAINER_TRAITS_H
#define QMCPLUSPLUS_CONTAINER_TRAITS_H

#include <sstream>

namespace qmcplusplus
{

namespace container_traits
{

template<typename CT, typename I>
inline void resize(CT& ref, I* n, int d)
{
  throw std::runtime_error("Unknown container, resizing is not available!");
}

template<typename CT>
inline size_t getSize(CT& ref)
{
  return ref.size();
}

// tempalte specialization for resize
template<typename CT, typename I>
inline void resize1D(CT& ref, I* n, int d)
{
  size_t nt = 1;
  for (int i = 0; i < d; ++i)
    nt *= n[i];
  ref.resize(nt);
}

template<typename T, class ALLOC, typename I>
inline void resize1D(std::vector<T, ALLOC>& ref, I* n, int d)
{
  resize1D(ref, n, d);
}

template<typename T, class ALLOC, typename I>
inline void resize1D(Vector<T, ALLOC>& ref, I* n, int d)
{
  resize1D(ref, n, d);
}

template<typename T, class ALLOC, typename I>
inline void resize(Matrix<T, ALLOC>& ref, I* n, int d)
{
  if (d != 2)
  {
    std::ostringstream err_msg;
    err_msg << "Matrix cannot be resized. Requested dimension = " << d << std::endl;
    throw std::runtime_error(err_msg.str());
  }
  ref.resize(n[0], n[1]);
}

template<typename T, unsigned D, typename I>
inline static void resize(Array<T, D>& ref, I* n, int d)
{
  if (d != D)
  {
    std::ostringstream err_msg;
    err_msg << "Array<T, " << D << "> cannot be resized. Requested dimension = " << d << std::endl;
    throw std::runtime_error(err_msg.str());
  }
  ref.resize(n);
}

} // namespace container_traits
} // namespace qmcplusplus
#endif
