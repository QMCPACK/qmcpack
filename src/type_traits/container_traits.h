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
template<typename CT>
struct container_traits
{
  static constexpr int DIM = -1;
  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    throw std::runtime_error("Unknown container, resizing is not available!");
  }
};

template<typename T, class ALLOC>
struct container_traits<std::vector<T, ALLOC>>
{
  static constexpr int DIM = 1;
  template<typename I>
  inline static void resize(std::vector<T, ALLOC>& ref, I* n, int d)
  {
    size_t nt = 1;
    for (int i = 0; i < d; ++i)
      nt *= n[i];
    ref.resize(nt);
  }
};

template<typename T, class ALLOC>
struct container_traits<Vector<T, ALLOC>>
{
  static constexpr int DIM = 1;
  template<typename I>
  inline static void resize(Vector<T, ALLOC>& ref, I* n, int d)
  {
    size_t nt = 1;
    for (int i = 0; i < d; ++i)
      nt *= n[i];
    ref.resize(nt);
  }
};

template<typename T, class ALLOC>
struct container_traits<Matrix<T, ALLOC>>
{
  static constexpr int DIM = 2;
  template<typename I>
  inline static void resize(Matrix<T, ALLOC>& ref, I* n, int d)
  {
    if (d != 2)
    {
      std::ostringstream err_msg;
      err_msg << "Matrix cannot be resized. Requested dimension = " << d << std::endl;
      throw std::runtime_error(err_msg.str());
    }
    ref.resize(n[0], n[1]);
  }
};

template<typename T, unsigned D>
struct container_traits<Array<T, D>>
{
  static constexpr int DIM = D;
  template<typename I>
  inline static void resize(Array<T, D>& ref, I* n, int d)
  {
    if (d != DIM)
    {
      std::ostringstream err_msg;
      err_msg << "Array<T, " + DIM + "> cannot be resized. Requested dimension = " << d << std::endl;
      throw std::runtime_error(err_msg.str());
    }
    ref.resize(n);
  }
};
} // namespace qmcplusplus
#endif
