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

#include <stdexcept>
#include <vector>

namespace qmcplusplus
{
template<typename CT>
struct container_traits
{
  /// the data type of elements
  using element_type = typename CT::value_type;

  /** resize container
   * @param n the size of all the dimensions
   * @param d the number of dimensions
   */
  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    throw std::runtime_error("Unknown container, resizing is not available!");
  }

  /// get the current linear storage size of a container
  inline static size_t getSize(const CT& ref) { return ref.size(); }

  /// get the linear storage pointer of a container
  inline static auto getElementPtr(CT& ref) { return ref.data(); }
};

// template specialization for std::vector
template<typename T, class ALLOC>
struct container_traits<std::vector<T, ALLOC>>
{
  using element_type = T;
  using CT           = std::vector<T, ALLOC>;

  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    size_t nt = d > 0 ? 1 : 0;
    for (int i = 0; i < d; ++i)
      nt *= n[i];
    ref.resize(nt);
  }

  inline static size_t getSize(const CT& ref) { return ref.size(); }

  inline static auto getElementPtr(CT& ref) { return ref.data(); }
};

} // namespace qmcplusplus
#endif
