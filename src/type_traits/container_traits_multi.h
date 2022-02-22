//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONTAINER_TRAITS_MULTI_H
#define QMCPLUSPLUS_CONTAINER_TRAITS_MULTI_H

#include <multi/array.hpp>
#include <multi/array_ref.hpp>
#include "type_traits/container_traits.h"

namespace qmcplusplus
{
// template specialization for functions in container_traits
template<typename T, boost::multi::dimensionality_type D, class Alloc>
struct container_traits<boost::multi::array<T, D, Alloc>>
{
  using element_type = T;
  using CT           = boost::multi::array<T, D, Alloc>;

  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    if (d != D)
    {
      std::ostringstream err_msg;
      err_msg << "boost::multi::array<T, " << D << ", Alloc> cannot be resized. Requested dimension = " << d
              << std::endl;
      throw std::runtime_error(err_msg.str());
    }
    std::array<I, 2> shape;
    for (int i = 0; i < d; ++i)
      shape[i] = n[i];
    ref.reextent({static_cast<boost::multi::size_t>(shape[0]), static_cast<boost::multi::size_t>(shape[1])});
  }

  inline static size_t getSize(const CT& ref) { return ref.num_elements(); }

  inline static auto getElementPtr(CT& ref) { return std::addressof(*ref.origin()); }
};

template<typename T, boost::multi::dimensionality_type D>
struct container_traits<boost::multi::array_ref<T, D>>
{
  using element_type = T;
  using CT           = boost::multi::array_ref<T, D>;

  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    throw std::runtime_error("Can not resize container_proxy<boost::multi::array_ref<T,D>>!\n");
  }

  inline static size_t getSize(const CT& ref) { return ref.num_elements(); }

  inline static auto getElementPtr(CT& ref) { return std::addressof(*ref.origin()); }
};

} // namespace qmcplusplus

#endif
