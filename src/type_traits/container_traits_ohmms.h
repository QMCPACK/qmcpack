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


#ifndef QMCPLUSPLUS_CONTAINER_TRAITS_OHMMS_H
#define QMCPLUSPLUS_CONTAINER_TRAITS_OHMMS_H

#include <sstream>
#include "type_traits/container_traits.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{
// template specialization for Ohmms Vector/Matrix/Array
template<typename T, class ALLOC>
struct container_traits<Vector<T, ALLOC>>
{
  using element_type = T;
  using CT           = Vector<T, ALLOC>;

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

template<typename T, class ALLOC>
struct container_traits<Matrix<T, ALLOC>>
{
  using element_type = T;
  using CT           = Matrix<T, ALLOC>;

  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    if (d != 2)
    {
      std::ostringstream err_msg;
      err_msg << "Matrix cannot be resized. Requested dimension = " << d << std::endl;
      throw std::runtime_error(err_msg.str());
    }
    ref.resize(n[0], n[1]);
  }

  inline static size_t getSize(const CT& ref) { return ref.size(); }

  inline static auto getElementPtr(CT& ref) { return ref.data(); }
};

template<typename T, unsigned D>
struct container_traits<Array<T, D>>
{
  using element_type = T;
  using CT           = Array<T, D>;

  template<typename I>
  inline static void resize(CT& ref, I* n, int d)
  {
    if (d != D)
    {
      std::ostringstream err_msg;
      err_msg << "Array<T, " << D << "> cannot be resized. Requested dimension = " << d << std::endl;
      throw std::runtime_error(err_msg.str());
    }
    ref.resize(n);
  }

  inline static size_t getSize(const CT& ref) { return ref.size(); }

  inline static auto getElementPtr(CT& ref) { return ref.data(); }
};

} // namespace qmcplusplus
#endif
