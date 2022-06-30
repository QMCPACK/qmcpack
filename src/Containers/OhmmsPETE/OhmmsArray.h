//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_PETE_ARRAY_H
#define OHMMS_PETE_ARRAY_H

/** A D-dimensional Array class based on PETE
 *
 */
#include <array>
#include <type_traits>
#include "OhmmsVector.h"

template<class T, unsigned D, typename ALLOC = std::allocator<T>>
class Array
{
public:
  using Type_t      = T;
  using Container_t = qmcplusplus::Vector<T, ALLOC>;
  using This_t      = Array<T, D, ALLOC>;

  //default constructor
  Array() { Length.fill(0); }

  //copy constructor
  Array(const Array& rhs)
  {
    resize(rhs.shape());
    std::copy(rhs.begin(), rhs.end(), X.begin());
  }

  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  Array(const std::array<SIZET, D>& dims)
  {
    resize(dims);
  }

  /// Provides specialized constructors with signature (size_1, ... , size_D) for  array dimension D
  template<typename... Args>
  Array(Args... sizes)
  {
    resize(sizes...);
  }

  inline unsigned dim() const { return D; }
  inline const std::array<size_t, D>& shape() const { return Length; }
  inline size_t size() const { return X.size(); }
  inline size_t size(int i) const { return Length[i]; }

  Container_t& storage() { return X; }

  /// Resize the container. For performance consideration, previous data may or may not get kept.
  /// Please avoid relying on previous data after resizing.
  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  void resize(const std::array<SIZET, D>& dims)
  {
    for (int i = 0; i < dims.size(); i++)
      Length[i] = dims[i];
    X.resize(full_size(Length));
  }

  /// Provides specialized resize(size_1, ... , size_D) functions for the array D
  template<typename... Args>
  void resize(Args... sizes)
  {
    static_assert(sizeof...(Args) == D, "resize arguments must match dimensionality of Array");
    resize({static_cast<std::size_t>(std::forward<Args>(sizes))...});
  }

  inline typename Container_t::iterator begin() { return X.begin(); }
  inline typename Container_t::iterator end() { return X.end(); }
  inline typename Container_t::const_iterator begin() const { return X.begin(); }
  inline typename Container_t::const_iterator end() const { return X.end(); }

  ///@{
  /// access the container data pointer
  inline Type_t* data() { return X.data(); }
  inline const Type_t* data() const { return X.data(); }
  template<typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  inline Type_t* device_data()
  {
    return X.device_data();
  }
  template<typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  inline const Type_t* device_data() const
  {
    return X.device_data();
  }
  ///@}

  ///@{
  /// access the data pointer at {index_1, ..., index_D}
  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  Type_t* data_at(const std::array<SIZET, D>& indices)
  {
    return X.data() + compute_offset(indices);
  }
  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  const Type_t* data_at(const std::array<SIZET, D>& indices) const
  {
    return X.data() + compute_offset(indices);
  }
  template<typename SIZET     = size_t,
           typename           = std::is_integral<SIZET>,
           typename Allocator = ALLOC,
           typename           = qmcplusplus::IsDualSpace<Allocator>>
  Type_t* device_data_at(const std::array<SIZET, D>& indices)
  {
    return X.device_data() + compute_offset(indices);
  }
  template<typename SIZET     = size_t,
           typename           = std::is_integral<SIZET>,
           typename Allocator = ALLOC,
           typename           = qmcplusplus::IsDualSpace<Allocator>>
  const Type_t* device_data_at(const std::array<SIZET, D>& indices) const
  {
    return X.device_data() + compute_offset(indices);
  }

  template<typename... Args>
  Type_t* data_at(Args... indices)
  {
    static_assert(sizeof...(Args) == D, "data arguments must match dimensionality of Array");
    return data_at({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  template<typename... Args>
  const Type_t* data_at(Args... indices) const
  {
    static_assert(sizeof...(Args) == D, "data arguments must match dimensionality of Array");
    return data_at({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  template<typename... Args, typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  Type_t* device_data_at(Args... indices)
  {
    static_assert(sizeof...(Args) == D, "device_data arguments must match dimensionality of Array");
    return device_data_at({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  template<typename... Args, typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  const Type_t* device_data_at(Args... indices) const
  {
    static_assert(sizeof...(Args) == D, "device_data arguments must match dimensionality of Array");
    return device_data_at({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  ///@}

  inline const Type_t* first_address() const { return X.data(); }

  inline const Type_t* last_address() const { return X.data() + X.size(); }

  inline Type_t* first_address() { return X.data(); }

  inline Type_t* last_address() { return X.data() + X.size(); }

  This_t& operator=(const T& rhs)
  {
    std::fill(X.begin(), X.end(), rhs);
    return *this;
  }

  This_t& operator=(const Array& rhs)
  {
    if (&rhs != this)
    {
      resize(rhs.shape());
      std::copy(rhs.begin(), rhs.end(), X.begin());
    }
    return *this;
  }

  template<typename TT, typename ALLOC2>
  This_t& operator=(const Array<TT, D, ALLOC2>& rhs)
  {
    resize(rhs.shape());
    std::copy(rhs.begin(), rhs.end(), X.begin());
    return *this;
  }

  ///@{
  /// access the element at {index_1, ..., index_D}
  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  Type_t& operator()(const std::array<SIZET, D>& indices)
  {
    return X[compute_offset(indices)];
  }
  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  const Type_t& operator()(const std::array<SIZET, D>& indices) const
  {
    return X[compute_offset(indices)];
  }
  template<typename... Args>
  Type_t& operator()(Args... indices)
  {
    static_assert(sizeof...(Args) == D, "operator() arguments must match dimensionality of Array");
    return operator()({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  template<typename... Args>
  const Type_t& operator()(Args... indices) const
  {
    static_assert(sizeof...(Args) == D, "operator() arguments must match dimensionality of Array");
    return operator()({static_cast<std::size_t>(std::forward<Args>(indices))...});
  }
  ///@}

  inline Type_t sum() const
  {
    Type_t s = 0;
    for (int i = 0; i < X.size(); ++i)
      s += X[i];
    return s;
  }

  // Abstract Dual Space Transfers
  template<typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  void updateTo()
  {
    X.updateTo();
  }
  template<typename Allocator = ALLOC, typename = qmcplusplus::IsDualSpace<Allocator>>
  void updateFrom()
  {
    X.updateFrom();
  }

private:
  std::array<size_t, D> Length;
  Container_t X;
  /// compute the full size of dims
  size_t full_size(const std::array<size_t, D>& dims) const
  {
    size_t total = dims[0];
    for (int i = 1; i < dims.size(); i++)
      total *= dims[i];
    return total;
  }

  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  SIZET compute_offset(const std::array<SIZET, D>& indices) const
  {
    SIZET offset = indices[0];
    for (int i = 1; i < indices.size(); i++)
      offset = offset * Length[i] + indices[i];
    return offset;
  }
};

template<class T, unsigned D, class Alloc>
bool operator==(const Array<T, D, Alloc>& lhs, const Array<T, D, Alloc>& rhs)
{
  static_assert(qmcplusplus::qmc_allocator_traits<Alloc>::is_host_accessible,
                "operator== requires host accessible Vector.");
  if (lhs.size() == rhs.size())
  {
    for (int i = 0; i < rhs.size(); i++)
      if (lhs(i) != rhs(i))
        return false;
    return true;
  }
  else
    return false;
}

template<class T, unsigned D, class Alloc>
bool operator!=(const Array<T, D, Alloc>& lhs, const Array<T, D, Alloc>& rhs)
{
  static_assert(qmcplusplus::qmc_allocator_traits<Alloc>::is_host_accessible,
                "operator== requires host accessible Vector.");
  return !(lhs == rhs);
}
#endif //OHMMS_PETE_ARRAY_H
