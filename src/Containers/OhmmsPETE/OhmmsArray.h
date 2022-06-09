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

  inline unsigned dim() const { return D; }
  inline const std::array<size_t, D>& shape() const { return Length; }
  inline size_t size() const { return X.size(); }
  inline size_t size(int i) const { return Length[i]; }

  Container_t& storage() { return X; }

  template<typename SIZET = size_t, typename = std::is_integral<SIZET>>
  void resize(const std::array<SIZET, D>& dims)
  {
    for (int i = 0; i < dims.size(); i++)
      Length[i] = dims[i];
    X.resize(full_size(Length));
  }

  void resize(size_t n) { resize(std::array<size_t, 1>{n}); }

  void resize(size_t m, size_t n) { resize({m, n}); }

  void resize(size_t l, size_t m, size_t n) { resize({l, m, n}); }

  inline typename Container_t::iterator begin() { return X.begin(); }
  inline typename Container_t::iterator end() { return X.end(); }
  inline typename Container_t::const_iterator begin() const { return X.begin(); }
  inline typename Container_t::const_iterator end() const { return X.end(); }

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

  inline const Type_t* first_address() const { return &(X[0]); }

  inline const Type_t* last_address() const { return &(X[0]) + X.size(); }

  inline Type_t* first_address() { return &(X[0]); }

  inline Type_t* last_address() { return &(X[0]) + X.size(); }

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

  // Get and Set Operations
  inline Type_t& operator()(size_t i) { return X[i]; }

  inline Type_t operator()(size_t i) const { return X[i]; }
  inline Type_t& operator()(size_t i, size_t j) { return X[j + Length[1] * i]; }
  inline Type_t operator()(size_t i, size_t j) const { return X[j + Length[1] * i]; }
  inline Type_t& operator()(size_t i, size_t j, size_t k) { return X[k + Length[2] * (j + Length[1] * i)]; }
  inline Type_t operator()(size_t i, size_t j, size_t k) const { return X[k + Length[2] * (j + Length[1] * i)]; }
  inline Type_t& operator()(size_t i, size_t j, size_t k, size_t l)
  {
    return X[l + Length[3] * (k + Length[2] * (j + Length[1] * i))];
  }
  inline Type_t operator()(size_t i, size_t j, size_t k, size_t l) const
  {
    return X[l + Length[3] * (k + Length[2] * (j + Length[1] * i))];
  }

  inline Type_t sum() const
  {
    Type_t s = 0;
    for (int i = 0; i < X.size(); ++i)
      s += X[i];
    return s;
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
};

template<class T, unsigned D, class Alloc>
bool operator==(const Array<T, D, Alloc>& lhs, const Array<T, D, Alloc>& rhs)
{
  static_assert(qmcplusplus::qmc_allocator_traits<Alloc>::is_host_accessible, "operator== requires host accessible Vector.");
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
  static_assert(qmcplusplus::qmc_allocator_traits<Alloc>::is_host_accessible, "operator== requires host accessible Vector.");
  return !(lhs == rhs);
}
#endif //OHMMS_PETE_ARRAY_H
