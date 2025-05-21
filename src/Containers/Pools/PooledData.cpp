//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "PooledData.h"

template<typename T>
void PooledData<T>::add(const std::complex<T>* first, const std::complex<T>* last)
{
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  size_type dn = 2 * (last - first);
  myData.insert(myData.end(), dn, T{});
  while (first != last)
  {
    myData[Current++] = (*first).real();
    myData[Current++] = (*first).imag();
    ++first;
  }
}

template<typename T>
void PooledData<T>::add(std::complex<T>* first, std::complex<T>* last)
{
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  size_type dn = 2 * (last - first);
  myData.insert(myData.end(), dn, T{});
  while (first != last)
  {
    myData[Current++] = (*first).real();
    myData[Current++] = (*first).imag();
    ++first;
  }
}

template<typename T>
void PooledData<T>::get(std::complex<T>* first, std::complex<T>* last)
{
  auto getptr = first;
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  while (getptr != last)
  {
    getptr->real(myData[Current++]);
    getptr->imag(myData[Current++]);
    ++getptr;
  }
}


template struct PooledData<float>;
template struct PooledData<double>;
