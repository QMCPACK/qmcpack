//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_RAW_POINTERS_DETAIL_HPP
#define AFQMC_RAW_POINTERS_DETAIL_HPP

#include <type_traits>
#include <complex>
#include "AFQMC/Utilities/type_conversion.hpp"

namespace qmcplusplus
{
namespace afqmc
{
template<class T, typename = typename std::enable_if_t<std::is_fundamental<T>::value>>
inline static T* to_address(T* p)
{
  return p;
}

template<class T>
inline static std::complex<T>* to_address(std::complex<T>* p)
{
  return p;
}

template<class T>
inline static std::complex<T> const* to_address(std::complex<T> const* p)
{
  return p;
}

//  template<class Q, class T,
//           typename = typename std::enable_if_t<std::is_fundamental<Q>::value>,
//           typename = typename std::enable_if_t<std::is_fundamental<T>::value>>
//  inline static Q* pointer_cast(T* p) { return reinterpret_cast<Q*>(p); }

//  template<class Q, class T>
//  inline static Q* pointer_cast(std::complex<T>* p) { return reinterpret_cast<Q*>(p); }

template<class Q, class T>
inline static Q* pointer_cast(T* p)
{
  return reinterpret_cast<Q*>(p);
}


/************* copy_n_cast ****************/
template<class T, class Q, class Size>
Q* copy_n_cast(T const* A, Size n, Q* B)
{
  for (Size i = 0; i < n; i++, ++A, ++B)
    *B = static_cast<Q>(*A);
  return B;
}

/************* inplace_cast ****************/
template<class T, class Q, class Size>
void inplace_cast(T* A, Size n)
{
  Q* B(reinterpret_cast<Q*>(A));
  if (sizeof(T) >= sizeof(Q))
  {
    for (Size i = 0; i < n; i++, ++A, ++B)
      *B = static_cast<Q>(*A);
  }
  else if (sizeof(T) < sizeof(Q))
  {
    assert(sizeof(T) * 2 <= sizeof(Q));
    A += (n - 1);
    B += (n - 1);
    for (; n > 0; n--, --A, --B)
      *B = static_cast<Q>(*A);
  }
}

/************* fill2D ****************/
template<typename T, typename Size>
void fill2D(Size N, Size M, T* y, Size lda, T const a)
{
  for (Size ip = 0; ip < N; ++ip)
    for (Size jp = 0; jp < M; ++jp)
    {
      y[ip * lda + jp] = a;
    }
}

/************* print ****************/
template<typename T>
void print(std::string str, T const* p, int n)
{
  std::cout << str << " ";
  for (int i = 0; i < n; i++)
    std::cout << *(p + i) << " ";
  std::cout << std::endl;
}

} // namespace afqmc
} // namespace qmcplusplus

#endif
