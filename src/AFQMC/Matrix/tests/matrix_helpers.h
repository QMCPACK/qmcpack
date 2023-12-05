//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SPARSE_MATRIX_HELPER_H
#define AFQMC_SPARSE_MATRIX_HELPER_H

#include <stdio.h>
#include <string>
#include <complex>
#include <type_traits>

#include "multi/array.hpp"

using std::complex;
using std::string;

namespace qmcplusplus
{
template<typename T>
void myCHECK(T const& a, T const& b)
{
  CHECK(a == Approx(b));
}

template<typename T>
void myCHECK(std::complex<T> const& a, std::complex<T> const& b)
{
  CHECK(a.real() == Approx(b.real()));
  CHECK(a.imag() == Approx(b.imag()));
}

template<class M1,
         class M2,
         typename = typename std::enable_if<(M1::dimensionality == 1)>::type,
         typename = typename std::enable_if<(M2::dimensionality == 1)>::type>
void verify_approx(M1 const& A, M2 const& B)
{
  // casting in case operator[] returns a fancy reference
  using element1 = typename std::decay<M1>::type::element;
  using element2 = typename std::decay<M2>::type::element;
  REQUIRE(std::get<0>(A.sizes()) == std::get<0>(B.sizes()));
  for (int i = 0; i < std::get<0>(A.sizes()); i++)
    myCHECK(element1(A[i]), element2(B[i]));
}

template<class M1,
         class M2,
         typename = typename std::enable_if<(M1::dimensionality > 1)>::type,
         typename = typename std::enable_if<(M2::dimensionality > 1)>::type,
         typename = void>
void verify_approx(M1 const& A, M2 const& B)
{
  REQUIRE(A.size() == B.size());
  for (int i = 0; i < A.size(); i++)
    verify_approx(A[i], B[i]);
}

} // namespace qmcplusplus

#endif
