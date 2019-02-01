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

using std::string;
using std::complex;

namespace qmcplusplus
{

template<typename T>
void myREQUIRE(T const& a, T const& b)
{
  REQUIRE(a == Approx(b));
}

template<typename T>
void myREQUIRE(std::complex<T> const& a, std::complex<T> const& b)
{
  REQUIRE(a.real() == Approx(b.real()));
  REQUIRE(a.imag() == Approx(b.imag()));
}

template<class M1,
         class M2,
         typename = typename std::enable_if<(M1::dimensionality == 1)>::type,
         typename = typename std::enable_if<(M2::dimensionality == 1)>::type
         >
void verify_approx(M1 const& A, M2 const& B)
{
  REQUIRE(A.shape()[0] == B.shape()[0]);
  for(int i=0; i<A.shape()[0]; i++)
      myREQUIRE(A[i],B[i]);
}

template<class M1,
         class M2,
         typename = typename std::enable_if<(M1::dimensionality > 1)>::type,
         typename = typename std::enable_if<(M2::dimensionality > 1)>::type,
         typename = void
         >
void verify_approx(M1 const& A, M2 const& B)
{
  REQUIRE(A.shape()[0] == B.shape()[0]);
  for(int i=0; i<A.shape()[0]; i++)
    verify_approx(A[i],B[i]);
}

}

#endif
