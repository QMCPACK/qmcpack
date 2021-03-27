//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CHECKMATRIX_HPP
#define QMCPLUSPLUS_CHECKMATRIX_HPP

#include "catch.hpp"

#include <complex>
#include <type_traits>
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus {

template<typename T>
struct implIsComplex : public std::false_type {};
template<typename T>
struct implIsComplex<std::complex<T>> : public std::true_type {};

template<typename T>
using IsComplex = std::enable_if_t< implIsComplex<T>::value, bool>;
template<typename T>
using IsReal = std::enable_if_t< std::is_floating_point<T>::value, bool>;

  
template<typename T, IsComplex<T> = true >
bool approxEquality(T val_a, T val_b) {
  //  static_assert(!IsComplex<T>::value, "ComplexApprox is for Complex!");
  return val_a == ComplexApprox(val_b);
}

template<typename T, IsReal<T>  = true >
bool approxEquality(T val_a, T val_b) {
  //  static_assert(!IsComplex<T>::value, "Approx is for Reals!");
  return val_a == Approx(val_b);
}

template<typename T1, typename ALLOC1, typename T2, typename ALLOC2>
void checkMatrix(const Matrix<T1, ALLOC1>& a, const Matrix<T2, ALLOC2>& b, const std::string & desc = "", int line = 0)
{
  REQUIRE(a.rows() >= b.rows());
  REQUIRE(a.cols() >= b.cols());
  auto matrixElementError = [line](int i, int j, auto& a, auto& b, const std::string& desc) -> std::string {
                      	          std::stringstream error_msg;
				  error_msg << "In " << desc << ":" << line <<  "\nbad element at " << i << ":" << j
					    <<"  " << a(i,j) << " != " << b(i,j) << '\n';
				  return error_msg.str();
			    };
  for (int i = 0; i < b.rows(); i++)
    for (int j = 0; j < b.cols(); j++)
    {
      bool approx_equality = approxEquality<T1>(a(i, j), b(i, j));
      CHECKED_ELSE( approx_equality ) {
	FAIL( matrixElementError(i,j,a,b,desc) );
	      }
    }
}

  extern template bool approxEquality<double>(double val_a, double val_b);
  extern template bool approxEquality<std::complex<double>>(std::complex<double> val_a, std::complex<double> val_b);
  
}

#endif
