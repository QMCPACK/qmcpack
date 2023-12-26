// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS nrm2"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/blas.hpp>
#include <multi/array.hpp>

#include <multi/adaptors/complex.hpp>

#include<complex>

namespace multi = boost::multi;

using complex = multi::complex<double>;
constexpr complex I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

BOOST_AUTO_TEST_CASE(multi_blas_nrm2){
	namespace blas = multi::blas;

	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) blas conventional name
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0}
	};
	BOOST_REQUIRE( blas::nrm2(A[1]) == std::sqrt(blas::dot(A[1], A[1])) );

	{
		multi::array<complex, 1> const x = {1.0 + 1.0*I, 3.0 + 2.0*I, 3.0 + 4.0*I};  // NOLINT(readability-identifier-length) blas conventional name
		BOOST_REQUIRE( blas::dot(x, x) == (1.0 + 1.0*I)*(1.0 + 1.0*I) + (3.0 + 2.0*I)*(3.0 + 2.0*I) + (3.0 + 4.0*I)*(3.0 + 4.0*I) );
		using std::sqrt;
		BOOST_REQUIRE( blas::nrm2(x) == sqrt(norm(1.0 + 1.0*I) + norm(3.0 + 2.0*I) + norm(3.0 + 4.0*I)) );
	}
}
