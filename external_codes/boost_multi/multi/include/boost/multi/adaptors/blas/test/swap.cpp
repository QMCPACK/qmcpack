// Copyright 2019-2024 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS swap"
// #include <boost/test/unit_test.hpp>

#include "../../blas.hpp"

#include <multi/array.hpp>

#include <cassert>
#include <complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(lapack_potrf, *boost::unit_test::tolerance(0.00001)) {
	{
		multi::array<double, 2> A = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};
		BOOST_REQUIRE( A[0][2] ==  3.0 );
		BOOST_REQUIRE( A[2][2] == 11.0 );

		multi::blas::swap(A[0], A[2]);  // blas swap
		BOOST_REQUIRE( A[0][2] == 11.0 );
		BOOST_REQUIRE( A[2][2] ==  3.0 );

		swap(A[0], A[2]);  // built-in swap
		BOOST_REQUIRE( A[0][2] ==  3.0 );
		BOOST_REQUIRE( A[2][2] == 11.0 );
	}
	{
		multi::array<double, 2> A = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};
		BOOST_REQUIRE( A[0][0] == 1.0 );
		BOOST_REQUIRE( A[0][3] == 4.0 );

		multi::blas::swap(rotated(A)[0], rotated(A)[3]);  // blas swap (deep)
		BOOST_REQUIRE( A[0][0] == 4.0 );
		BOOST_REQUIRE( A[0][3] == 1.0 );

		swap(rotated(A)[0], rotated(A)[3]);  // built-in swap (deep)
		BOOST_REQUIRE( A[0][0] == 1.0 );
		BOOST_REQUIRE( A[0][3] == 4.0 );
	}
	{
		using complex = std::complex<double>;
		complex const I{0, 1};
		multi::array<complex, 2> A = {
			{1.0 + 2. * I,  2.0,  3.0, 4.0 + 3.0 * I},
			{         5.0,  6.0,  7.0,           8.0},
			{         9.0, 10.0, 11.0,          12.0},
		};
		BOOST_REQUIRE( A[0][0] == 1.0 + 2.0*I );
		BOOST_REQUIRE( A[0][3] == 4.0 + 3.0*I );

		multi::blas::swap(rotated(A)[0], rotated(A)[3]);  // blas swap (deep)
		BOOST_REQUIRE( A[0][0] == 4.0 + 3.0*I );
		BOOST_REQUIRE( A[0][3] == 1.0 + 2.0*I );

		swap(rotated(A)[0], rotated(A)[3]);  // built-in swap (deep)
		BOOST_REQUIRE( A[0][0] == 1.0 + 2.0*I );
		BOOST_REQUIRE( A[0][3] == 4.0 + 3.0*I );
	}
	{
		multi::array<double, 2> A = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};
		BOOST_REQUIRE( A[0][2] ==  3.0 );
		BOOST_REQUIRE( A[2][2] == 11.0 );

		auto it = multi::blas::swap(begin(A[0]), end(A[0]) - 1, begin(A[2]));  // blas swap
		BOOST_REQUIRE( it == end(A[2]) - 1 );
		BOOST_REQUIRE( A[0][2] == 11.0 );
		BOOST_REQUIRE( A[2][2] ==  3.0 );

		using std::swap_ranges;
		swap_ranges(begin(A[0]), end(A[0]), begin(A[2]));  // built-in swap
		BOOST_REQUIRE( A[0][2] ==  3.0 );
		BOOST_REQUIRE( A[2][2] == 11.0 );
	}
}
