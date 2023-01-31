#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS swap"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas.hpp"

#include "../../../array.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(lapack_potrf, *boost::unit_test::tolerance(0.00001) ){
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		BOOST_REQUIRE( A[0][2] == 3. );
		BOOST_REQUIRE( A[2][2] == 11. );

		multi::blas::swap(A[0], A[2]); // blas swap
		BOOST_REQUIRE( A[0][2] == 11. );
		BOOST_REQUIRE( A[2][2] == 3. );

		swap(A[0], A[2]); // built-in swap
		BOOST_REQUIRE( A[0][2] == 3. );
		BOOST_REQUIRE( A[2][2] == 11. );
	}
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		BOOST_REQUIRE( A[0][0] == 1. );
		BOOST_REQUIRE( A[0][3] == 4. );

		multi::blas::swap(rotated(A)[0], rotated(A)[3]); // blas swap (deep)
		BOOST_REQUIRE( A[0][0] == 4. );
		BOOST_REQUIRE( A[0][3] == 1. );

		             swap(rotated(A)[0], rotated(A)[3]); // built-in swap (deep)
		BOOST_REQUIRE( A[0][0] == 1. );
		BOOST_REQUIRE( A[0][3] == 4. );
	}
	{
		using complex = std::complex<double>; complex const I{0, 1};
		multi::array<complex, 2> A = {
			{1.+ 2.*I,  2.,  3.,  4. + 3.*I},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		BOOST_REQUIRE( A[0][0] == 1.+ 2.*I );
		BOOST_REQUIRE( A[0][3] == 4. + 3.*I );
		multi::blas::swap(rotated(A)[0], rotated(A)[3]); // blas swap (deep)
		BOOST_REQUIRE( A[0][0] == 4. + 3.*I );
		BOOST_REQUIRE( A[0][3] == 1.+ 2.*I );
		             swap(rotated(A)[0], rotated(A)[3]); // built-in swap (deep)
		BOOST_REQUIRE( A[0][0] == 1.+ 2.*I );
		BOOST_REQUIRE( A[0][3] == 4. + 3.*I );
	}
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		BOOST_REQUIRE( A[0][2] == 3. );
		BOOST_REQUIRE( A[2][2] == 11. );

		auto it = multi::blas::swap(begin(A[0]), end(A[0]) - 1, begin(A[2])); // blas swap
		BOOST_REQUIRE( it == end(A[2]) - 1 );
		BOOST_REQUIRE( A[0][2] == 11. );
		BOOST_REQUIRE( A[2][2] == 3. );
		using std::swap_ranges;
		      swap_ranges(begin(A[0]), end(A[0]), begin(A[2])); // built-in swap
		BOOST_REQUIRE( A[0][2] == 3. );
		BOOST_REQUIRE( A[2][2] == 11. );
	}
}

