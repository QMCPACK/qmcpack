// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/test/unit_test.hpp>

#include <boost/multi/adaptors/blas/scal.hpp>

#include <boost/multi/array.hpp>

namespace multi = boost::multi;
namespace blas  = multi::blas;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_n) {
	multi::array<double, 2> arr = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( (arr[0][2] == 3.0) && (arr[2][2] == 11.0) );

	blas::scal_n(2.0, arr[2].begin(), arr[2].size());
	BOOST_REQUIRE( arr[0][2] == 3. && arr[2][2] == 11.0*2.0 );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_it) {
	multi::array<double, 2> arr = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] == 11.0 );

	blas::scal(2.0, arr[2].begin(), arr[2].end());
	BOOST_REQUIRE( arr[0][2] == 3.0 );
	BOOST_REQUIRE(arr[2][2] == 11.0*2.0 );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_real) {
	multi::array<double, 2> arr = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] == 11.0 );

	BOOST_REQUIRE(  blas::scal(1.0, arr[2]) ==  arr[2] );
	BOOST_REQUIRE( &blas::scal(1.0, arr[2]) == &arr[2] );
	BOOST_REQUIRE( +blas::scal(1.0, arr[2]) ==  arr[2] );

	blas::scal(2.0, arr[2]);
	BOOST_REQUIRE( arr[0][2] == 3.0 && arr[2][2] == 11.0*2.0 );

	BOOST_REQUIRE( &blas::scal(1.0, arr[2]) == &arr[2] );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_real_2D) {
	multi::array<double, 2> arr = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] == 11.0 );

	blas::scal(2.0, arr.elements());

	BOOST_REQUIRE( arr[0][2] ==  6.0 );
	BOOST_REQUIRE( arr[2][2] == 22.0 );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_complex_2D) {
	auto const I = std::complex<double>(0.0, 1.0);  // NOLINT(readability-identifier-length) blas conventional name
	multi::array<std::complex<double>, 2> arr = {
		{1.0 + 0.0*I,  2.0 + 0.0*I,  3.0 + 0.0*I,  4.0 + 0.0*I},
		{5.0 + 0.0*I,  6.0 + 0.0*I,  7.0 + 0.0*I,  8.0 + 0.0*I},
		{9.0 + 0.0*I, 10.0 + 0.0*I, 11.0 + 0.0*I, 12.0 + 0.0*I},
	};
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] == 11.0 );

	blas::scal(2.0, arr.elements());

	BOOST_REQUIRE( arr[0][2] ==  6.0 );
	BOOST_REQUIRE( arr[2][2] == 22.0 );
}
