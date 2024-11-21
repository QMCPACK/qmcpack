// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <complex>

// Suppress warnings from boost.test
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(comparison_complex) {
	using complex = std::complex<double>;
	{
		multi::array<double, 1>  arr  = {1.0, 2.0, 3.0};
		multi::array<complex, 1> arr2 = {
			{1.0, 0.0},
			{2.0, 0.0},
			{3.0, 0.0},
		};

		BOOST_REQUIRE( arr[1] == arr2[1] );
		BOOST_REQUIRE( arr == arr2 );
		BOOST_REQUIRE( ! (arr != arr2) );
		BOOST_REQUIRE( arr2 == arr );
		BOOST_REQUIRE( ! (arr2 != arr) );
	}
	{
		multi::array<double, 2> const arr = {
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
		};

		multi::array<complex, 2> const arr2 = {
			{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}},
			{{4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}},
		};

		BOOST_REQUIRE( arr[1][1] == arr2[1][1] );
		BOOST_REQUIRE( arr == arr2 );
		BOOST_REQUIRE( ! (arr != arr2) );
		BOOST_REQUIRE( arr2 == arr );
		BOOST_REQUIRE( ! (arr2 != arr) );
		BOOST_REQUIRE( std::equal(arr[1].begin(), arr[1].end(), begin(arr2[1]), end(arr2[1])) );
	}
}

BOOST_AUTO_TEST_CASE(multi_comparisons_swap) {
	multi::array<double, 3> arr = {
		{ {1.2, 1.1},  {2.4, 1.0}},
		{{11.2, 3.0}, {34.4, 4.0}},
		{ {1.2, 1.1},  {2.4, 1.0}},
	};
	BOOST_REQUIRE( arr[0] < arr[1] );

	swap(arr[0], arr[1]);
	BOOST_REQUIRE( arr[1] < arr[0] );

	swap(arr[0], arr[1]);
	BOOST_REQUIRE( arr[0] < arr[1] );
}

BOOST_AUTO_TEST_CASE(comparisons_equality) {
	multi::array<double, 3> arr = {
		{ {1.2, 1.1},  {2.4, 1.0}},
		{{11.2, 3.0}, {34.4, 4.0}},
		{ {1.2, 1.1},  {2.4, 1.0}},
	};

	multi::array_ref<double, 3>  ref(arr.data_elements(), extensions(arr));
	multi::array_cref<double, 3> cref(data_elements(arr), extensions(arr));

	BOOST_REQUIRE( arr ==  arr );
	BOOST_REQUIRE( ! (arr !=  arr) );
	BOOST_REQUIRE( ref ==  arr );
	BOOST_REQUIRE( ! (ref !=  arr) );
	BOOST_REQUIRE( ref == cref );
	BOOST_REQUIRE( ! (ref != cref) );

	BOOST_REQUIRE( arr[0] ==  arr[2] );
	BOOST_REQUIRE( ref[0] ==  arr[2] );
	BOOST_REQUIRE( ref[0] == cref[2] );

	BOOST_REQUIRE( ! ( arr[0] != arr[2]) );
	BOOST_REQUIRE( ! ( ref[0] != ref[2]) );

	BOOST_REQUIRE( ! ( arr[0] != arr[2]) );
	BOOST_REQUIRE( ! ( ref[0] != ref[2]) );
}

BOOST_AUTO_TEST_CASE(comparisons_ordering) {
	multi::array<double, 3> arr = {
		{ {12, 11},  {24, 10}},
		{{112, 30}, {344, 40}},
		{ {12, 11},  {24, 10}},
	};

	multi::array_ref<double, 3> ref(arr.data_elements(), extensions(arr));

	multi::array_cref<double, 3> cref(data_elements(arr), extensions(arr));

	BOOST_REQUIRE(  arr[0]    <=  arr[1] );
	BOOST_REQUIRE(  ref[0]    <=  arr[1] );
	BOOST_REQUIRE( cref[0]    <= cref[1] );

	BOOST_REQUIRE(  arr[0][0] <= arr[0][1] );
	BOOST_REQUIRE(  ref[0][0] <= arr[0][1] );

	BOOST_REQUIRE(  arr[1][0][0] == 112 );
	BOOST_REQUIRE(  ref[1][0][0] == 112 );
	BOOST_REQUIRE( cref[1][0][0] == 112 );

	BOOST_REQUIRE(  arr[0][0][0] == 12 );
	BOOST_REQUIRE(  ref[0][0][0] == 12 );
	BOOST_REQUIRE( cref[0][0][0] == 12 );

	swap(ref[0], ref[1]);

	BOOST_REQUIRE(  begin(arr) <  end(arr) );
	BOOST_REQUIRE( cbegin(arr) < cend(arr) );

	BOOST_REQUIRE(  end(arr) -  begin(arr) == size(arr) );
}
