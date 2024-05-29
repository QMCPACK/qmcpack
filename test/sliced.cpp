// Copyright 2021-2023 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <numeric>  // std::iota

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
// #  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
// #  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
#  pragma clang diagnostic pop
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_sliced_empty) {
	multi::array<double, 2> const arr({0, 0}, 99.0);
	BOOST_REQUIRE( arr.sliced(0, 0).is_empty() );
	// BOOST_REQUIRE( arr.sliced(1, 1).is_empty() );  // this results in offsetting nullptr
}

BOOST_AUTO_TEST_CASE(multi_array_sliced) {
	multi::array<int, 4> arr({10, 20, 30, 40}, 99);
	std::iota(arr.elements().begin(), arr.elements().end(), 0);

	static_assert( decltype( arr.sliced(0, 5) )::rank::value == 4);
	static_assert( decltype( arr.sliced(0, 5) )::rank{} == 4);
	static_assert( decltype( arr.sliced(0, 5) )::rank_v == 4);

	BOOST_REQUIRE(  arr.sliced( 0, 5)[1][2][3][4] ==  arr[1][2][3][4] );
	BOOST_REQUIRE( &arr.sliced( 0, 5)[1][2][3][4] == &arr[1][2][3][4] );

	BOOST_REQUIRE(  arr.sliced( 0, 5)[1] ==  arr[1] );
	BOOST_REQUIRE( &arr.sliced( 0, 5)[1] == &arr[1] );

	BOOST_REQUIRE(  arr.sliced( 0,  0).empty() );
	BOOST_REQUIRE(  arr.sliced( 1,  1).empty() );
	BOOST_REQUIRE(  arr.sliced( 0, 10).size() == 10 );

	BOOST_REQUIRE(  arr[1].sliced(0, 5)[2][3][4] ==  arr[1][2][3][4] );
	BOOST_REQUIRE( &arr[1].sliced(0, 5)[2][3][4] == &arr[1][2][3][4] );

	BOOST_REQUIRE(  arr[1].sliced(0, 5)[2] ==  arr[1][2] );
	BOOST_REQUIRE( &arr[1].sliced(0, 5)[2] == &arr[1][2] );

	BOOST_REQUIRE( arr[1].sliced(0,  0).is_empty() );
	BOOST_REQUIRE( arr[1].sliced(1,  1).is_empty() );
	BOOST_REQUIRE( arr[1].sliced(0, 20).size() == 20 );

	BOOST_REQUIRE(  (arr.rotated()).sliced(0, 5)[1][2][3][4] ==  (arr.rotated())[1][2][3][4] );
	BOOST_REQUIRE( &(arr.rotated()).sliced(0, 5)[1][2][3][4] == &(arr.rotated())[1][2][3][4] );
}

BOOST_AUTO_TEST_CASE(multi_array_stride) {
	multi::array<int, 2> arr = {
		{ 10,  20,  30,  40},
		{ 50,  60,  70,  80},
		{ 90, 100, 110, 120},
		{130, 140, 150, 160},
	};
	BOOST_REQUIRE((
		arr.strided(2) == multi::array<int, 2>{
			{ 10,  20,  30,  40},
			{ 90, 100, 110, 120},
		}
	));
}

BOOST_AUTO_TEST_CASE(multi_array_take) {
	multi::array<int, 2> arr = {
		{ 10,  20,  30,  40},
		{ 50,  60,  70,  80},
		{ 90, 100, 110, 120},
		{130, 140, 150, 160},
	};
	BOOST_REQUIRE((
		arr.taked(2) == multi::array<int, 2>{
			{ 10,  20,  30,  40},
			{ 50,  60,  70,  80},
		}
	));
}

BOOST_AUTO_TEST_CASE(drop) {
	multi::array<double, 2> arr = {
		{ 10,  20,  30,  40},
		{ 50,  60,  70,  80},
		{ 90, 100, 110, 120},
		{130, 140, 150, 160},
	};
	BOOST_REQUIRE((
		arr.dropped(2) == multi::array<double, 2>{
			{ 90, 100, 110, 120},
			{130, 140, 150, 160},
		}
	));
}
