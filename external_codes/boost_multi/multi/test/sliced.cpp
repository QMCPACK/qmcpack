// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2021-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi slice"  // test tile NOLINT(cppcoreguidelines-macro-usage)
#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include<numeric>  // std::iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_sliced_empty) {
	multi::array<double, 2> arr({0, 0}, 99.0);
	BOOST_REQUIRE( arr.sliced(0, 0).is_empty() );
	BOOST_REQUIRE( arr.sliced(1, 1).is_empty() );
}

BOOST_AUTO_TEST_CASE(multi_array_sliced) {
	multi::array<double, 4> arr({10, 20, 30, 40}, 99.0);
	std::iota(arr.elements().begin(), arr.elements().end(), 0.);

	static_assert( decltype( arr.sliced(0, 5) )::rank::value == decltype(arr)::rank::value);  // NOLINT(misc-redundant-expression)
	static_assert( decltype( arr.sliced(0, 5) )::rank{} == decltype(arr)::rank{});  // NOLINT(misc-redundant-expression)
	static_assert( decltype( arr.sliced(0, 5) )::rank_v == decltype(arr)::rank_v);  // NOLINT(misc-redundant-expression)

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
	multi::array<double, 2> arr = {
		{ 1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0},
		{ 9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0},
	};
	BOOST_REQUIRE((
		arr.strided(2) == multi::array<double, 2>{
			{ 1.0,  2.0,  3.0,  4.0},
			{ 9.0, 10.0, 11.0, 12.0},
		}
	));
}

BOOST_AUTO_TEST_CASE(multi_array_take) {
	multi::array<double, 2> arr = {
		{ 1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0},
		{ 9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0},
	};
	BOOST_REQUIRE((
		arr.take(2) == multi::array<double, 2>{
			{ 1.0,  2.0,  3.0,  4.0},
			{ 5.0,  6.0,  7.0,  8.0},
		}
	));
}

BOOST_AUTO_TEST_CASE(drop) {
	multi::array<double, 2> arr = {
		{ 1.,  2.,  3.,  4.},
		{ 5.,  6.,  7.,  8.},
		{ 9., 10., 11., 12.},
		{13., 14., 15., 16.},
	};
	BOOST_REQUIRE((
		arr.drop(2) == multi::array<double, 2>{
			{ 9., 10., 11., 12.},
			{13., 14., 15., 16.},
		}
	));
}
