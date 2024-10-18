// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_reversed_3d) {
	multi::array<double, 3> arr({30, 40, 50});

	BOOST_TEST_REQUIRE( arr.reversed().size() == 50 );

	BOOST_REQUIRE( & arr.reversed()[3][5][7] == &arr[7][5][3] );
}

template<class Array>
auto flatted_last(Array&& arr) {
	return reversed(flatted(transposed(reversed(std::forward<Array>(arr)))));
}

template<class Array>
auto partitioned_last(Array&& arr, multi::size_type n) {
	return reversed(transposed(partitioned(reversed(std::forward<Array>(arr)), n)));
}

BOOST_AUTO_TEST_CASE(multi_reversed_4d) {
	multi::array<double, 4> arr({13, 5, 7, 11});

	BOOST_TEST_REQUIRE( arr.reversed().size() == 11 );

	BOOST_REQUIRE( &arr.reversed()[1][2][3][4] == &arr[4][3][2][1] );


	BOOST_REQUIRE( std::get<0>( arr.reversed().transposed().flatted().reversed().sizes() ) == 13 );
	BOOST_REQUIRE( std::get<1>( arr.reversed().transposed().flatted().reversed().sizes() ) ==  5 );
	BOOST_REQUIRE( std::get<2>( arr.reversed().transposed().flatted().reversed().sizes() ) == 77 );

	BOOST_REQUIRE(( sizes(arr.reversed().transposed().flatted().reversed()) == decltype(sizes(arr.reversed().transposed().flatted().reversed())){13, 5, 77} ));

	BOOST_REQUIRE( &arr.reversed().transposed().flatted().reversed()[1][2][ 5] == & arr[1][2][0][ 5] );
	BOOST_REQUIRE( &arr.reversed().transposed().flatted().reversed()[1][2][10] == & arr[1][2][0][10] );
	BOOST_REQUIRE( &arr.reversed().transposed().flatted().reversed()[1][2][11] == & arr[1][2][1][ 0] );
	BOOST_REQUIRE( &arr.reversed().transposed().flatted().reversed()[1][2][12] == & arr[1][2][1][ 1] );

	BOOST_REQUIRE( & flatted_last(arr)[1][2][12] == & arr[1][2][1][1] );
}

BOOST_AUTO_TEST_CASE(multi_reversed_4d_partition_last) {
	multi::array<double, 4> arr({11, 5, 7, 12});

	BOOST_REQUIRE( arr.reversed().size() == 12 );

	BOOST_REQUIRE( & arr.reversed()[1][2][3][4] == &arr[4][3][2][1] );

	BOOST_REQUIRE( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][0][1] == & arr[1][2][3][1] );
	BOOST_REQUIRE( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][1][0] == & arr[1][2][3][4] );
	BOOST_REQUIRE( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][1][1] == & arr[1][2][3][5] );

	BOOST_REQUIRE( & partitioned_last(arr, 3)[1][2][3][1][1] == & arr[1][2][3][5] );
}
