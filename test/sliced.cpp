// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2021-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi slice"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<numeric> // std::iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_sliced_empty) {
	multi::array<double, 2> A({0, 0}, 99.);
	BOOST_REQUIRE( A.sliced(0, 0).is_empty() );
	BOOST_REQUIRE( A.sliced(1, 1).is_empty() );
}

BOOST_AUTO_TEST_CASE(multi_array_sliced) {
	multi::array<double, 4> A({10, 20, 30, 40}, 99.);
	std::iota(A.elements().begin(), A.elements().end(), 0.);

	static_assert( decltype( A.sliced(0, 5) )::rank_v == decltype(A)::rank_v , "!"); //NOLINT(misc-redundant-expression)

	BOOST_REQUIRE(  A.sliced( 0, 5)[1][2][3][4] ==  A[1][2][3][4] );
	BOOST_REQUIRE( &A.sliced( 0, 5)[1][2][3][4] == &A[1][2][3][4] );

	BOOST_REQUIRE(  A.sliced( 0, 5)[1] ==  A[1] );
	BOOST_REQUIRE( &A.sliced( 0, 5)[1] == &A[1] );

	BOOST_REQUIRE(  A.sliced( 0,  0).empty() );
	BOOST_REQUIRE(  A.sliced( 1,  1).empty() );
	BOOST_REQUIRE(  A.sliced( 0, 10).size() == 10 );

	BOOST_REQUIRE(  A[1].sliced(0, 5)[2][3][4] ==  A[1][2][3][4] );
	BOOST_REQUIRE( &A[1].sliced(0, 5)[2][3][4] == &A[1][2][3][4] );

	BOOST_REQUIRE(  A[1].sliced(0, 5)[2] ==  A[1][2] );
	BOOST_REQUIRE( &A[1].sliced(0, 5)[2] == &A[1][2] );

	BOOST_REQUIRE( A[1].sliced(0,  0).is_empty() );
	BOOST_REQUIRE( A[1].sliced(1,  1).is_empty() );
	BOOST_REQUIRE( A[1].sliced(0, 20).size() == 20 );

	BOOST_REQUIRE(  (A.rotated()).sliced(0, 5)[1][2][3][4] ==  (A.rotated())[1][2][3][4] );
	BOOST_REQUIRE( &(A.rotated()).sliced(0, 5)[1][2][3][4] == &(A.rotated())[1][2][3][4] );
}

BOOST_AUTO_TEST_CASE(multi_array_stride) {
	multi::array<double, 2> A = {
		{ 1.,  2.,  3.,  4.},
		{ 5.,  6.,  7.,  8.},
		{ 9., 10., 11., 12.},
		{13., 14., 15., 16.},
	};
	BOOST_REQUIRE((
		A.strided(2) == multi::array<double, 2>{
			{ 1.,  2.,  3.,  4.},
			{ 9., 10., 11., 12.},
		}
	));
}
