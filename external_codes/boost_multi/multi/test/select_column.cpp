// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi select range"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_range_section_1D) {
	multi::array<double, 1> A = {00., 01., 02.}; (void)A;
	BOOST_REQUIRE( A == A(multi::ALL) );
	BOOST_REQUIRE( size(A( 1 <= multi::ALL )) == 2 );
	BOOST_REQUIRE( A( 1 <= multi::ALL )[0] == 1. );
	BOOST_REQUIRE( size(A( multi::ALL < 2 )) == 2 );
	BOOST_REQUIRE( A( multi::ALL < 2 )[1] == 1. );
}

BOOST_AUTO_TEST_CASE(multi_array_range_section_part1) {
	multi::array<double, 2> A = {
		{00., 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.},
		{30., 31., 32.},
	};

	using multi::_;
	using multi::U;

	BOOST_REQUIRE( size( A(      multi::ALL     , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(      multi::_       , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(     *multi::_       , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(      multi::U       , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(      multi::ALL     , 2) ) == 4 );
	BOOST_REQUIRE( size( A(      multi::ALL < 2 , 2) ) == 2 );
	BOOST_REQUIRE( size( A( 1 <= multi::ALL     , 2) ) == 3 );
	BOOST_REQUIRE( size( A( 1 <= multi::ALL < 3 , 2) ) == 2 );

	BOOST_REQUIRE( size( A(      multi::_       , 2) ) == 4 );
	BOOST_REQUIRE( size( A(      multi::_   < 2 , 2) ) == 2 );
	BOOST_REQUIRE( size( A( 1 <= multi::_       , 2) ) == 3 );
	BOOST_REQUIRE( size( A( 1 <= multi::_   < 3 , 2) ) == 2 );

	BOOST_REQUIRE( size( A(             _       , 2) ) == 4 );
	BOOST_REQUIRE( size( A(             _ < 2   , 2) ) == 2 );
	BOOST_REQUIRE( size( A( 1 <=        _       , 2) ) == 3 );
	BOOST_REQUIRE( size( A( 1 <=        _ < 3   , 2) ) == 2 );
}

BOOST_AUTO_TEST_CASE(multi_array_range_section_part2) {
	multi::array<double, 2> A = {
		{00., 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.},
		{30., 31., 32.},
	};

	BOOST_REQUIRE( size( A(A.extension(), 2) ) == size(A) );

	auto&& col2( A(A.extension(), 2) ); // select column #2
	// same as A(extesion(A), 2)
	// same as A(A.extension(0), 2);
	// same as rotated(A)[2];
//  BOOST_REQUIRE( col2.size(0) == size(A) );

	BOOST_REQUIRE( dimensionality(col2) == 1 );
	BOOST_REQUIRE( size(col2) == size(A) );
	BOOST_REQUIRE( col2.size() == size(A) );
	BOOST_REQUIRE( col2.stride() == 3 );
	BOOST_REQUIRE( col2[0] == 02. );
	BOOST_REQUIRE( col2[1] == 12. );
	BOOST_REQUIRE(( col2 == multi::array<double, 1>{02., 12., 22., 32.} ));
	BOOST_REQUIRE(( col2 == multi::array<double, 1>(rotated(A)[2]) ));
	BOOST_REQUIRE(( col2 == rotated(A)[2] ));
	BOOST_REQUIRE(( col2 == A(A.extension(), 2) ));
}

BOOST_AUTO_TEST_CASE(multi_array_range_section_syntax) {
	multi::array<double, 2> A = {
		{00., 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.},
		{30., 31., 32.},
	};

	using multi::_;
	BOOST_REQUIRE( size( A(       _       , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(      *_       , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(      (_)      , 2) ) == size(A) );

	using multi::U;
	BOOST_REQUIRE( size( A(       U       , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(       U       , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(       U       , 2) ) == size(A) );

	using multi::V;
	BOOST_REQUIRE( size( A(       V       , 2) ) == size(A) );
	BOOST_REQUIRE( size( A(       V       , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(       V       , 2) ) == size(A) );

//	using multi::A;
//	BOOST_REQUIRE( size( A(       A       , 2) ) == size(A) );
//	BOOST_REQUIRE( size( A(       A       , 2) ) == size(A) );

//	BOOST_REQUIRE( size( A(       A       , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(       _  < 2  , 2) ) == 2 );
	BOOST_REQUIRE( size( A(      *_  < 2  , 2) ) == 2 );
	BOOST_REQUIRE( size( A(       U  < 2  , 2) ) == 2 );

	BOOST_REQUIRE( size( A( 1 <=  _       , 2) ) == 3 );
	BOOST_REQUIRE( size( A( 1 <= *_       , 2) ) == 3 );
	BOOST_REQUIRE( size( A( 1 <=  U       , 2) ) == 3 );

	BOOST_REQUIRE( size( A( 1 <=  _ < 3   , 2) ) == 2 );
	BOOST_REQUIRE( size( A( 1 <= *_ < 3   , 2) ) == 2 );
	BOOST_REQUIRE( size( A( 1 <=  U < 3   , 2) ) == 2 );

	BOOST_REQUIRE( size( A(      *_ < 2   , 2) ) == 2 );
	BOOST_REQUIRE( size( A(       U < 2   , 2) ) == 2 );
}

