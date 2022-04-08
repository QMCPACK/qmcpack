#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi slice"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<numeric> // iota

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

	BOOST_REQUIRE(  A.sliced(0, 5)[1][2][3][4] ==  A[1][2][3][4] );
	BOOST_REQUIRE( &A.sliced(0, 5)[1][2][3][4] == &A[1][2][3][4] );

	BOOST_REQUIRE(  A.sliced(0, 5)[1] ==  A[1] );
	BOOST_REQUIRE( &A.sliced(0, 5)[1] == &A[1] );

	BOOST_REQUIRE( A.sliced(0,  0).is_empty() );
	BOOST_REQUIRE( A.sliced(1,  1).is_empty() );
	BOOST_REQUIRE( A.sliced(0, 10).size() == 10 );

	BOOST_REQUIRE(  A[1].sliced(0, 5)[2][3][4] ==  A[1][2][3][4] );
	BOOST_REQUIRE( &A[1].sliced(0, 5)[2][3][4] == &A[1][2][3][4] );

	BOOST_REQUIRE(  A[1].sliced(0, 5)[2] ==  A[1][2] );
	BOOST_REQUIRE( &A[1].sliced(0, 5)[2] == &A[1][2] );

	BOOST_REQUIRE( A[1].sliced(0,  0).is_empty() );
	BOOST_REQUIRE( A[1].sliced(1,  1).is_empty() );
	BOOST_REQUIRE( A[1].sliced(0, 20).size() == 20 );

	BOOST_REQUIRE(  (A << 1).sliced(0, 5)[1][2][3][4] ==  (A << 1)[1][2][3][4] );
	BOOST_REQUIRE( &(A << 1).sliced(0, 5)[1][2][3][4] == &(A << 1)[1][2][3][4] );
}

