#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi select range"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_range_section_1D){
	multi::array<double, 1> A = {00., 01., 02.}; (void)A;
	BOOST_REQUIRE( A == A(multi::all) );
	BOOST_REQUIRE( size(A( 1 <= multi::all )) == 2 );
	BOOST_REQUIRE( A( 1 <= multi::all )[0] == 1. );
	BOOST_REQUIRE( size(A( multi::all < 2 )) == 2 );
	BOOST_REQUIRE( A( multi::all < 2 )[1] == 1. );
}

BOOST_AUTO_TEST_CASE(multi_array_range_section){
	multi::array<double, 2> A = {
		{00., 01., 02.},
		{10., 11., 12.},
		{20., 21., 22.},
		{30., 31., 32.},
	};

	BOOST_REQUIRE( size( A( multi::all, 2) ) == size(A) );
	BOOST_REQUIRE( size( A( multi::_ , 2) ) == size(A) );
	BOOST_REQUIRE( size( A( *multi::_ , 2) ) == size(A) );
	BOOST_REQUIRE( size( A( multi::U , 2) ) == size(A) );

	BOOST_REQUIRE( size( A(   multi::all  , 2) ) == 4 );
	BOOST_REQUIRE( size( A(   multi::all<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(1<=multi::all  , 2) ) == 3 );
	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=multi::all<3, 2) ) == 2 );

	BOOST_REQUIRE( size( A(   multi::_  , 2) ) == 4 );
	BOOST_REQUIRE( size( A(   multi::_<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(1<=multi::_  , 2) ) == 3 );
	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=multi::_<3, 2) ) == 2 );

	using multi::_;

	BOOST_REQUIRE( size( A(   _  , 2) ) == 4 );
	BOOST_REQUIRE( size( A(   _<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(1<=_  , 2) ) == 3 );
	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=_<3, 2) ) == 2 );

	using multi::U;
	BOOST_REQUIRE( size( A(_, 2) ) == size(A) );
	BOOST_REQUIRE( size( A(*_, 2) ) == size(A) );

	BOOST_REQUIRE( size( A(_<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(*_<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(U<2, 2) ) == 2 );

	BOOST_REQUIRE( size( A(1<=_, 2) ) == 3 );
	BOOST_REQUIRE( size( A(1<=*_, 2) ) == 3 );
	BOOST_REQUIRE( size( A(1<=U, 2) ) == 3 );

	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=_<3, 2) ) == 2 );
	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=*_<3, 2) ) == 2 );
	// cppcheck-suppress compareBoolExpressionWithInt ; because DSL
	BOOST_REQUIRE( size( A(1<=U<3, 2) ) == 2 );

	BOOST_REQUIRE( size( A(*_<2, 2) ) == 2 );
	BOOST_REQUIRE( size( A(U<2, 2) ) == 2 );

	BOOST_REQUIRE( size( A(A.extension(), 2) ) == size(A) );

	auto&& col2( A(A.extension(0), 2) ); // select column #2
	// same as A(extesion(A), 2)
	// same as A(A.extension(0), 2);
	// same as rotated(A)[2];
	BOOST_REQUIRE( col2.size(0) == size(A) );

	BOOST_REQUIRE( dimensionality(col2) == 1 );
	BOOST_REQUIRE( size(col2) == size(A) );
	BOOST_REQUIRE( col2.size() == size(A) );
	BOOST_REQUIRE( stride(col2) == 3 );
	BOOST_REQUIRE( col2[0] == 02. );
	BOOST_REQUIRE( col2[1] == 12. );
	BOOST_REQUIRE(( col2 == multi::array<double, 1>{02., 12., 22., 32.} ));
	BOOST_REQUIRE(( col2 == multi::array<double, 1>(rotated(A)[2]) ));
	BOOST_REQUIRE(( col2 == rotated(A)[2] ));
	BOOST_REQUIRE(( col2 == A(A.extension(0), 2) ));
}

