#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi rotate"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_rotate_3d){
	multi::array<double, 3> A({3, 4, 5});
	BOOST_REQUIRE(( sizes(A) == decltype(sizes(A)){3, 4, 5} ));

	auto&& RA = rotated(A);
	BOOST_REQUIRE(( sizes(RA) == decltype(sizes(RA)){4, 5, 3} ));
	BOOST_REQUIRE(  &A[0][1][2] == &RA[1][2][0] );

	auto&& UA = unrotated(A);
	BOOST_REQUIRE(( sizes(UA) == decltype(sizes(UA)){5, 3, 4} ));
	BOOST_REQUIRE( &A[0][1][2] == &UA[2][0][1] );

	auto&& RRA = rotated(RA);
	BOOST_REQUIRE(( sizes(RRA) == decltype(sizes(RRA)){5, 3, 4} ));
	BOOST_REQUIRE( &A[0][1][2] == &RRA[2][0][1] );
}

BOOST_AUTO_TEST_CASE(multi_rotate_4d){
	multi::array<double, 4> original({14, 14, 7, 4});

	auto&& unrotd = original.unrotated();
	BOOST_REQUIRE(( sizes(unrotd) == decltype(sizes(unrotd)){4, 14, 14, 7} ));
	BOOST_REQUIRE( &original[0][1][2][3] == &unrotd[3][0][1][2] );

	auto&& unrotd2 = original.unrotated(2);
	BOOST_REQUIRE(( sizes(unrotd2) == decltype(sizes(unrotd2)){7, 4, 14, 14} ));
	BOOST_REQUIRE( &original[0][1][2][3] == &unrotd2[2][3][0][1] );
}

BOOST_AUTO_TEST_CASE(multi_rotate_4d_op){
	multi::array<double, 4> original({14, 14, 7, 4});

	auto&& unrotd = (original >> 1);
	BOOST_REQUIRE(( sizes(unrotd) == decltype(sizes(unrotd)){4, 14, 14, 7} ));
	BOOST_REQUIRE( &original[0][1][2][3] == &unrotd[3][0][1][2] );

	auto&& unrotd2 = (original >> 2);
	BOOST_REQUIRE(( sizes(unrotd2) == decltype(sizes(unrotd2)){7, 4, 14, 14} ));
	BOOST_REQUIRE( &original[0][1][2][3] == &unrotd2[2][3][0][1] );
}

BOOST_AUTO_TEST_CASE(multi_rotate){
	double a[4][5] {
		{ 0,  1,  2,  3,  4}, 
		{ 5,  6,  7,  8,  9}, 
		{10, 11, 12, 13, 14}, 
		{15, 16, 17, 18, 19}
	};
	double b[4][5];
	multi::array_ref<double, 2> A(&a[0][0], {4, 5});
	multi::array_ref<double, 2> B(&b[0][0], {4, 5});
	rotated(B) = rotated(A);
	BOOST_REQUIRE( B[1][1] == 6  );
	BOOST_REQUIRE( B[2][1] == 11 );
	BOOST_REQUIRE( B[1][2] == 7  );
	BOOST_REQUIRE( (B <<1) == (A <<1) );
	BOOST_REQUIRE( (B<<1)[2][1] == 7 );

	{
		multi::array<double, 2> a = {
			{00, 01},
			{10, 11}
		};
		BOOST_REQUIRE(       a[1][0] == 10 );
		BOOST_REQUIRE( (a <<1)[0][1] == 10 );
		BOOST_REQUIRE( &     a[1][0] == &(a <<1)[0][1] );

		BOOST_REQUIRE( a.transposed()[0][1] == 10 );
		BOOST_REQUIRE( transposed(a)[0][1] == 10 );
		BOOST_REQUIRE( (~a)[0][1] == 10 );
		BOOST_REQUIRE( &a[1][0] == &a.transposed()[0][1] );

		(a<<1)[0][1] = 100;
		BOOST_REQUIRE( a[1][0] == 100 );
	}
	{
		multi::array<double, 3> a({11, 13, 17});
		BOOST_REQUIRE( &a[3][5][7] == &a.transposed()[5][3][7] );
		BOOST_REQUIRE( &a[3][5][7] == &transposed(a)[5][3][7] );
		BOOST_REQUIRE( &a[3][5][7] == &(~a)[5][3][7] );

		BOOST_REQUIRE( & a[3][5][7] == & a[3].transposed()[7][5] );	
		BOOST_REQUIRE( & a[3][5][7] == & (~a[3])[7][5] );

		BOOST_REQUIRE( & a[3][5] == & (~a)[5][3] );
		BOOST_REQUIRE( ~~a == a );
		BOOST_REQUIRE( &~~a == &a );
		BOOST_REQUIRE( &(a <<3) == &a);
		BOOST_REQUIRE( &a == &(a<<3) );
		BOOST_REQUIRE( &(a <<1) != &a);
		BOOST_REQUIRE( (a >>1 <<1) == a );
		BOOST_REQUIRE( &(a >>1 <<1) == &a );
	}
	{
		multi::array<double, 2> const a = {
			{00, 01},
			{10, 11}
		};
		BOOST_REQUIRE( (a<<1)[0][1] == 10 );
		BOOST_REQUIRE( &(a<<1)[1][0] == &a[0][1] );
		BOOST_REQUIRE( &(~a)[1][0] == &a[0][1] );
	}
	{
		multi::array<double, 3> const A({3, 5, 7});
	}

}

