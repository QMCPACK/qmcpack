#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi range selection"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#if __cplusplus>=201703L
#include<tuple> // std::apply
#else
#include<experimental/tuple> // std::experimental::apply
#endif

#include<numeric>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_range_section){
{
	multi::array<double, 4> A({10, 20, 30, 40}, 99.);
	std::iota(data_elements(A), data_elements(A) + num_elements(A), 0.);

#if defined(__cpp_lib_apply) and __cpp_lib_apply>=201603
	using std::apply;
#else
	using std::experimental::apply;
#endif
	{
		BOOST_REQUIRE( A({0, 10}, {0, 20}, {0, 30}, {0, 40}).dimensionality == 4 );
		BOOST_REQUIRE( A( 5, {0, 20}, {0, 30}, {0, 40}).dimensionality == 3 );
		BOOST_REQUIRE( A({0, 10}, 10, {0, 30}, {0, 40}).dimensionality == 3 );
		BOOST_REQUIRE( A({0, 10}, {0, 20}, 15, {0, 40}).dimensionality == 3 );
		BOOST_REQUIRE( A({0, 10}, {0, 20}, {0, 30}, 20).dimensionality == 3 );

		BOOST_REQUIRE( A( 5, 6, {0, 30}, {0, 40}).dimensionality == 2 );
		BOOST_REQUIRE( A({0, 10}, 6, 15, {0, 40}).dimensionality == 2 );
		BOOST_REQUIRE( A({0, 10}, {0, 30},15, 20).dimensionality == 2 );
	}
	{
		auto&& all = A({0, 10}, {0, 20}, {0, 30}, {0, 40});
		BOOST_REQUIRE( &A[1][2][3][4] == &all[1][2][3][4] );
		BOOST_REQUIRE( &A[1][2][3][4] == &A({0, 10}, {0, 20}, {0, 30}, {0, 40})[1][2][3][4] );
//		BOOST_REQUIRE( &A[1][2][3][4] == &apply(A, extensions(A))[1][2][3][4] );
	}
	{
		using multi::_;
		auto&& all = A( {0, 10} , {0, 20} );//, *_ , *_ );
		BOOST_REQUIRE( &A[1][2][3][4] == &all[1][2][3][4] );
	}
	{
		BOOST_REQUIRE( &A(0, 0, 0, 0) == &A[0][0][0][0] );
	}
	{
		auto&& sub = A({0, 5}, {0, 10}, {0, 15}, {0, 20});
		BOOST_REQUIRE( &sub[1][2][3][4] == &A[1][2][3][4] );
	}
	{
		BOOST_REQUIRE(  A[1][2][3][4] ==  A(1, 2, 3, 4) );
		BOOST_REQUIRE( &A[1][2][3][4] == &A(1, 2, 3, 4) );
		BOOST_REQUIRE( &A[1][2][3][4] == &apply(A, std::array<int, 4>{1, 2, 3, 4}) );
		BOOST_REQUIRE( &A[1][2][3][4] == &apply(A, std::make_tuple(1, 2, 3, 4))    );
	}
}
{
	multi::array<double, 2> A = {
		{1., 2., 3., 4.},
		{5., 6., 7., 8.},
		{9., 0., 1., 2.},
		{3., 4., 5., 6.}
	};
	multi::array<double, 2> B = {
		{91., 92., 93., 94.},
		{95., 96., 97., 98.},
		{99., 90., 91., 92.},
		{93., 94., 95., 96.}
	};

	A({0, 2}, {0, 2}) = B({0, 2}, {0, 2});
	BOOST_REQUIRE( A != B );
	BOOST_REQUIRE( A({0, 2}, {0, 2}) == B({0, 2}, {0, 2}) );
	BOOST_REQUIRE( A[1][1] == 96. );

}
}

