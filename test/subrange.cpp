// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi subrange selection"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<numeric> // iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_range_section) {
{
	multi::array<double, 4> A({10, 20, 30, 40}, 99.);
	std::iota(A.elements().begin(), A.elements().end(), 0.);

	{
		static_assert( decltype( A({0, 10}, {0, 20}, {0, 30}, {0, 40}) )::rank_v == 4 , "!");
		static_assert( decltype( A(      5, {0, 20}, {0, 30}, {0, 40}) )::rank_v == 3 , "!");
		static_assert( decltype( A({0, 10},      10, {0, 30}, {0, 40}) )::rank_v == 3 , "!");
		static_assert( decltype( A({0, 10}, {0, 20},      15, {0, 40}) )::rank_v == 3 , "!");
		static_assert( decltype( A({0, 10}, {0, 20}, {0, 30},      20) )::rank_v == 3 , "!");

		static_assert( decltype( A(      5,       6, {0, 30}, {0, 40}) )::rank_v == 2 , "!");
		static_assert( decltype( A({0, 10},       6,      15, {0, 40}) )::rank_v == 2 , "!");
		static_assert( decltype( A({0, 10}, {0, 20},      15,      20) )::rank_v == 2 , "!");
	}
	{
		auto&& all = A({0, 10}, {0, 20}, {0, 30}, {0, 40});
		BOOST_REQUIRE( &A[1][2][3][4] == &all[1][2][3][4] );
		BOOST_REQUIRE( &A[1][2][3][4] == &A({0, 10}, {0, 20}, {0, 30}, {0, 40})[1][2][3][4] );
	}
	{
		using multi::_;
		auto&& all = A( {0, 10} , {0, 20} );
		BOOST_REQUIRE( &A[1][2][3][4] == &all[1][2][3][4] );
	}
	{
		BOOST_REQUIRE( &A(0, 0, 0, 0) == &A[0][0][0][0] );
	}
	 {
		auto&& sub = A({0, 5}, {0, 10}, {0, 15}, {0, 20});
		BOOST_REQUIRE( &sub[1][2][3][4] == &A[1][2][3][4] );
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

BOOST_AUTO_TEST_CASE(subrange_assignment) {
	multi::array<double, 2> const A = {
		{1., 2., 3., 4.},
		{5., 6., 7., 8.},
		{9., 0., 1., 2.},
		{3., 4., 5., 6.}
	};
	{
		multi::array<double, 2> B = {
			{9., 9., 9.},
			{9., 9., 9.},
			{9., 9., 9.}
		};
		B({0, 3}, {0, 3}) = A({0, 3}, {0, 3});
		BOOST_REQUIRE( B[1][2] == A[1][2] );
	}
	{
		multi::array<double, 2> B = {
			{9., 9., 9.},
			{9., 9., 9.},
			{9., 9., 9.}
		};
		B() = A({0, 3}, {0, 3});
		BOOST_REQUIRE( B[1][2] == A[1][2] );
		BOOST_REQUIRE( B() == A({0, 3}, {0, 3}) );
	}
	 {
		multi::array<double, 2> B = {
			{9., 9., 9.},
			{9., 9., 9.},
			{9., 9., 9.}
		};
		B = A({0, 3}, {0, 3});
		BOOST_REQUIRE( B[1][2] == A[1][2] );
		BOOST_REQUIRE( B == A({0, 3}, {0, 3}) );
	}
}

