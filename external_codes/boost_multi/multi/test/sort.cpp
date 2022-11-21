// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "Unit Tests for Multi sort"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<algorithm> // stable_sort
#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_stable_sort) {
	std::vector<double> vec = {1., 2., 3.};
	BOOST_REQUIRE( std::is_sorted(begin(vec), end(vec)) );

	multi::array<double, 2> d2D = {
		{150, 16, 17, 18, 19},
		{ 30,  1,  2,  3,  4},
		{100, 11, 12, 13, 14},
		{ 50,  6,  7,  8,  9}
	};
	BOOST_REQUIRE( not std::is_sorted(begin(d2D), end(d2D) ) );

	std::stable_sort( begin(d2D), end(d2D) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D), end(d2D) ) );

	BOOST_REQUIRE((
		d2D == decltype(d2D){
			{30, 1, 2, 3, 4},
			{50, 6, 7, 8, 9},
			{100, 11, 12, 13, 14},
			{150, 16, 17, 18, 19}
		}
	));

	BOOST_REQUIRE( not std::is_sorted( begin(d2D.rotated()), end(d2D.rotated()) ) );

	std::stable_sort( begin(d2D.rotated()), end(d2D.rotated()) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D.rotated()), end(d2D.rotated()) ) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D          ), end(d2D          ) ) );

	BOOST_REQUIRE((
		d2D == decltype(d2D){
			{1, 2, 3, 4, 30},
			{6, 7, 8, 9, 50},
			{11, 12, 13, 14, 100},
			{16, 17, 18, 19, 150}
		}
	));
}

BOOST_AUTO_TEST_CASE(multi_array_ref_stable_sort) {
	std::vector<double> vec = {1., 2., 3.};
	BOOST_REQUIRE( std::is_sorted(begin(vec), end(vec)) );

	std::array<std::array<double, 5>, 4> d2D {{
		{{150, 16, 17, 18, 19}},
		{{ 30,  1,  2,  3,  4}},
		{{100, 11, 12, 13, 14}},
		{{ 50,  6,  7,  8,  9}}
	}};
	auto&& d2D_ref = *multi::array_ptr<double, 2>(&d2D[0][0], {4, 5});  // NOLINT(readability-container-data-pointer) test access

	BOOST_REQUIRE( not std::is_sorted(begin(d2D_ref), end(d2D_ref) ) );
	std::stable_sort( begin(d2D_ref), end(d2D_ref) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref), end(d2D_ref) ) );

	BOOST_REQUIRE( not std::is_sorted( begin(d2D_ref.rotated()), end(d2D_ref.rotated()) ) );
	std::stable_sort( begin(d2D_ref.rotated()), end(d2D_ref.rotated()) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref.rotated()), end(d2D_ref.rotated()) ) );
}
