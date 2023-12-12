// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_flatted_3d) {
	multi::array<double, 3> arr({13, 4, 5});

	BOOST_REQUIRE( arr.size() == 13 );
	BOOST_REQUIRE( arr.rotated().is_flattable() );

	{
		auto&& arrRFU = arr.rotated().flatted().unrotated();
		BOOST_REQUIRE( &arrRFU[11][7] == &arr[11][1][2] );
	}
	{
		auto&& arrRFU = (arr.rotated()).flatted().unrotated();
		BOOST_REQUIRE( &arrRFU[11][7] == &arr[11][7/5][7%5] );
	}
}

BOOST_AUTO_TEST_CASE(array_flatted_3d_bis) {
	multi::array<double, 3> const arr({13, 4, 5});
	BOOST_REQUIRE( arr.size() == 13 );
}
