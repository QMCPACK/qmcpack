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
	BOOST_REQUIRE( arr.is_flattable() );
	BOOST_REQUIRE( arr.flatted().size() == 52 );
}

BOOST_AUTO_TEST_CASE(empty_array_3D_flatted) {
	multi::array<double, 3> const arr;
	BOOST_REQUIRE( arr.is_flattable() );
	BOOST_REQUIRE( arr.flatted().size() == 0 );
}

BOOST_AUTO_TEST_CASE(empty_array_2D_flatted) {
	multi::array<double, 2> const arr;
	BOOST_REQUIRE( arr.is_flattable() );
	BOOST_REQUIRE( arr.flatted().size() == 0 );
}
