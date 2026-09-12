// Copyright 2018-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(array_flatted_2d)
	{
		multi::array<int, 2> arr = {
			{0, 1, 2},
			{3, 4, 5},
		};

		BOOST_TEST( arr.flatted()[1] == 1 );
		BOOST_TEST( arr.flatted()[4] == 4 );

		arr.flatted()[4] = 44;

		BOOST_TEST( arr.flatted()[4] == 44 );
	}

	// BOOST_AUTO_TEST_CASE(array_flatted_3d)
	{
		multi::array<double, 3> arr({13, 4, 5});

		BOOST_TEST( arr.size() == 13 );
		// BOOST_TEST( arr.rotated().is_flattable() );

		{
			auto&& arrRFU = arr.rotated().flatted().unrotated();  // TODO(correaa) remove flatted?
			BOOST_TEST( &arrRFU[11][7] == &arr[11][1][2] );
		}
		{
			auto&& arrRFU = (arr.rotated()).flatted().unrotated();
			BOOST_TEST( &arrRFU[11][7] == &arr[11][7/5][7%5] );
		}
	}

	// BOOST_AUTO_TEST_CASE(array_flatted_3d_bis)
	{
		multi::array<double, 3> const arr({13, 4, 5});
		BOOST_TEST( arr.size() == 13 );
		// BOOST_TEST( arr.is_flattable() );
		BOOST_TEST( arr.flatted().size() == 13L*4L );

		auto const& arr2 = arr.rotated();
		BOOST_TEST( arr2.flatted().size() == 20 );
	}

	// BOOST_AUTO_TEST_CASE(empty_array_3D_flatted)
	{
		multi::array<double, 3> const arr;
		// BOOST_TEST( arr.is_flattable() );
		BOOST_TEST( arr.flatted().size() == 0 );
	}

	// BOOST_AUTO_TEST_CASE(empty_array_2D_flatted)
	{
		multi::array<double, 2> const arr;
		// BOOST_TEST( arr.is_flattable() );
		BOOST_TEST( arr.flatted().size() == 0 );
	}

	// flatted strided
	{
		multi::array<int, 2> arr({4, 10}, 99);
		BOOST_TEST( &arr.strided(2).flatted()[5] != &arr[0][4] );
		BOOST_TEST( &arr.strided(2).flatted()[5] == &arr[0][5] );
	}

	return boost::report_errors();
}
