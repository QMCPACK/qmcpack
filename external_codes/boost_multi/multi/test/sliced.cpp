// Copyright 2021-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <numeric>  // for std::iota
#include <vector>   // for std::vector
// IWYU pragma: no_include <algorithm>  // for copy

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_array_sliced_empty)
	{
		multi::array<double, 2> const arr({0, 0}, 99.0);
		BOOST_TEST( arr.sliced(0, 0).is_empty() );
		// BOOST_TEST( arr.sliced(1, 1).is_empty() );  // this results in offsetting nullptr
	}

	// BOOST_AUTO_TEST_CASE(multi_array_sliced)
	{
		multi::array<int, 4> arr({10, 20, 30, 40}, 99);
		std::iota(arr.elements().begin(), arr.elements().end(), 0);

		static_assert(decltype(arr.sliced(0, 5))::dimensionality == 4);

		BOOST_TEST(  arr.sliced( 0, 5)[1][2][3][4] ==  arr[1][2][3][4] );
		BOOST_TEST( &arr.sliced( 0, 5)[1][2][3][4] == &arr[1][2][3][4] );

		BOOST_TEST(  arr.sliced( 0, 5)[1] ==  arr[1] );
		BOOST_TEST( &arr.sliced( 0, 5)[1] == &arr[1] );

		BOOST_TEST(  arr.sliced( 0,  0).empty() );
		BOOST_TEST(  arr.sliced( 1,  1).empty() );
		BOOST_TEST(  arr.sliced( 0, 10).size() == 10 );

		BOOST_TEST(  arr[1].sliced(0, 5)[2][3][4] ==  arr[1][2][3][4] );
		BOOST_TEST( &arr[1].sliced(0, 5)[2][3][4] == &arr[1][2][3][4] );

		BOOST_TEST(  arr[1].sliced(0, 5)[2] ==  arr[1][2] );
		BOOST_TEST( &arr[1].sliced(0, 5)[2] == &arr[1][2] );

		BOOST_TEST( arr[1].sliced(0,  0).is_empty() );
		BOOST_TEST( arr[1].sliced(1,  1).is_empty() );
		BOOST_TEST( arr[1].sliced(0, 20).size() == 20 );

		BOOST_TEST(  (arr.rotated()).sliced(0, 5)[1][2][3][4] ==  (arr.rotated())[1][2][3][4] );
		BOOST_TEST( &(arr.rotated()).sliced(0, 5)[1][2][3][4] == &(arr.rotated())[1][2][3][4] );
	}

	// BOOST_AUTO_TEST_CASE(multi_array_stride)
	{
		multi::array<int, 2> arr = {
			{ 10,  20,  30,  40},
			{ 50,  60,  70,  80},
			{ 90, 100, 110, 120},
			{130, 140, 150, 160},
		};
		BOOST_TEST((
		arr.strided(2) == multi::array<int, 2>{
			{ 10,  20,  30,  40},
			{ 90, 100, 110, 120},
		}
	));
	}

	// BOOST_AUTO_TEST_CASE(multi_array_take)
	{
		multi::array<int, 2> arr = {
			{ 10,  20,  30,  40},
			{ 50,  60,  70,  80},
			{ 90, 100, 110, 120},
			{130, 140, 150, 160},
		};
		BOOST_TEST((
		arr.taked(2) == multi::array<int, 2>{
			{ 10,  20,  30,  40},
			{ 50,  60,  70,  80},
		}
	));
	}

	// BOOST_AUTO_TEST_CASE(drop)
	{
		multi::array<double, 2> arr = {
			{ 10,  20,  30,  40},
			{ 50,  60,  70,  80},
			{ 90, 100, 110, 120},
			{130, 140, 150, 160},
		};
		BOOST_TEST((
		arr.dropped(2) == multi::array<double, 2>{
			{ 90, 100, 110, 120},
			{130, 140, 150, 160},
		}
	));
	}

	// slicing 1D array
	{
		std::vector<int> VV = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

		multi::array_ref<int, 1, int*, multi::contiguous_layout<>> const AA({static_cast<multi::ssize_t>(VV.size())}, VV.data());

		BOOST_TEST( AA.nelems() == static_cast<multi::ssize_t>(VV.size()) );
		BOOST_TEST( !AA.is_empty() );

		BOOST_TEST( AA.sliced(2, 9).size() == 7 );
		BOOST_TEST( AA.sliced(2, 9)[0] == 2 );
		BOOST_TEST( AA.sliced(2, 9)[6] == 8 );
	}

	return boost::report_errors();
}
