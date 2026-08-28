// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for range, layout_t, get, extents_t

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		multi::array<int, 2> const arr = {
			{ 0,  1,  2,  3,  4},
			{ 5,  6,  7,  8,  9},
			{10, 11, 12, 13, 14}
		};

		auto const& barr = arr.flattened();

		BOOST_TEST( barr.size() == 15 );

		auto it = barr.begin();

		BOOST_TEST( &(*(it + 0)) == &arr[0][0] );
		BOOST_TEST( &(*(it + 1)) == &arr[0][1] );
		BOOST_TEST( &(*(it + 2)) == &arr[0][2] );
		BOOST_TEST( &(*(it + 3)) == &arr[0][3] );
		BOOST_TEST( &(*(it + 4)) == &arr[0][4] );

		BOOST_TEST( &(*(it + 5)) == &arr[1][0] );
		BOOST_TEST( &(*(it + 6)) == &arr[1][1] );

		BOOST_TEST( &(*(it + 10)) == &arr[2][0] );
		BOOST_TEST( &(*(it + 11)) == &arr[2][1] );

		BOOST_TEST( &(*it) == &arr[0][0] );
		++it;
		BOOST_TEST( &(*it) == &arr[0][1] );
		++it;
		BOOST_TEST( &(*it) == &arr[0][2] );
		it += 2;
		BOOST_TEST( &(*it) == &arr[0][4] );

		++it;
		BOOST_TEST( &(*(it)) == &arr[1][0] );
		++it;
		BOOST_TEST( &(*(it)) == &arr[1][1] );

		it += 7;
		BOOST_TEST( &(*(it)) == &arr[2][3] );

		auto it2 = barr.begin();
		++it2;
		it2 += 12;
		BOOST_TEST( &(*(it2)) == &arr[2][3] );
	}

	// BOOST_AUTO_TEST_CASE(layout_3)
	{
		multi::array<double, 2> arr({50, 50});
		BOOST_TEST( arr.size() == 50 );

		BOOST_TEST( arr[0].sliced(10, 20).size() == 10 );
		BOOST_TEST( size(arr[0].sliced(10, 20))  == 10 );

		static_assert(decltype(arr(0, {10, 20}))::rank_v == 1);

		BOOST_TEST( size(arr(0, {10, 20})) == 10 );

		BOOST_TEST(      arr.layout() == arr.layout()  );
		BOOST_TEST( !(arr.layout() <  arr.layout()) );

		auto const& barr = arr.flattened();
		BOOST_TEST( &barr[10] == &arr[0][10] );
	}

	{
		multi::array<double, 2> arr({3, 5});

		auto&& barr = arr.flattened();

		BOOST_TEST( &barr [0] == &arr[0][0] );

		BOOST_TEST( &*barr.begin() == &barr[0] );
		BOOST_TEST( &*(barr.begin() + 1) == &barr[1] );
	}

	{
		multi::array<int, 2> arr = {
			{ 0,  1,  2,  3,  4},
			{ 5,  6,  7,  8,  9},
			{10, 11, 12, 13, 14}
		};

		BOOST_TEST( arr.size() == 3 );

		auto&& barr = arr().flattened();

		BOOST_TEST( &barr [0] == &arr[0][0] );

		BOOST_TEST( &*barr.begin() == &barr[0] );
		BOOST_TEST( &*(barr.begin() + 1) == &barr[1] );

		auto it = barr.begin();
		BOOST_TEST( &*(it + 0) == &arr[0][0] );
		BOOST_TEST( &*(it + 1) == &arr[0][1] );
		BOOST_TEST( &*(it + 2) == &arr[0][2] );
		BOOST_TEST( &*(it + 3) == &arr[0][3] );
		BOOST_TEST( &*(it + 4) == &arr[0][4] );

		BOOST_TEST( &*(it + 5) == &arr[1][0] );
		BOOST_TEST( &*(it + 6) == &arr[1][1] );

		BOOST_TEST( &*it == &arr[0][0] );
		++it;
		BOOST_TEST( &*it == &arr[0][1] );
		++it;
		BOOST_TEST( &*it == &arr[0][2] );
		++it;
		BOOST_TEST( &*it == &arr[0][3] );
		++it;
		BOOST_TEST( &*it == &arr[0][4] );
		++it;
		BOOST_TEST( &*it == &arr[1][0] );
	}
	{
		multi::array<int, 2> arr_original = {
			{ 0,  1,  2,  3,  4, 99},
			{ 5,  6,  7,  8,  9, 99},
			{10, 11, 12, 13, 14, 99},
			{99, 99, 99, 99, 99, 99}
		};

		auto&& arr = arr_original({0, 3}, {0, 5});

		multi::array<int, 2>::iterator const it2(arr.base(), arr.layout().sub(), arr.layout().stride());
		BOOST_TEST(( *it2 == multi::array<int, 1>{0, 1, 2, 3, 4} ));

		multi::array<int, 1>::iterator const it1(arr.base(), arr.layout().sub().sub(), arr.layout().sub().stride());

		// multi::array_iterator<int, 1, multi::array_iterator<int, 1, int*> > const bit({arr.base(), arr.base(), arr.base()}, arr.layout().sub(), arr.stride());

		// multi::detail::what(arr.begin());

		BOOST_TEST( arr.size() == 3 );

		auto&& barr = arr().flattened();

		BOOST_TEST( &barr [0] == &arr[0][0] );

		BOOST_TEST( &*barr.begin() == &barr[0] );
		BOOST_TEST( &*(barr.begin() + 1) == &barr[1] );

		auto it = barr.begin();
		BOOST_TEST( &*(it + 0) == &arr[0][0] );
		BOOST_TEST( &*(it + 1) == &arr[0][1] );
		BOOST_TEST( &*(it + 2) == &arr[0][2] );
		BOOST_TEST( &*(it + 3) == &arr[0][3] );
		BOOST_TEST( &*(it + 4) == &arr[0][4] );

		BOOST_TEST( &*(it + 5) == &arr[1][0] );
		BOOST_TEST( &*(it + 6) == &arr[1][1] );

		BOOST_TEST( &*it == &arr[0][0] );
		++it;
		BOOST_TEST( &*it == &arr[0][1] );
		++it;
		BOOST_TEST( &*it == &arr[0][2] );
		++it;
		BOOST_TEST( &*it == &arr[0][3] );
		++it;
		BOOST_TEST( &*it == &arr[0][4] );
		++it;
		BOOST_TEST( &*it == &arr[1][0] );
	}
	{
		multi::array<int, 1> const arr({5});

		// auto const& barr = arr.flattened();  // compilation error, good
	}
	{
		multi::array<int, 3> const arr({3, 5, 7});

		auto const& barr = arr.flattened();

		BOOST_TEST( barr.size() == 15 );

		auto it = barr.begin();

		BOOST_TEST( &(*(it + 0))[0] == &arr[0][0][0] );
		BOOST_TEST( &(*(it + 1))[0] == &arr[0][1][0] );
		BOOST_TEST( &(*(it + 2))[0] == &arr[0][2][0] );
		BOOST_TEST( &(*(it + 3))[0] == &arr[0][3][0] );
		BOOST_TEST( &(*(it + 4))[0] == &arr[0][4][0] );

		BOOST_TEST( &(*(it + 5))[0] == &arr[1][0][0] );

		BOOST_TEST( &(*(it + 10))[0] == &arr[2][0][0] );
	}
	{
		multi::array<int, 2> const arr = {
			{ 0,  1,  2,  3,  4, 55},
			{ 5,  6,  7,  8,  9, 99},
			{10, 11, 12, 13, 14, 44},
			{44, 55, 66, 77, 88, 99}
		};

		auto const& barr = arr({0, 3}, {0, 5}).flattened();

		BOOST_TEST( barr.size() == 15 );

		{
			auto it = barr.begin();

			BOOST_TEST( &(*(it + 0)) == &arr[0][0] );
			BOOST_TEST( &(*(it + 1)) == &arr[0][1] );
			BOOST_TEST( &(*(it + 2)) == &arr[0][2] );
			BOOST_TEST( &(*(it + 3)) == &arr[0][3] );
			BOOST_TEST( &(*(it + 4)) == &arr[0][4] );

			BOOST_TEST( &(*(it + 5)) == &arr[1][0] );
			BOOST_TEST( &(*(it + 6)) == &arr[1][1] );

			BOOST_TEST( &(*(it + 10)) == &arr[2][0] );
			BOOST_TEST( &(*(it + 11)) == &arr[2][1] );

			// auto&& sgm = it.segment();
		}
		{
			auto it = barr.begin();

			BOOST_TEST( &(*it) == &arr[0][0] );

			++it;
			BOOST_TEST( &(*it) == &arr[0][1] );
			++it;
			BOOST_TEST( &(*it) == &arr[0][2] );
			++it;
			BOOST_TEST( &(*it) == &arr[0][3] );
			++it;
			BOOST_TEST( &(*it) == &arr[0][4] );

			++it;
			BOOST_TEST( &(*it) == &arr[1][0] );

			++it;
			BOOST_TEST( &(*it) == &arr[1][1] );

			++it;
			BOOST_TEST( &(*it) == &arr[1][2] );
		}
	}
	{
		multi::array<double, 2> arr({6, 10});

		auto const& barr = arr.strided(2).flattened();

		BOOST_TEST( &barr [0] == &arr[0][0] );
		BOOST_TEST( &barr [1] == &arr[0][1] );
		// ...
		BOOST_TEST( &barr [9] == &arr[0][9] );

		BOOST_TEST( &barr[10] == &arr[2][0] );
		BOOST_TEST( &barr[11] == &arr[2][1] );
		BOOST_TEST( &barr[12] == &arr[2][2] );
		// ...
		BOOST_TEST( &barr[19] == &arr[2][9] );

		BOOST_TEST( &barr[20] == &arr[4][0] );
		BOOST_TEST( &barr[21] == &arr[4][1] );
		BOOST_TEST( &barr[22] == &arr[4][2] );
		// ...
		BOOST_TEST( &barr[29] == &arr[4][9] );

		BOOST_TEST( arr.num_elements() == 60 );
		BOOST_TEST( barr.size() == 30 );
	}
	{
		multi::array<double, 2> arr({6, 10});

		auto const& barr = arr.strided(2).transposed().strided(2).transposed().flattened();

		BOOST_TEST( &barr [0] == &arr[0][0] );
		BOOST_TEST( &barr [1] == &arr[0][2] );
		BOOST_TEST( &barr [2] == &arr[0][4] );
		BOOST_TEST( &barr [3] == &arr[0][6] );
		BOOST_TEST( &barr [4] == &arr[0][8] );

		BOOST_TEST( &barr [5] == &arr[2][0] );
		BOOST_TEST( &barr [6] == &arr[2][2] );
		BOOST_TEST( &barr [7] == &arr[2][4] );
		BOOST_TEST( &barr [8] == &arr[2][6] );
		BOOST_TEST( &barr [9] == &arr[2][8] );

		BOOST_TEST( &barr [10] == &arr[4][0] );
		BOOST_TEST( &barr [11] == &arr[4][2] );
		BOOST_TEST( &barr [12] == &arr[4][4] );
		BOOST_TEST( &barr [13] == &arr[4][6] );
		BOOST_TEST( &barr [14] == &arr[4][8] );

		BOOST_TEST( arr.num_elements() == 60 );
		BOOST_TEST( barr.size() == 15 );
	}
	{
		multi::array<int, 2> const arr_original = {
			{ 0,  1,  2,  3,  4, 99},
			{ 5,  6,  7,  8,  9, 99},
			{10, 11, 12, 13, 14, 99},
			{99, 99, 99, 99, 99, 99}
		};

		auto&& arr = arr_original({0, 3}, {0, 5});
		BOOST_TEST( arr.elements()[2] == 2 );
		// BOOST_TEST( arr.elements()[3] = 22 );  // doesn't compile, good
	}
	{
		multi::array<double, 3> arr = {
			{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}},
			{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}},
			{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}},
			{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}
		};

		auto [ni, nj, nk] = arr.sizes();
		BOOST_TEST( ni == 4 );
		BOOST_TEST( nj == 2 );
		BOOST_TEST( nk == 3 );

		BOOST_TEST( arr({0, 3}).flattened().size() == 6 );
		BOOST_TEST( (*arr({0, 3}).flattened().begin()).size() == 3 );
		BOOST_TEST( arr({0, 3}).flattened()[0].size() == 3 );
	}

	return boost::report_errors();
}
