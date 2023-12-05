// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi allocators"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

//#include "../../multi/memory/stack.hpp" //TODO(correaa) test custom allocator

#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(std_vector_of_arrays) {
	std::vector<multi::array<double, 2>> va;
	std::transform(
		begin(multi::iextension(3)), end(multi::iextension(3)),
		std::back_inserter(va),
		[](auto idx){return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(va[0]) == 0 );
	BOOST_REQUIRE( size(va[1]) == 1 );
	BOOST_REQUIRE( size(va[2]) == 2 );
	BOOST_REQUIRE( va[1] [0][0] == 1 );
	BOOST_REQUIRE( va[2] [0][0] == 2 );

	std::vector<multi::array<double, 2>> wa = {
		multi::array<double, 2>({0, 0}, 0.),
		multi::array<double, 2>({1, 1}, 1.),
		multi::array<double, 2>({2, 2}, 2.),
	};
	BOOST_REQUIRE( size(va) == size(wa) );
	BOOST_REQUIRE( va == wa );

	std::vector<multi::array<double, 2>> ua(3);
	auto iex = multi::iextension(static_cast<multi::size_type>(ua.size()));
	std::transform(
		begin(iex), end(iex),
		begin(ua),
		[](auto idx){return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);
	BOOST_REQUIRE( ua == va );
}

BOOST_AUTO_TEST_CASE(array1d_of_arrays2d) {
	multi::array<multi::array<double, 2>, 1> arr(multi::extensions_t<1>{multi::iextension{10}}, multi::array<double, 2>{});
	BOOST_REQUIRE( size(arr) == 10 );

	std::transform(
		begin(extension(arr)), end(extension(arr)), begin(arr),
		[](auto idx) {return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(arr[0]) == 0 );
	BOOST_REQUIRE( size(arr[1]) == 1 );
	BOOST_REQUIRE( size(arr[8]) == 8 );
	BOOST_REQUIRE( arr[8][4][4] == 8 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d)  {
	multi::array<multi::array<double, 3>, 2> AA({10, 20}, multi::array<double, 3>{});
	for(int i = 0; i != 10; ++i) {
		for(int j = 0; j != 20; ++j) {
			AA[i][j] = multi::array<double, 3>({i+j, i+j, i+j}, 99.);
		}
	}
	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99. );
}

BOOST_AUTO_TEST_CASE(array_3d_with_hint_int) {
	multi::array<double, 2> const arr({3, 4});
	multi::array<int, 3> arr_hint({3, 4, 5}, arr.cbase());

	arr_hint[1][2][3] = 4;
	BOOST_REQUIRE( size(arr_hint) == 3 );
	BOOST_REQUIRE( arr_hint[1][2][3] == 4 );

	multi::array<int, 3> arr2({3, 4, 5}, 0);
	BOOST_REQUIRE( size(arr2) == 3 );

	multi::array<int, 3> arr3({3, 4, 5}, 99);
	BOOST_REQUIRE( size(arr3) == 3 );
}

BOOST_AUTO_TEST_CASE(array_3d_with_hint_size_t) {
	multi::array<double, 2> const arr({3, 4});
	multi::array<size_t, 3> arr_hint({3, 4, 5}, arr.cbase());

	arr_hint[1][2][3] = 4;
	BOOST_REQUIRE( size(arr_hint) == 3 );
	BOOST_REQUIRE( arr_hint[1][2][3] == 4 );

	multi::array<size_t, 3> arr3({3, 4, 5}, 0UL);
	BOOST_REQUIRE( size(arr3) == 3 );

	multi::array<size_t, 3> arr4({3, 4, 5}, 99);
	BOOST_REQUIRE( size(arr4) == 3 );
}
