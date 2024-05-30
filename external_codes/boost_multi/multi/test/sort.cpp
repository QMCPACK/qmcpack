// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <algorithm>  // for std::stable_sort
#include <array>
#include <vector>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_1D_partial_order_syntax) {
	multi::array<int, 1> const tt = {1, 1, 1};
	multi::array<int, 1> const uu = {2, 2, 2};

	BOOST_REQUIRE(     tt <  uu   );
	BOOST_REQUIRE( !  (tt >  uu)  );
	BOOST_REQUIRE(     tt <= uu   );
	BOOST_REQUIRE( !  (tt >= uu)  );
	BOOST_REQUIRE( !  (tt == uu)  );
	BOOST_REQUIRE(    (tt != uu)  );
	BOOST_REQUIRE( ! (uu <  tt)  );
	BOOST_REQUIRE(    (uu >  tt)  );
	BOOST_REQUIRE( !  (uu <= tt)  );
	BOOST_REQUIRE(    (uu >= tt)  );
}

#if defined(__cpp_lib_ranges)
BOOST_AUTO_TEST_CASE(sort_2D) {
	multi::array<int, 2> A = {
		{3, 3, 3},
		{2, 2, 2},
		{1, 1, 1},
	};
	BOOST_REQUIRE(! std::ranges::is_sorted(A));

	std::ranges::sort(A);

	BOOST_REQUIRE(  std::ranges::is_sorted(A));

	static_assert(std::permutable<boost::multi::array_iterator<int, 2, int *>>);
}

BOOST_AUTO_TEST_CASE(sort_strings) {
	auto A = multi::array<char, 2>{
		{'S', 'e', 'a', 'n', ' ', ' '},
		{'A', 'l', 'e', 'x', ' ', ' '},
		{'B', 'j', 'a', 'r', 'n', 'e'},
	};
	BOOST_REQUIRE(! std::ranges::is_sorted(A));

	std::ranges::sort(A);

	BOOST_REQUIRE(  std::ranges::is_sorted(A));

	BOOST_REQUIRE((
		A == multi::array<char, 2>{
			{'A', 'l', 'e', 'x', ' ', ' '},
			{'B', 'j', 'a', 'r', 'n', 'e' },
			{'S', 'e', 'a', 'n', ' ', ' '},
		}
	));

	std::ranges::sort(~A);
	BOOST_REQUIRE(  std::ranges::is_sorted(~A));

	static_assert(std::permutable<boost::multi::array_iterator<int, 2, int *>>);
}
#endif

BOOST_AUTO_TEST_CASE(multi_array_stable_sort) {
	std::vector<double> vec = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( std::is_sorted(begin(vec), end(vec)) );

	multi::array<double, 2> d2D = {
		{150.0, 16.0, 17.0, 18.0, 19.0},
		{ 30.0,  1.0,  2.0,  3.0,  4.0},
		{100.0, 11.0, 12.0, 13.0, 14.0},
		{ 50.0,  6.0,  7.0,  8.0,  9.0},
	};
	BOOST_REQUIRE( ! std::is_sorted(begin(d2D), end(d2D) ) );

	std::stable_sort(begin(d2D), end(d2D));
	BOOST_REQUIRE( std::is_sorted( begin(d2D), end(d2D) ) );

	BOOST_REQUIRE((
		d2D == decltype(d2D){
			{ 30.0,  1.0,  2.0,  3.0,  4.0},
			{ 50.0,  6.0,  7.0,  8.0,  9.0},
			{100.0, 11.0, 12.0, 13.0, 14.0},
			{150.0, 16.0, 17.0, 18.0, 19.0},
		}
	));

	BOOST_REQUIRE( ! std::is_sorted( begin(d2D.rotated()), end(d2D.rotated()) ) );

	std::stable_sort(begin(d2D.rotated()), end(d2D.rotated()));
	BOOST_REQUIRE( std::is_sorted( begin(d2D.rotated()), end(d2D.rotated()) ) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D          ), end(d2D          ) ) );

	BOOST_REQUIRE((
		d2D == decltype(d2D){
			{ 1.0,  2.0,  3.0,  4.0,  30.0},
			{ 6.0,  7.0,  8.0,  9.0,  50.0},
			{11.0, 12.0, 13.0, 14.0, 100.0},
			{16.0, 17.0, 18.0, 19.0, 150.0},
		}
	));
}

BOOST_AUTO_TEST_CASE(multi_array_ref_stable_sort) {
	std::vector<double> vec = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( std::is_sorted(begin(vec), end(vec)) );

	// clang-format off
	std::array<std::array<double, 5>, 4> d2D {{
		{{150.0, 16.0, 17.0, 18.0, 19.0}},
		{{ 30.0,  1.0,  2.0,  3.0,  4.0}},
		{{100.0, 11.0, 12.0, 13.0, 14.0}},
		{{ 50.0,  6.0,  7.0,  8.0,  9.0}}
	}};
	// clang-format on

	auto&& d2D_ref = *multi::array_ptr<double, 2>(&d2D[0][0], {4, 5});  // NOLINT(readability-container-data-pointer) test access

	BOOST_REQUIRE( ! std::is_sorted(begin(d2D_ref), end(d2D_ref) ) );
	std::stable_sort(begin(d2D_ref), end(d2D_ref));
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref), end(d2D_ref) ) );

	BOOST_REQUIRE( ! std::is_sorted( begin(d2D_ref.rotated()), end(d2D_ref.rotated()) ) );
	std::stable_sort(begin(d2D_ref.rotated()), end(d2D_ref.rotated()));
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref.rotated()), end(d2D_ref.rotated()) ) );
}

BOOST_AUTO_TEST_CASE(lexicographical_compare) {
	multi::array<char, 1> const name1 = {'a', 'b', 'c'};
	multi::array<char, 1> const name2 = {'a', 'c', 'c'};
	BOOST_REQUIRE(name1 != name2 );
	BOOST_REQUIRE(name1 < name2);
	BOOST_REQUIRE(name1 <= name2);
	BOOST_REQUIRE(!(name1 > name2));
	BOOST_REQUIRE(!(name1 > name2));
}

BOOST_AUTO_TEST_CASE(lexicographical_compare_offset) {
	multi::array<char, 1> const name1 = {'a', 'b', 'c'};
	multi::array<char, 1> name2({{1, 4}}, '\0');

	BOOST_REQUIRE(  name2.size() == 3 );
	BOOST_REQUIRE(( name2.extension() == multi::extension_t<multi::index>{1, 4} ));
	BOOST_REQUIRE(( name2.extension() == multi::extension_t{multi::index{1}, multi::index{4}} ));

	// BOOST_REQUIRE(( name2.extension() == multi::extension_t{1L, 4L} ));

	BOOST_REQUIRE(( name2.extension() == multi::extension_t<>{1, 4} ));
	// BOOST_REQUIRE(( name2.extension() == multi::extension_t{1 , 4 } )); TODO(correaa) solve ambiguity

	name2[1] = 'a';
	name2[2] = 'b';
	name2[3] = 'c';

	BOOST_REQUIRE(  name2 != name1 );
	BOOST_REQUIRE(!(name2 == name1));

	BOOST_REQUIRE(  name2 <  name1 );
	BOOST_REQUIRE(  name2 <= name1 );

	BOOST_REQUIRE(!(name2 >  name1));
	BOOST_REQUIRE(!(name2 >= name1));

	// BOOST_REQUIRE(!(name1 > name2));
	// BOOST_REQUIRE(!(name1 > name2));
}

BOOST_AUTO_TEST_CASE(lexicographical_compare_offset_2d) {
	multi::array<char, 2> const name1 = {{'a', 'b'}, {'b', 'c'}, {'c', 'd'}};
	multi::array<char, 2> name2({{1, 4}, {0, 2}}, '\0');

	BOOST_REQUIRE(  name2.size() == 3 );
	BOOST_REQUIRE(( name2.extension() == multi::extension_t<multi::index>{1, 4} ));
	BOOST_REQUIRE(( name2.extension() == multi::extension_t<>{1, 4} ));
	// BOOST_REQUIRE(( name2.extension() == multi::extension_t{1 , 4 } )); TODO(correaa) solve ambiguity

	name2[1][0] = 'a';  name2[1][1] = 'a';
	name2[2][0] = 'b';  name2[2][1] = 'a';
	name2[3][0] = 'c';  name2[3][1] = 'a';

	BOOST_REQUIRE(  name2 != name1 );
	BOOST_REQUIRE(!(name2 == name1));

	BOOST_REQUIRE(  name2 <  name1 );
	BOOST_REQUIRE(  name2 <= name1 );

	// BOOST_REQUIRE(!(name2 >  name1));
	// BOOST_REQUIRE(!(name2 >= name1));

	BOOST_REQUIRE( name1 > name2 );
	BOOST_REQUIRE(!(name1 < name2));
}
