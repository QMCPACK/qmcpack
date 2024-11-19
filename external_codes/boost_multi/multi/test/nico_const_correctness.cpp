// Copyright 2022-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for subarray, array, range, operator!=

#include <algorithm>    // for fill, copy, for_each
#include <iostream>     // for operator<<, basic_ostream::opera...
#include <iterator>     // for begin, end, ostream_iterator
#include <type_traits>  // for decay_t
#include <utility>      // for forward

namespace multi = boost::multi;

template<class Array1D>
void print(Array1D const& coll) {
	// *coll.begin() = 99;  // doesn't compile "assignment of read-only location"

	std::copy(std::begin(coll), std::end(coll), std::ostream_iterator<typename Array1D::value_type>(std::cout, ", "));
	std::cout << '\n';
}

template<class Array1D>
auto fill_99(Array1D&& col) -> Array1D&& {
	std::fill(std::begin(col), std::end(col), 99);
	return std::forward<Array1D>(col);
}

template<class Array2D>
void print_2d(Array2D const& coll) {
	// *(coll.begin()->begin()) = 99;  // doesn't compile "assignment of read-only location"

	std::for_each(std::begin(coll), std::end(coll), [](auto const& row) {
		std::copy(std::begin(row), std::end(row), std::ostream_iterator<typename Array2D::element_type>(std::cout, ", "));
		std::cout << '\n';
	});
}

template<class Array1D>
auto fill_2d_99(Array1D&& coll) -> Array1D&& {
	// for(auto const& row : coll) {  // does not work because it would make it const
	std::for_each(std::begin(coll), std::end(coll), [](typename std::decay_t<Array1D>::reference row) {
		std::fill(std::begin(row), std::end(row), 99);
	});
	// std::transform(coll.begin(), coll.end(), coll.begin(), [](auto&& row) {
	//  std::fill(row.begin(), row.end(), 99);
	//  return std::forward<decltype(row)>(row);
	// });
	return std::forward<Array1D>(coll);
}

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(const_views) {
		multi::array<int, 1> coll1 = {0, 8, 15, 47, 11, 42};
		print(coll1);  // prints "0, 8, 15, 47, 11, 42"

		print(coll1({0, 3}));  // similar to coll1 | take(3) // prints "0, 8, 15"

		auto&& coll1_take3 = coll1({0, 3});
		print(coll1_take3);  // prints "0, 8, 15"
	}

	BOOST_AUTO_TEST_CASE(mutating_views) {
		multi::array<int, 1> coll1 = {0, 8, 15, 47, 11, 42};

		fill_99(coll1);
		fill_99(coll1({0, 3}));

		auto&& coll1_take3 = coll1({0, 3});
		fill_99(coll1_take3);

		auto const& coll2 = coll1;
		// fill_99( coll2 );  // doesn't compile because coll2 is const ("assignment of read-only" inside fill_99)
		// fill_99( coll2({0, 3}) );  // similar to coll2 | take(3) doesn't compile ("assignment of read-only")

		auto const& coll1_take3_const = coll1({0, 3});
		// fill_99( coll1_take3_const );  // doesn't compile because coll1_take3_const is const ("assignment of read-only")

		(void)coll2, (void)coll1_take3_const, (void)coll1_take3;
	}

	BOOST_AUTO_TEST_CASE(const_views_2d) {
		multi::array<int, 2> coll1 = {
			{0, 8, 15, 47, 11, 42},
			{0, 8, 15, 47, 11, 42},
		};

		print_2d(coll1);  // prints "0, 8, 15, 47, 11, 42"

		print_2d(coll1({0, 2}, {0, 3}));  // similar to coll1 | take(3) // prints "0, 8, 15"

		auto&& coll1_take3 = coll1({0, 2}, {0, 3});
		print_2d(coll1_take3);  // prints "0, 8, 15"
	}

	BOOST_AUTO_TEST_CASE(mutating_views_2d) {
		multi::array<int, 2> coll1 = {
			{0, 8, 15, 47, 11, 42},
			{0, 8, 15, 47, 11, 42},
		};

		fill_2d_99(coll1);
		fill_2d_99(coll1({0, 2}, {0, 3}));

		auto&& coll1_take3 = coll1({0, 2}, {0, 3});
		fill_2d_99(coll1_take3);

		auto const& coll2 = coll1;
		// fill_99( coll2 );  // doesn't compile because coll2 is const ("assignment of read-only" inside fill_99)
		// fill_99( coll2({0, 3}) );  // similar to coll2 | take(3) doesn't compile ("assignment of read-only")

		auto const& coll1_take3_const = coll1({0, 2}, {0, 3});
		// fill_99( coll1_take3_const );  // doesn't compile because coll1_take3_const is const ("assignment of read-only")

		(void)coll2, (void)coll1_take3_const, (void)coll1_take3;
	}

	{
		multi::array<int, 1> arr1d = {1, 2, 3};

		// multi::array<int, 1>::const_iterator cfirst = arr1d.cbegin();
		// *cfirst.base() = 5;
		// *cfirst = 5;  // correctly fails to compile
		// cfirst[0] = 5;  // correctly fails to compile

		BOOST_TEST( arr1d[0] == 1 );
	}
	{
		multi::array<int, 1> const arr1d = {1, 2, 3};

		// multi::array<int, 1>::iterator cfirst = arr1d.begin();  // correctly fails to compile

		BOOST_TEST( arr1d[0] == 1 );
	}

	return boost::report_errors();
}
