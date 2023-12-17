// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi views constness"  // test title NOLINT(cppcoreguidelines-macro-usage)
#include <boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<iostream>

namespace multi = boost::multi;

template<class Array1D>
void print(Array1D const& coll) {
	// *coll.begin() = 99;  // doesn't compile "assignment of read-only location"

	std::copy(coll.begin(), coll.end(), std::ostream_iterator<typename Array1D::value_type>{std::cout, ", "});
	std::cout<<std::endl;
}

BOOST_AUTO_TEST_CASE(const_views) {
	multi::array<int, 1> coll1 = {0, 8, 15, 47, 11, 42};
	print(coll1);  // prints "0, 8, 15, 47, 11, 42"

	print(coll1({0, 3}));  // similar to coll1 | take(3) // prints "0, 8, 15"

	auto&& coll1_take3 = coll1({0, 3});
	print(coll1_take3);  // prints "0, 8, 15"
}

template<class Array1D>
void fill_99(Array1D&& coll) {
	std::fill(coll.begin(), coll.end(), 99);
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

template<class Array2D>
void print_2d(Array2D const& coll) {
	// *(coll.begin()->begin()) = 99;  // doesn't compile "assignment of read-only location"

	std::for_each(coll.begin(), coll.end(), [](auto const& row) {
		std::copy(row.begin(), row.end(), std::ostream_iterator<typename Array2D::element_type>(std::cout, ", "));
		std::cout<<std::endl;
	});
}

BOOST_AUTO_TEST_CASE(const_views_2d) {
	multi::array<int, 2> coll1 = {
		{0, 8, 15, 47, 11, 42},
		{0, 8, 15, 47, 11, 42}
	};

	print_2d(coll1);  // prints "0, 8, 15, 47, 11, 42"

	print_2d(coll1({0, 2}, {0, 3}));  // similar to coll1 | take(3) // prints "0, 8, 15"

	auto&& coll1_take3 = coll1({0, 2}, {0, 3});
	print_2d(coll1_take3);  // prints "0, 8, 15"
}

template<class Array1D>
void fill_2d_99(Array1D&& coll) {
	// for(auto const& row : coll) {  // does not work because it would make it const
	std::for_each(coll.begin(), coll.end(), [](auto&& row) {
		std::fill(row.begin(), row.end(), 99);
	});
	// std::transform(coll.begin(), coll.end(), coll.begin(), [](auto&& row) {
	//  std::fill(row.begin(), row.end(), 99);
	//  return std::forward<decltype(row)>(row);
	// });
}

BOOST_AUTO_TEST_CASE(mutating_views_2d) {
	multi::array<int, 2> coll1 = {
		{0, 8, 15, 47, 11, 42},
		{0, 8, 15, 47, 11, 42}
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
