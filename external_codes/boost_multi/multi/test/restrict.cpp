// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/restriction.hpp>  // for array, dynamic_array, num_elements

#include <boost/core/lightweight_test.hpp>

// IWYU pragma: no_include <utility>  // for forward, declval, move

namespace multi = boost::multi;

auto main() -> int {
	{
		auto const R1 = multi::restriction({2, 3}, [](auto ii, auto jj) { return ii + jj; });
		auto const R2 = [](auto ii, auto jj) { return ii + jj; } ^ multi::extents_t{2, 3};

		auto const R3 = multi::restricted([](auto ii, auto jj) { return ii + jj; }, multi::extents_t<2>{2, 3});
		auto const R4 = multi::restricted([](auto ii, auto jj) { return ii + jj; }, multi::extents_t{2, 3});
		auto const R5 = multi::restricted<2>([](auto ii, auto jj) { return ii + jj; }, {2, 3});

		// auto const R6 = multi::restrict([](auto ii, auto jj) { return ii + jj; }, {2, 3});

		BOOST_TEST( R1[1][1] == 2 );
		BOOST_TEST( R2[1][1] == 2 );
		BOOST_TEST( R3[1][1] == 2 );
		BOOST_TEST( R4[1][1] == 2 );
		BOOST_TEST( R5[1][1] == 2 );
		// BOOST_TEST( R6[1][1] == 2 );
	}

	return boost::report_errors();
}
