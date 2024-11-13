// Copyright 2023-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>

#include <algorithm>  // for std::ranges::fold_left
#include <cmath>      // IWYU pragma: keep  for std::abs
// IWYU pragma: no_include <cstdlib>                          // for abs
// IWYU pragma: no_include <stdlib.h>                          // for abs

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(broadcast_as_fill) {
	multi::array<int, 1> bb = { 10, 11 };

	multi::array<int, 2> BB({ 10, 2 });

	// std::fill  (BB.begin(), BB.end(), bb);                                       // canonical way
	std::fill_n(BB.begin(), BB.size(), bb);  // canonical way

	// std::copy_n(bb.broadcasted().begin(), BB.size(), BB.begin());                // doesn't work because faulty implementation of copy_n
	// thrust::copy_n(bb.broadcasted().begin(), BB.size(), BB.begin());                // equivalent, using broadcast

	// std::copy_n(bb.broadcasted().begin(), bb.broadcasted().size(), BB.begin());  // incorrect, undefined behavior, no useful size()
	// std::copy  (bb.broadcasted().begin(), bb.broadcasted().end(), BB.begin());   // incorrect, undefined behavior, non-terminating loop (end is not reacheable)
	// BB = bb.broadcasted();

	BOOST_TEST( BB[0] == bb );
	BOOST_TEST( BB[1] == bb );

	BOOST_TEST( std::all_of(BB.begin(), BB.end(), [&bb](auto const& row) { return row == bb; }) );

	multi::array<double, 0> const one{1.0};

	BOOST_TEST( one == 1.0 );

	auto const& ones = one.broadcasted();
	BOOST_TEST( std::abs( *ones.begin() - 1.0 ) < 1.0e-8 );
}

return boost::report_errors();}
