// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/elementwise.hpp>
#include <boost/multi/io.hpp>
// IWYU pragma: no_include "boost/multi/restriction.hpp"  // for restriction, operator!=
#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <iostream>
#include <numeric>
#include <utility>  // for move

namespace multi = boost::multi;

// NOLINTBEGIN(readability-identifier-length)
auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		multi::array<int, 2> A({3, 4});
		std::iota(A.elements().begin(), A.elements().end(), 0);

		BOOST_TEST((
			A == multi::array<int, 2>{
				{0, 1, 2, 3},
				{4, 5, 6, 7},
				{8, 9, 10, 11}
			}
		));

		multi::array<int, 2> B({3, 1});
		std::iota(B.elements().begin(), B.elements().end(), 100);

		BOOST_TEST((
			B == multi::array<int, 2>{
				{100},
				{101},
				{102}
			}
		));

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement ;

		auto const& C1 = ~(~A + (~B)[0]);

		std::cout << "C1 = " << C1 << '\n';

		multi::array<int, 1> b({3}, 0);
		std::iota(b.elements().begin(), b.elements().end(), 100);

		BOOST_TEST((
			b == multi::array<int, 1>{100, 101, 102}
		));

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement ;

		auto const& C2 = ~(~A + b);

		std::cout << "C2 = " << C2 << '\n';
	}
	{
		multi::array<int, 1> A({4}, int{});
		BOOST_TEST( A.size() == 4 );

		auto rep2 = std::move(A).repeated(2);
		BOOST_TEST( rep2.num_elements() == 8 );

		auto [rep2n, rep2m] = rep2.sizes();
		BOOST_TEST( rep2n == 2 );
		BOOST_TEST( rep2m == 4 );
	}

	return boost::report_errors();
}
// NOLINTEND(readability-identifier-length)
