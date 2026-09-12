// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	/* drop */
	{
		multi::array<int, 2> AA({
									{1, 5},
									{1, 5}
        },
								0);
		AA[1][1] = 12;
		AA[2][1] = 13;

		BOOST_TEST( AA.dropped(1)[1][1] == 13 );
	}

	return boost::report_errors();
}
