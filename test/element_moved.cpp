// Copyright 2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 10.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for array, apply, array_types<>::ele...

// IWYU pragma: no_include <type_traits>  // for remove_reference<>::type
// IWYU pragma: no_include <map>
// IWYU pragma: no_include <set>
// IWYU pragma: no_include <stack>
// IWYU pragma: no_include <algorithm>  // for fill_n
// IWYU pragma: no_include <cstdlib>
// IWYU pragma: no_include <pstl/glue_algorithm_defs.h>
#include <utility>  // for move, swap
#include <vector>   // for vector, operator==, vector<>::va...

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// explicit_move_subarray_vector_2d_assign
	{
		multi::array<std::vector<double>, 2> arrA({10, 10}, std::vector<double>(5, {}, {}));

		BOOST_TEST( arrA[2][2].size() == 5 );
		{
			using std::move;
			multi::array<std::vector<double>, 2> arrB(arrA().element_moved());  // (arrA.extensions());
			// arrB = arrA().element_moved();

			BOOST_TEST( arrA.size() == 10 );
			BOOST_TEST( arrB.size() == 10 );
			BOOST_TEST( arrA[2][2].empty() );
			BOOST_TEST( arrB[2][2].size() == 5 );
		}
	}

	// explicit_move_subarray_vector_2d_ctor
	// {
	//  multi::array<std::vector<double>, 2> arrA({10, 10}, std::vector<double>(5));

	//  BOOST_TEST( arrA[2][2].size() == 5 );
	//  {
	//      using std::move;
	//      multi::array<std::vector<double>, 2> arrB{arrA().element_moved()};

	//      BOOST_TEST( arrA.size() == 10 );
	//      BOOST_TEST( arrB.size() == 10 );
	//      BOOST_TEST( arrA[2][2].size() == 0 );
	//      BOOST_TEST( arrB[2][2].size() == 5 );
	//  }
	// }

	// BOOST_AUTO_TEST_CASE(explicit_move_subarray_vector_1d)
	// {
	//  multi::array<std::vector<double>, 1> arrA(10, std::vector<double>(5));

	//  BOOST_TEST( arrA[2].size() == 5 );
	//  {
	//      using std::move;
	//      multi::array<std::vector<double>, 1> arrB(arrA.extensions());
	//      arrB() = arrA().element_moved();
	//      // std::copy(arrA().element_moved().begin(), arrA().element_moved().end(), arrB.begin());

	//      BOOST_TEST( arrA.size() == 10 );
	//      BOOST_TEST( arrB.size() == 10 );
	//      BOOST_TEST( arrA[2].size() == 0 );
	//      BOOST_TEST( arrB[2].size() == 5 );
	//  }
	// }



	return boost::report_errors();
}
