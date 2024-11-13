// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

// IWYU pragma: no_include <type_traits>                      // for add_const_t
#include <utility>  // for as_const

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	/* subarray_assignment */
	{
		multi::array<int, 3> A1({3, 4, 5}, 99);
		A1[2][1][1] = 88;

		auto constA2 = std::as_const(A1)[2];
		BOOST_TEST( constA2[1][1] == 88 );

		auto A2 = A1[2];
		BOOST_TEST( A2[1][1] == 88 );

		A2[1][1] = 77;
		BOOST_TEST( A2[1][1] == 77 );
	}

	/* subarray_base */
	{
		multi::array<int, 3> A1({3, 4, 5}, 99);

		auto&& Asub  = A1();
		*Asub.base() = 88;

		BOOST_TEST( A1[0][0][0] == 88 );

		*A1().base() = 77;

		BOOST_TEST( A1[0][0][0] == 77 );

		// *std::as_const(Asub).base() = 66;  // should not compile, read-only
	}

	/* test ref(begin, end)*/
	{
		multi::array<int, 2> A2D = {
			{1, 2},
			{3, 4}
		};
		BOOST_TEST( A2D[0][0] == 1 );

		multi::const_subarray<int, 2> R2D(A2D.begin(), A2D.end());
		BOOST_TEST( R2D.addressof()== A2D.addressof() );
	}

	/* test ref(begin, end)*/
	// {
	//  multi::array<int, 2> A2D = {
	//      {1, 2},
	//         {3, 4}
	//  };
	//  BOOST_TEST( A2D[0][0] == 1 );

	//  multi::const_subarray<int, 2> R2D(A2D.home(), A2D.sizes());
	//  BOOST_TEST( R2D.addressof()== A2D.addressof() );
	// }

	/* test ref(begin, end)*/
	{
		multi::array<int, 2> A2D = {
			{1, 2},
			{3, 4}
		};
		BOOST_TEST( A2D[0][0] == 1 );

		multi::subarray<int, 2> R2D(A2D.begin(), A2D.end());
		BOOST_TEST( R2D.addressof()== A2D.addressof() );
		R2D[0][0] = 77;
	}

	// {
	//  multi::array<double, 2> A2D({10000, 10000}, 55.5);
	//  auto const&             A2D_block = A2D({1000, 9000}, {1000, 9000});

	//  multi::array<double, 2> B2D({10000, 10000}, 66.6);
	//  auto const&             B2D_block = B2D({1000, 9000}, {1000, 9000});

	//  *B2D_block.begin() = *A2D_block.begin();  // doesn't compile, CORRECT
	// }

	return boost::report_errors();
}
