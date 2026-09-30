// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && !defined(__clang__)  // for gcc 14
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for equal
#include <array>      // for array
// #include <iterator>   // for size, begin, end

namespace multi = boost::multi;

// NOLINTNEXTLINE(readability-function-cognitive-complexity,bugprone-exception-escape)
auto main() -> int {
	//  one_based_1D
	{
		// clang-format off
		multi::array<double, 1> const Ac({{0, 10}}, 0.0);
		// clang-format on

		BOOST_TEST( Ac.size() == 10 );

		//  multi::array<double, 1> Af({{1, 1 + 10}}, 0.);
		//  Af[1] = 1.;
		//  Af[2] = 2.;
		//  Af[3] = 3.;

		//  BOOST_TEST( Af[1] = 1. );
		//  BOOST_TEST( *Af.data_elements() == 1. );
		//  BOOST_TEST( size(Af) == 10 );
		//  BOOST_TEST( extension(Af).start() == 1 );
		//  BOOST_TEST( extension(Af).finish() == 11 );

		//  auto Af1 = multi::array<double, 1>(multi::extents_t<1>{multi::iextension{10}}, 0.).reindex(1);

		//  BOOST_TEST( size(Af1) == 10 );
		//  BOOST_TEST( Af1[10] == 0. );

		//  multi::array<double, 1> B({{0, 10}}, 0.);
		//  B[0] = 1.;
		//  B[1] = 2.;
		//  B[2] = 3.;

		//  BOOST_TEST( size(B) == 10 );
		//  BOOST_TEST( B != Af );
		//  BOOST_TEST( std::equal(begin(Af), end(Af), begin(B), end(B) ) );

		//  BOOST_TEST( Af.reindexed(0) == B );
	}

	// BOOST_AUTO_TEST_CASE(reindex_from_0_to_indexed_from_1)
	{
		multi::array<int, 2> A1(
			{
				{1, 1 + 10},
				{1, 1 + 20}
        },
			0
		);

		BOOST_TEST( A1.extent().front() ==  1);
		BOOST_TEST( A1.extent().back()  == 10);

		A1[1][1]   = 10;
		A1[2][2]   = 20;
		A1[3][3]   = 30;
		A1[10][20] = 990;

		auto const&& A0 = A1.reindexed(0, 0);

		BOOST_TEST( A0.extent().front() == 0);
		BOOST_TEST( A0.extent().back()  == 9);

		BOOST_TEST( A0[0][0]   == 10 );
		BOOST_TEST( A0[1][1]   == 20 );
		BOOST_TEST( A0[2][2]   == 30 );
		BOOST_TEST( A0[9][19]  == 990);
	}

	// BOOST_AUTO_TEST_CASE(one_based_2D)
	{
		multi::array<int, 2> const Ac({
										  {0, 10},
										  {0, 20}
        },
									  0);
		BOOST_TEST( Ac.size() == 10 );

		multi::array<int, 2> Af({
									{1, 1 + 10},
									{1, 1 + 20}
        },
								0);

		Af[1][1]   = 10;
		Af[2][2]   = 20;
		Af[3][3]   = 30;
		Af[10][20] = 990;

		BOOST_TEST( Af[1][1] == 10 );
		BOOST_TEST( Af[10][20] == 990 );  // cppcheck-suppress knownConditionTrueFalse ; for test
		BOOST_TEST( *Af.elements().begin() == 10 );
		BOOST_TEST( Af.elements()[Af.num_elements()-1] == 990 );
		BOOST_TEST( size(Af) == 10 );

		BOOST_TEST( Af.extent().first()  ==  1 );
		BOOST_TEST( Af.extent().last() == 11 );

		multi::array<int, 2> BB({
									{0, 10},
									{0, 20}
        },
								0);
		BB[0][0]  = 10;
		BB[1][1]  = 20;
		BB[2][2]  = 30;
		BB[9][19] = 990;

		BOOST_TEST( BB.size() == 10 );
		BOOST_TEST( BB != Af );
		BOOST_TEST( Af.reindexed(0, 0) == BB );

		BOOST_TEST( Af.reindexed(0, 0)[0][0] == BB[0][0] );

		BOOST_TEST( std::equal(begin(Af.reindexed(0, 0)), end(Af.reindexed(0, 0)), BB.begin(), BB.end()) );
		//  BOOST_TEST( std::equal(begin(Af), end(Af), begin(B.reindexed(1, 1)), end(B.reindexed(1, 1)) ) );
		//  BOOST_TEST( std::equal(begin(Af), end(Af), begin(B.reindexed(0, 1)), end(B.reindexed(0, 1)) ) );

		//  BOOST_TEST( Af.reindexed(0, 0) == B );

		//  B = Af;  // TODO(correaa) implement assignment for 1-based arrays
		//  BOOST_TEST( B[1][1] = 1. );
		//  BOOST_TEST( B[10][20] == 99.0 );
		//  BOOST_TEST( B == Af );
	}

	// BOOST_AUTO_TEST_CASE(one_base_2D_ref)
	{
		// clang-format off
	std::array<std::array<int, 5>, 3> arr = {{
		{{  10,  20,  30,  40,  50 }},
		{{  60,  70,  80,  90, 100 }},
		{{ 110, 120, 130, 140, 150 }},
	}};
		// clang-format on

		BOOST_TEST( arr[0][0] == 10 );

		multi::array_ref<int, 2> const& Ar = *&multi::array_ref<int, 2>({3, 5}, arr[0].data());
		BOOST_TEST( &Ar[1][3] == &arr[1][3] );

		// clang-format off
		multi::array_ref<int, 2> const& Ar2 = *&multi::array_ref<int, 2>(
			{
				{1, 1 + 3},
				{1, 1 + 5},
			},
			arr[0].data()
		);
		// clang-format on

		BOOST_TEST( sizes(Ar) == sizes(Ar2) );
		BOOST_TEST( &Ar2[1][1] == arr[0].data() );
		BOOST_TEST( &Ar2[2][4] == &arr[1][3] );

		BOOST_TEST( Ar2.extents() != Ar.extents() );
		BOOST_TEST( !(Ar2 == Ar) );
		BOOST_TEST( Ar2 != Ar );
	}

	{
		multi::array<int, 2> AA = {
			{'a', 'b'},
			{'c', 'd'},
		};

		BOOST_TEST(AA[0][0] == 'a' );
		BOOST_TEST(AA[1][1] == 'd' );

		auto&& Aone = AA.reindexed(1);

		BOOST_TEST(Aone[1][0] == 'a' );
		BOOST_TEST(Aone[2][1] == 'd' );
	}
	{
		multi::array<int, 2> AA = {
			{'a', 'b'},
			{'c', 'd'},
		};

		BOOST_TEST(AA[0][0] == 'a' );
		BOOST_TEST(AA[1][1] == 'd' );

		auto&& Aone = AA.reindexed(1, 1);

		BOOST_TEST(Aone[1][1] == 'a' );
		BOOST_TEST(Aone[2][2] == 'd' );
	}

	return boost::report_errors();
}
