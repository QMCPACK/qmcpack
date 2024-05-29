// Copyright 2019-2023 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>

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

BOOST_AUTO_TEST_CASE(one_based_1D) {
	multi::array<double, 1> const Ac({{0, 10}}, 0.0);
	BOOST_REQUIRE( Ac.size() == 10 );

//  multi::array<double, 1> Af({{1, 1 + 10}}, 0.);
//  Af[1] = 1.;
//  Af[2] = 2.;
//  Af[3] = 3.;

//  BOOST_REQUIRE( Af[1] = 1. );
//  BOOST_REQUIRE( *Af.data_elements() == 1. );
//  BOOST_REQUIRE( size(Af) == 10 );
//  BOOST_REQUIRE( extension(Af).start() == 1 );
//  BOOST_REQUIRE( extension(Af).finish() == 11 );

//  auto Af1 = multi::array<double, 1>(multi::extensions_t<1>{multi::iextension{10}}, 0.).reindex(1);

//  BOOST_REQUIRE( size(Af1) == 10 );
//  BOOST_REQUIRE( Af1[10] == 0. );

//  multi::array<double, 1> B({{0, 10}}, 0.);
//  B[0] = 1.;
//  B[1] = 2.;
//  B[2] = 3.;

//  BOOST_REQUIRE( size(B) == 10 );
//  BOOST_REQUIRE( B != Af );
//  BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B), end(B) ) );

//  BOOST_REQUIRE( Af.reindexed(0) == B );
}

BOOST_AUTO_TEST_CASE(one_based_2D) {
	multi::array<double, 2> const Ac({{0, 10}, {0, 20}}, 0.0);
	BOOST_REQUIRE( Ac.size() == 10 );

//  multi::array<double, 2> Af({{1, 1 + 10}, {1, 1 + 20}}, 0.);
//  Af[1][1] = 1.;
//  Af[2][2] = 2.;
//  Af[3][3] = 3.;
//  Af[10][20] = 99.;

//  BOOST_REQUIRE( Af[1][1] = 1. );
//  BOOST_REQUIRE( Af[10][20] == 99. );
//  BOOST_REQUIRE( *Af.data_elements() == 1. );
//  BOOST_REQUIRE( Af.data_elements()[Af.num_elements()-1] == 99. );
//  BOOST_REQUIRE( size(Af) == 10 );
//  BOOST_REQUIRE( extension(Af).start()  ==  1 );
//  BOOST_REQUIRE( extension(Af).finish() == 11 );

//  auto Af1 = multi::array<double, 2>({10, 10}, 0.).reindex(1, 1);

//  BOOST_REQUIRE( size(Af1) == 10 );
//  BOOST_REQUIRE( Af1[10][10] == 0. );

//  multi::array<double, 2> B({{0, 10}, {0, 20}}, 0.);
//  B[0][0] = 1.;
//  B[1][1] = 2.;
//  B[2][2] = 3.;
//  B[9][19] = 99.;

//  BOOST_REQUIRE( size(B) == 10 );
//  BOOST_REQUIRE( B != Af );
//  BOOST_REQUIRE( std::equal(begin(Af.reindexed(0, 0)), end(Af.reindexed(0, 0)), begin(B), end(B)) );
//  BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B.reindexed(1, 1)), end(B.reindexed(1, 1)) ) );
//  BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B.reindexed(0, 1)), end(B.reindexed(0, 1)) ) );

//  BOOST_REQUIRE( Af.reindexed(0, 0) == B );

//  B = Af;  // TODO(correaa) implement assignment for 1-based arrays
//  BOOST_REQUIRE( B[1][1] = 1. );
//  BOOST_REQUIRE( B[10][20] == 99. );
//  BOOST_REQUIRE( B == Af );
}

BOOST_AUTO_TEST_CASE(one_base_2D_ref) {
	std::array<std::array<double, 5>, 3> arr = {{
		{{ 1.0,  2.0,  3.0,  4.0,  5.0}},
		{{ 6.0,  7.0,  8.0,  9.0, 10.0}},
		{{11.0, 12.0, 13.0, 14.0, 15.0}},
	}};

	BOOST_REQUIRE( arr[0][0] == 1.0 );

//  multi::array_ref<double, 2> const& Ar = *multi::array_ptr<double, 2>(&arr[0][0], {3, 5});
//  BOOST_REQUIRE( &Ar[1][3] == &arr[1][3] );

//  multi::array_ref<double, 2> const& Ar2 = *multi::array_ptr<double, 2>(&arr[0][0], {{1, 1+3}, {1, 1+5}});
//  BOOST_REQUIRE( sizes(Ar) == sizes(Ar2) );
//  BOOST_REQUIRE( &Ar2[1][1] == &arr[0][0] );
//  BOOST_REQUIRE( &Ar2[2][4] == &arr[1][3] );

//  BOOST_REQUIRE( Ar2.extensions() != Ar.extensions() );
//  BOOST_REQUIRE( not(Ar2 == Ar) );
//  BOOST_REQUIRE( Ar2 != Ar );
//  BOOST_REQUIRE( extensions(Ar2.reindexed(0, 0)) == extensions(Ar) );
//  BOOST_REQUIRE( Ar2.reindexed(0, 0) == Ar );

//  static_assert( not std::is_assignable<decltype(Ar2.reindexed(0, 0)[0][0]), double>{}, "!" );
}
