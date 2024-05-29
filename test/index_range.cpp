// Copyright 2021-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array_ref.hpp>

#include <numeric>  // for accumulate

// Suppress warnings from boost
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
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4244)
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#include <boost/iterator/transform_iterator.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_range) {
#if defined(__cpp_deduction_guides) && __cpp_deduction_guides && ! defined(__NVCC__)
	BOOST_REQUIRE(( multi::range{5, 5}.empty() ));
#else
	BOOST_REQUIRE(( multi::range<std::ptrdiff_t>{5, 5}.empty() ));
#endif
	{
		auto drng = multi::range<std::ptrdiff_t>{5, 10};
		std::vector<double> vec(drng.begin(), drng.end());  // testing std::vector NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( vec[1] == 6 );
	}
	{
		auto drng = multi::range<std::ptrdiff_t>{5, 10};

		auto fun = [](auto idx) { return idx + 1; };

		std::vector<double> vec(  // testing std::vector NOLINT(fuchsia-default-arguments-calls)
			boost::make_transform_iterator(drng.begin(), fun),
			boost::make_transform_iterator(drng.end(), fun)
		);
		BOOST_REQUIRE( vec[1] == 7 );
	}
}

BOOST_AUTO_TEST_CASE(crazy_range) {
	// auto trng = multi::range(
	//    multi::detail::tuple<int, int>{5, 3},
	//    multi::detail::tuple<int, int>{5, 9},
	//  [](auto t , int d) {return std::get<1>(t) + d;}
	//  [](auto t1, auto t2) {return std::get<1>(t1) - std::get<1>(t2);}
	// );

	// BOOST_REQUIRE( trng[0] == (std::tuple{5, 3}) );
	// BOOST_REQUIRE( trng[1] == (std::tuple{5, 4}) );

	// BOOST_REQUIRE( *trng.begin() == (std::tuple{5, 3}) );
	// BOOST_REQUIRE( *(trng.begin() + 1) == (std::tuple{5, 4}) );
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr) {
	// BOOST_REQUIRE( multi::extension_t<int>{5} == 5 );
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> const irng{5, 12};

	BOOST_REQUIRE( irng.contains(6) );
	BOOST_REQUIRE( ! irng.contains(12) );

	BOOST_REQUIRE( * irng.begin()      ==  5 );
	BOOST_REQUIRE( *(irng.begin() + 1) ==  6 );

	BOOST_REQUIRE(   irng.first()      ==  5 );
	BOOST_REQUIRE(   irng.last()       == 12 );

	BOOST_REQUIRE(   irng.front()      ==  5 );
	BOOST_REQUIRE(   irng.back ()      == 11 );

	std::vector<int> vec = {5, 6, 7, 8, 9, 10, 11};  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)

	BOOST_REQUIRE(std::equal(irng.begin(), irng.end(), vec.begin(), vec.end()));  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)

	auto sum = std::accumulate(irng.begin(), irng.end(), 0);
	BOOST_REQUIRE( sum == 5 + 6 + 7 + 8 + 9 + 10 + 11 );
}

BOOST_AUTO_TEST_CASE(multi_range2) {
	multi::index_extension const iex(10);

	BOOST_REQUIRE( *begin(iex) == 0 );
	BOOST_REQUIRE( size(iex) == 10 );
	BOOST_REQUIRE( iex[0] == 0 );
	BOOST_REQUIRE( iex[1] == 1 );
	BOOST_REQUIRE( iex[9] == 9 );

	auto const xbeg = begin(iex);
	BOOST_REQUIRE( xbeg[0] == iex[0] );
	BOOST_REQUIRE( xbeg[1] == iex[1] );

	BOOST_REQUIRE( std::accumulate( begin(iex), end(iex), static_cast<multi::index_extension::value_type>(0U)) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

	{
		multi::iextensions<3> const ies({
			{0, 3},
			{0, 4},
			{0, 5},
		});

		auto const ies0 = std::get<0>(ies);
		auto const ies1 = std::get<1>(ies);
		auto const ies2 = std::get<2>(ies);

		BOOST_REQUIRE( ies0.size() == 3 );
		BOOST_REQUIRE( ies1.size() == 4 );
		BOOST_REQUIRE( ies2.size() == 5 );

		BOOST_REQUIRE( std::get<0>(ies).size() == 3 );
		BOOST_REQUIRE( std::get<1>(ies).size() == 4 );
		BOOST_REQUIRE( std::get<2>(ies).size() == 5 );

	#ifndef _MSC_VER  // doesn't work in MSVC 14.3 in c++17 mode
		auto [eyes, jays, kays] = ies;
		BOOST_REQUIRE( eyes.size() == 3 );
		BOOST_REQUIRE( jays.size() == 4 );
		BOOST_REQUIRE( kays.size() == 5 );
	#endif
	}
}
