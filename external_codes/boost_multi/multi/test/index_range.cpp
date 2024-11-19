// Copyright 2021-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array_ref.hpp>

#include <algorithm>  // for equal
#include <numeric>    // for accumulate
#include <vector>     // for vector

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(multi_range_in_constexpr) {
	BOOST_TEST(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> const irng{ 5, 12 };

	BOOST_TEST(  irng.contains(6) );
	BOOST_TEST( !irng.contains(12) );

	BOOST_TEST( * irng.begin()      ==  5 );
	BOOST_TEST( *(irng.begin() + 1) ==  6 );

	BOOST_TEST(   irng.first()       ==  5 );
	BOOST_TEST(   irng.last()       == 12 );

	BOOST_TEST(   irng.front()      ==  5 );
	BOOST_TEST(   irng.back ()      == 11 );

	std::vector<int> vec = { 5, 6, 7, 8, 9, 10, 11 };  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)

	BOOST_TEST(std::equal(irng.begin(), irng.end(), vec.begin(), vec.end()));  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)

	auto sum = std::accumulate(irng.begin(), irng.end(), 0);
	BOOST_TEST( sum == 5 + 6 + 7 + 8 + 9 + 10 + 11 );
}

BOOST_AUTO_TEST_CASE(multi_range2) {
	multi::index_extension const iex(10);

	BOOST_TEST( *begin(iex) == 0 );
	BOOST_TEST( size(iex) == 10 );
	BOOST_TEST( iex[0] == 0 );
	BOOST_TEST( iex[1] == 1 );
	BOOST_TEST( iex[9] == 9 );

	auto const xbeg = begin(iex);
	BOOST_TEST( xbeg[0] == iex[0] );
	BOOST_TEST( xbeg[1] == iex[1] );

	BOOST_TEST( std::accumulate( begin(iex), end(iex), static_cast<multi::index_extension::value_type>(0U)) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

	{
		multi::iextensions<3> const ies({
			{0, 3},
			{0, 4},
			{0, 5},
		});

		using std::get;

		auto const ies0 = get<0>(ies);
		auto const ies1 = get<1>(ies);
		auto const ies2 = get<2>(ies);

		BOOST_TEST( ies0.size() == 3 );
		BOOST_TEST( ies1.size() == 4 );
		BOOST_TEST( ies2.size() == 5 );

		BOOST_TEST( get<0>(ies).size() == 3 );
		BOOST_TEST( get<1>(ies).size() == 4 );
		BOOST_TEST( get<2>(ies).size() == 5 );

#ifndef _MSC_VER  // doesn't work in MSVC 14.3 in c++17 mode
		auto const [eyes, jays, kays] = ies;
		BOOST_TEST( eyes.size() == 3 );
		BOOST_TEST( jays.size() == 4 );
		BOOST_TEST( kays.size() == 5 );
#endif
	}
}
return boost::report_errors();}
