// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/elementwise.hpp>
#include <boost/multi/restriction.hpp>  // for restriction, operator!=

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <numeric>
// IWYU pragma: no_include <utility>  // for forward

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	multi::array<int, 1> arr(13);
	std::iota(arr.begin(), arr.end(), 0);

	multi::array<int, 1> brr(11);
	std::iota(brr.begin(), brr.end(), 0);

	// this method uses a restriction directly
	auto const& M1 = multi::restricted(
		[&arr, &brr](auto ii, auto jj) { return arr[ii] * brr[jj]; },
		multi::extents_t<2>(arr.extent(), brr.extent())
	);

	BOOST_TEST( M1[5][7] == 5*7 );

#if (!defined(__GNUC__) || (__GNUC__ > 10)) || defined(__clang__)
#if !defined(__CUDACC_VER_MAJOR__) || (__CUDACC_VER_MAJOR__ > 11)
	// this method brings both arrays to the right (final) dimensionality and then does elementwise mult
	using multi::elementwise::operator*;  // cppcheck-suppress constStatement

	auto const& M2 = (~arr.repeated(brr.size())) * brr.repeated(arr.size());

	BOOST_TEST( M2[5][7] == 5*7 );

	// it is actually only necessary to bring one to the final dimension
	using multi::elementwise::operator*;  // cppcheck-suppress constStatement

	auto const& M3 = ~arr.repeated(brr.size()) * brr;

	BOOST_TEST( M3[5][7] == 5*7 );

	// if one insists in using row and column 2D arrays (assignments are for brevity)
	multi::array<int, 2> Aarr({arr.size(), 1});
	(~Aarr)[0] = arr;
	multi::array<int, 2> Barr({1, brr.size()});
	Barr[0] = brr;

	using multi::elementwise::operator*;  // cppcheck-suppress constStatement

	auto const& M4 = (~((~Aarr)[0].repeated(Barr[0].size()))) * Barr[0];

	BOOST_TEST( M4[5][7] == 5*7 );
#endif
#endif

	return boost::report_errors();
}
