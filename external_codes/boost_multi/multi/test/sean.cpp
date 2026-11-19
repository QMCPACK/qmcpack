// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include "boost/multi/restriction.hpp"  // for restriction, operator==

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <algorithm>  // IWYU pragma: keep  // for std::equal
#include <iterator>   // IWYU pragma: keep
#include <tuple>      // IWYU pragma: keep  // for std::tuple  // NOLINT(misc-include-cleaner)
// IWYU pragma: no_include <utility>  // for declval, forward, move

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#include <concepts>     // for constructible_from, defau...
#include <ranges>       // IWYU pragma: keep
#include <type_traits>  // for is_constructible_v
#endif

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep  // for _GLIBCXX_RELEASE, __GLIBC...
#endif

namespace multi = boost::multi;

// Observed failure: clang-15 (as its own toolchain, not later clang versions) paired with
// GCC-11's libstdc++ (_GLIBCXX_RELEASE 11) fails to instantiate std::ranges::ref_view<T> for
// *any* range T (confirmed even for std::vector/std::list, not specific to boost::multi) when
// piped through std::views::reverse and similar adaptors. See test/broadcast_softmax.cpp for
// the same guard.
#if defined(__clang__) && (__clang_major__ == 15) && defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE < 12)
#define BOOST_MULTI_STDRANGES_PIPE_BROKEN 1
#endif

auto main() -> int {  // NOLINT(bugprone-exception-escape,readability-function-cognitive-complexity)
#ifdef __NVCC__
	auto fun = [](auto ii, auto jj) noexcept { return static_cast<int>((10 * ii) + jj); };
	auto rst = fun ^ multi::extents_t(5, 5);
#else
	auto rst = [](auto ii, auto jj) noexcept { return static_cast<int>((10 * ii) + jj); } ^ multi::extents_t(5, 5);
#endif
	multi::array<int, 2> const AA = rst;

	BOOST_TEST( AA.size() == rst.size() );
	BOOST_TEST( AA.extents() == rst.extents() );

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
#ifndef BOOST_MULTI_STDRANGES_PIPE_BROKEN
	multi::array<int, 2> const BB(rst | std::ranges::views::reverse);

	BOOST_TEST( AA[0] == BB[4] );  // as A[0][0] == B[4][0] && A[0][1] == B[4][1] ...
	BOOST_TEST( AA[1] == BB[3] );  // as A[1][0] == B[3][0] && A[1][1] == B[3][1] ...
	BOOST_TEST( AA[2] == BB[2] );  // ...
	BOOST_TEST( AA[3] == BB[1] );
	BOOST_TEST( AA[4] == BB[0] );
#endif
#endif

	auto rstT = rst.transposed();

	using std::get;
	BOOST_TEST( get<0>(rstT.extents()) == get<1>(rst.extents()) );
	BOOST_TEST( get<1>(rstT.extents()) == get<0>(rst.extents()) );

	BOOST_TEST( rstT[1][2] == rst[2][1] );

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
	static_assert(std::weakly_incrementable<decltype(rstT.begin())>);
	static_assert(std::input_or_output_iterator<decltype(rstT.begin())>);
	BOOST_TEST( rstT.begin() == std::ranges::begin(rstT) );

	static_assert(std::constructible_from<decltype(rstT.end())>);
	static_assert(std::default_initializable<decltype(rstT.end())>);
	static_assert(std::is_constructible_v<decltype(rstT.end())>);
	static_assert(std::semiregular<decltype(rstT.end())>);
	BOOST_TEST( rstT.end() == std::ranges::end(rstT) );

	// static_assert(std::ranges::viewable_range<decltype(rstT)>);
	// auto rstTR = rstT | std::ranges::views::reverse;

	// BOOST_TEST( rstTR.back()[0] == rstT.front()[0] );
	// BOOST_TEST( rstTR.front()[0] == rstT.back()[0] );

	// auto rstTR2 = rstT.reversed();

	// BOOST_TEST( rstTR2[3][4] == rstTR[3][4] );
#endif

	return boost::report_errors();
}
