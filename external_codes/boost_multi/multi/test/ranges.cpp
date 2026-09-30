// Copyright 2023-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, subarray, static_array  // IWYU pragma: keep

#include <boost/core/lightweight_test.hpp>

#include <algorithm>   // for std::ran  // IWYU pragma: keep  //  NOLINT(misc-include-cleaner)
#include <functional>  // for plus<>  // IWYU pragma: keep  // NOLINT(misc-include-cleaner)
#include <memory>      // for allocator  // IWYU pragma: keep  // NOLINT(misc-include-cleaner)

#if defined(__cplusplus) && (__cplusplus >= 202002L) && __has_include(<ranges>)
// clang 16 c++23 crashed when including ranges
#include <iterator>  // for random_access_iterator  // IWYU pragma: keep
#include <ranges>    // IWYU pragma: keep
#endif

#if defined(__cpp_lib_ranges_fold) && (__cpp_lib_ranges_fold >= 202207L)
#include <complex>      // for complex, real, operator==, imag  // IWYU pragma: keep
#include <numeric>      // for iota  // IWYU pragma: keep
#include <type_traits>  // for is_same_v  // IWYU pragma: keep
#include <utility>      // for pair  // IWYU pragma: keep
#endif

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		multi::array<int, 1> const arr({10}, 99);
		for(auto const& elem : arr) {  // NOLINT(altera-unroll-loops) test plain loop
			BOOST_TEST( elem == 99 );
		}
	}

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)
	// multi iterators must model the standard library ranges iterator/range concepts
	{
		static_assert(std::random_access_iterator<multi::array<int, 1>::iterator>);
		static_assert(std::random_access_iterator<multi::array<int, 1>::const_iterator>);
		static_assert(std::random_access_iterator<multi::array<int, 2>::iterator>);
		static_assert(std::random_access_iterator<multi::array<int, 2>::const_iterator>);
		static_assert(std::random_access_iterator<multi::array<int, 3>::iterator>);
		static_assert(std::random_access_iterator<multi::array<int, 3>::const_iterator>);

		static_assert(std::ranges::random_access_range<multi::array<int, 1>>);
		static_assert(std::ranges::random_access_range<multi::array<int, 2>>);
		static_assert(std::ranges::random_access_range<multi::array<int, 3>>);
		static_assert(std::ranges::random_access_range<multi::array<int, 2> const>);
	}
#endif

	// range accumulate
	{
#if defined(__cpp_lib_ranges_fold) && (__cpp_lib_ranges_fold >= 202207L)
		namespace multi = boost::multi;

		static constexpr auto accumulate = [](auto const& R) {
			return std::ranges::fold_left(R, 0, std::plus<>{});
		};

		auto const values = multi::array<int, 2>{
			{2, 0, 2, 2},
			{2, 2, 0, 4},
			{2, 2, 0, 4},
			{2, 2, 0, 0},
			{2, 7, 0, 2},
			{2, 2, 4, 4},
		};

		boost::multi::array<int, 1, std::allocator<int>> aaa = {1, 2, 3};

		constexpr auto rowOddSum = [](auto const& arr) {
			return std::ranges::find_if(arr, [](auto const& row) { return (accumulate(row) & 1) == 1; });
		};

		auto const result = rowOddSum(values);

		BOOST_TEST( result - values.begin() == 4 );

		{
			auto const lazy1D = [](auto ii) noexcept { return static_cast<float>(ii); } ^ multi::extents_t(6);
			auto const sum    = std::ranges::fold_left(lazy1D, float{}, std::plus<>{});
			BOOST_TEST( std::abs(sum - (0.0F + 1.0F + 2.0F + 3.0F + 4.0F + 5.0F)) < 1e-12F );
		}
#endif
	}

	// BOOST_AUTO_TEST_CASE(range_find)
	{
#if defined(__cpp_lib_ranges_fold) && (__cpp_lib_ranges_fold >= 202207L)
		using Array2D = multi::array<int, 2>;

		Array2D const a = {
			{1, 2},
			{3, 4},
		};
		{
			auto const needle = std::ranges::find_if(a, [](auto const& row) { return row[0] == 9; });
			BOOST_TEST(needle == a.end());
		}
		{
			std::ranges::equal_to eto;

			auto a1     = a[1];
			auto a1_val = +a[1];

			bool const res = eto(a1_val, a1);
			BOOST_TEST( res );
		}
		{
			auto&& a1 = a[1];

			auto const needle = std::ranges::find(a, a1);
			BOOST_TEST(needle != a.end());
			BOOST_TEST( *needle == a1 );
			BOOST_TEST( *needle == a[1] );
		}
		[&] {
			auto const needle = std::ranges::find(a, a[1]);
			BOOST_TEST(needle != a.end());
			BOOST_TEST( *needle == a[1] );
		}();
#endif
	}

	// #if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)
	// BOOST_AUTO_TEST_CASE(range_copy_n_1D) {
	//  namespace multi = boost::multi;

	//  multi::array<int, 1> const X1 = {1, 2, 3};
	//  multi::array<int, 1> X2(X1.extensions());

	//  std::ranges::copy_n(X1.begin(), 10, X2.begin());

	//  BOOST_TEST( X1 == X2 );
	// }

	// BOOST_AUTO_TEST_CASE(range_copy_n) {
	//  namespace multi = boost::multi;

	//  multi::array<int, 2> const X1({ 10, 10 }, 99);
	//  multi::array<int, 2> X2(X1.extensions());

	//  std::ranges::copy_n(X1.begin(), 10, X2.begin());
	//  BOOST_TEST( X1 == X2 );
	// }
	// #endif

#if defined(__cpp_lib_ranges_zip) && (__cpp_lib_ranges_zip >= 202110L) && !defined(_MSC_VER) && !defined(__NVCOMPILER)
	{
		multi::array<int, 2> A = {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		};

		multi::array<int, 1> const V = {10, 11, 12};

		multi::array<int, 1> const R{std::ranges::views::zip_transform(std::plus<>{}, A[0], V)};

		BOOST_TEST( R[0] == 11 );
		BOOST_TEST( R[1] == 13 );
		BOOST_TEST( R[2] == 15 );
	}
	{
		multi::array<int, 2> A = {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		};

		multi::array<int, 1> const V = {10, 11, 12};

		multi::array<int, 1> const R = std::ranges::views::zip_transform(std::plus<>{}, A[0], V);

		BOOST_TEST( R[0] == 11 );
		BOOST_TEST( R[1] == 13 );
		BOOST_TEST( R[2] == 15 );
	}
	{
		static_assert(std::ranges::viewable_range<boost::multi::const_subarray<int, 1, int*>>);
		multi::array<int, 2> const A = {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		};

		BOOST_TEST( A.begin() == std::ranges::begin(A) );
		BOOST_TEST( A.end() == std::ranges::end(A) );

		BOOST_TEST( A().begin() == std::ranges::begin(A()) );
		BOOST_TEST( A().end() == std::ranges::end(A()) );

		BOOST_TEST( A[1].begin() == std::ranges::begin(A[1]) );
		BOOST_TEST( A[1].end() == std::ranges::end(A[1]) );

		multi::array<int, 1> const V = {10, 11, 12};
		multi::array<int, 1> const R = std::ranges::views::zip_transform(std::plus<>{}, A[0], V);

		BOOST_TEST( R[0] == 11 );
		BOOST_TEST( R[1] == 13 );
		BOOST_TEST( R[2] == 15 );

		// R = std::ranges::views::zip_transform(std::plus<>{}, A[0], V);
	}
	{
		multi::array<int, 2> A = {
			{2, 2, 2, 2},
			{2, 0, 0, 2},
			{2, 0, 0, 2},
			{2, 2, 2, 2}
		};

		std::ranges::fill(A({1, 3}, {1, 3}).elements(), 1);

		BOOST_TEST( A[1][1] == 1 );
		BOOST_TEST( A[2][2] == 1 );
	}
#endif

	return boost::report_errors();
}
