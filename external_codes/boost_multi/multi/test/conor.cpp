// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// vvv this has no effect, needs to be passed directly from compilation line "-Wno-psabi"
// #ifdef __GNUC__
// #pragma GCC diagnostic ignored "-Wpsabi"  // for ranges backwards compatibility message
// #endif

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <functional>  // IWYU pragma: keep
#include <iostream>    // IWYU pragma: keep
#include <limits>      // IWYU pragma: keep

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep  // for __cpp_lib_ranges_fold
#endif

#if (__cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)) && defined(__cpp_lib_ranges_fold) && (__cpp_lib_ranges_fold >= 202207L)

#include <boost/multi/array.hpp>

#include <algorithm>  // IWYU pragma: keep  // for std::equal
#include <cmath>      // for std::abs
#include <concepts>   // for constructible_from  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep
#include <iterator>   // IWYU pragma: keep
#include <ranges>     // IWYU pragma: keep
#include <tuple>      // for std::get  // IWYU pragma: keep

namespace stdr = std::ranges;
namespace stdv = std::views;

auto printR2(auto const& lbl, auto const& arr2D) {
	// return fmt::print("{} = \n[{}]\n\n", lbl, fmt::join(arr2D, ",\n "));
	std::cout << lbl << " = \n";
	for(auto const& row : arr2D) {
		for(auto const& elem : row)
			std::cout << elem << ", ";
		std::cout << '\n';
	}
	std::cout << '\n';
}

constexpr auto maxR1 = []<class R, class V = stdr::range_value_t<R>>(R const& row, V low = std::numeric_limits<V>::lowest()) {
	return stdr::fold_left(row, low, stdr::max);
};

constexpr auto sumR1 = []<class R, class V = stdr::range_value_t<R>>(R const& rng, V zero = {}) {
	return stdr::fold_left(rng, zero, std::plus<>{});
};

#define FWD(var) std::forward<decltype(var)>(var)

auto softmax(auto&& matrix) noexcept {
	return           //
		FWD(matrix)  //
		|
		stdv::transform([](auto&& row) {
			auto max = maxR1(row);
			return        //
				FWD(row)  //
				|
				stdv::transform([=](auto ele) noexcept { return std::exp(ele - max); });
		})  //
		|
		stdv::transform([](auto&& nums) {
			auto den = sumR1(nums);
			return         //
				FWD(nums)  //
				|
				stdv::transform([=](auto num) noexcept { return num / den; });
		});
}

namespace multi = boost::multi;

auto main() -> int {
	auto const matrix =
		([](auto ii) noexcept { return static_cast<float>(ii); } ^
		 multi::extents_t(6))
			.partitioned(2);

	printR2("matrix", matrix);

	printR2("softmax", softmax(matrix));

	auto const alloc_matrix = multi::array<float, 2>{
		{0.0F, 1.0F, 2.0F},
		{3.0F, 4.0F, 5.0F}
	};

	printR2("softmax", softmax(alloc_matrix));

	auto const sofmax_copy = multi::array<float, 2>(softmax(alloc_matrix));

	BOOST_TEST( std::abs(sumR1(sofmax_copy[1]) - 1.0F) < 1e-12F );

	// auto softmax_copy = +softmax(alloc_matrix);

	return boost::report_errors();
}
#else
auto main() -> int {
	return boost::report_errors();
}
#endif
