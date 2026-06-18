// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #include <fmt/ranges.h>

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#if !defined(__clang_major__) || (__clang_major__ != 16)
#include <ranges>  // IWYU pragma: keep
#endif
#endif

#if defined(__cpp_lib_ranges) && !defined(_MSC_VER) && !defined(__circle_build__)

#include <boost/multi/array.hpp>  // from https://github.com/correaa/boost-multi
#include <boost/multi/elementwise.hpp>
#include <boost/multi/restriction.hpp>  // for operator^, restriction

#include <algorithm>   // for max
#include <cmath>       // for exp, __cpp_lib_ranges
#include <functional>  // IWYU pragma: keep  // for plus
#include <iostream>
#include <limits>
#include <numeric>
#include <string>  // for basic_string, operator<<
#include <utility>

namespace stdr = std::ranges;
namespace stdv = std::views;

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace {
void printR2(std::string const& lbl, auto const& arr2D) {  // NOLINT(readability-identifier-naming)
	//  fmt::print("\n{} = \n[{}]\n\n", lbl, fmt::join(arr2D, ",\n "));
	std::cout << lbl << "=\n";
	for(auto const& row : arr2D) {     // NOLINT(altera-unroll-loops)
		for(auto const& elem : row) {  // NOLINT(altera-unroll-loops)
			std::cout << elem << ' ';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}
}  // namespace

template<class R, class V = stdr::range_value_t<R>>
constexpr auto maxR1(R const& rng) noexcept {  // NOLINT(readability-identifier-naming,misc-use-internal-linkage)
	// fmt::print("M");
	std::cout << 'M';
#ifdef __cpp_lib_ranges_fold
	return stdr::fold_left(rng, std::numeric_limits<V>::lowest(), stdr::max);
#else
	return std::accumulate(rng.begin(), rng.end(), std::numeric_limits<V>::lowest(), stdr::max);
#endif
}

constexpr auto sumR1 = []<class R, class V = stdr::range_value_t<R>>(R const& rng, V zero = {}) noexcept {  // NOLINT(fuchsia-default-arguments-declarations)
	// fmt::print("S");
	std::cout << 'S';
#ifdef __cpp_lib_ranges_fold
	return stdr::fold_left(rng, zero, std::plus<>{});
#else
	return std::accumulate(rng.begin(), rng.end(), zero);
#endif
};

#define FWD(var) std::forward<decltype(var)>(var)

namespace multi = boost::multi;

namespace {

template<class M>
class ret_t {
	M mat_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	template<class MM>
	BOOST_MULTI_HD constexpr explicit ret_t(MM&& mat) : mat_{std::forward<MM>(mat)} {}  // NOLINT(bugprone-forwarding-reference-overload)

	BOOST_MULTI_HD constexpr auto operator()(multi::index irow) const {
		using multi::elementwise::operator-;
		using multi::elementwise::exp;

		auto mati = mat_[irow];
		return exp(std::move(mati) - maxR1(mati));
	}
};

auto softmax2(auto&& mat) noexcept {  // -> decltype(auto) {
	using multi::elementwise::operator-;
	using multi::elementwise::exp;
	using multi::elementwise::operator/;

	auto ret = [mat = FWD(mat)](multi::index irow) { auto mati = mat[irow]; return exp(std::move(mati) - maxR1(mati)); } ^ multi::extents_t<1>{2};

	// auto ret = ret_t<decltype(mat)>{FWD(mat)} ^ multi::extents_t<1>{2};

	return
		[ret = std::move(ret)](auto irow) {
			auto reti = ret[irow];
			return std::move(reti) / sumR1(reti);
		} ^
		multi::extents_t<1>{2};
}

auto softmax(auto&& matrix) noexcept {
	return FWD(matrix)  //
		 |              //
		   stdv::transform([](auto&& row) {
			   auto max = maxR1(row);
			   return FWD(row) |
					  stdv::transform([=](auto ele) noexcept { return std::exp(ele - max); });
		   })  //
		 |     //
		   stdv::transform([](auto&& nums) {
			   auto den = sumR1(nums);
			   return FWD(nums) |
					  stdv::transform([=](auto elem) noexcept { return elem / den; });
		   });
}
}  // namespace

namespace multi = boost::multi;

struct iden_t {
	BOOST_MULTI_HD constexpr auto operator()(multi::index idx) const -> float {
		return static_cast<float>(idx);
	}
};

constexpr iden_t iden;

auto main() -> int {
	auto const lazy_matrix =
		(iden ^ multi::extents_t(6))
			.partitioned(2);

	printR2("lazy matrix", lazy_matrix);

	printR2("softmax of lazy array", softmax(lazy_matrix));
	printR2("softmax2 of lazy array", softmax2(lazy_matrix));

	multi::array<float, 2> alloc_matrix = {
		{0.0F, 1.0F, 2.0F},
		{3.0F, 4.0F, 5.0F}
	};

	printR2("softmax of alloc array", softmax(alloc_matrix));
	printR2("softmax2 of alloc array", softmax2(alloc_matrix));

	// materialize
	multi::array<float, 2> const sofmax_copy(softmax(alloc_matrix));

	// printR2("materialized softmax", sofmax_copy);

	//    assert(std::abs(sumR1(sofmax_copy[1]) - 1.0F) < 1e-12F);

	return boost::report_errors();
}
#else
auto main() -> int { return boost::report_errors(); }
#endif
