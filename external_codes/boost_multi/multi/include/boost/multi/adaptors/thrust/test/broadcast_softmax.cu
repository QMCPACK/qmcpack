// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
#include <boost/multi/array.hpp>  // from https://github.com/correaa/boost-multi
#include <boost/multi/elementwise.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include <algorithm>  // for max
#include <cmath>      // for exp, __cpp_lib_ranges
#include <iostream>
#include <limits>
#include <numeric>
#include <string>  // for basic_string, operator<<
#include <utility>
#include <ranges>

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

template<class R, class V = typename R::value_type>
BOOST_MULTI_HD constexpr auto maxR1(R const& rng) noexcept {  // NOLINT(readability-identifier-naming,misc-use-internal-linkage)
	printf("M");
	return std::accumulate(rng.begin(), rng.end(), std::numeric_limits<V>::lowest(), [](auto const& a, auto const& b) { return std::max(a, b); });
}

template<class R, class V = typename R::value_type>
BOOST_MULTI_HD constexpr auto sumR1(R const& rng, V zero = {}) noexcept {
	printf("S");
	return std::accumulate(rng.begin(), rng.end(), zero);
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

}  // namespace

namespace multi = boost::multi;

struct iden_t {
	BOOST_MULTI_HD constexpr auto operator()(multi::index idx) const -> float {
		return static_cast<float>(idx);
	}
};

constexpr iden_t iden;

int main() {
	auto const lazy_matrix =
		(iden ^ multi::extents_t(6))
			.partitioned(2);

	printR2("lazy matrix", lazy_matrix);

	// printR2("softmax of lazy array", softmax(lazy_matrix));
	printR2("softmax2 of lazy array", softmax2(lazy_matrix));

	multi::array<float, 2> alloc_matrix = {
		{0.0F, 1.0F, 2.0F},
		{3.0F, 4.0F, 5.0F}
	};

	// printR2("softmax of alloc array", softmax(alloc_matrix));
	printR2("softmax2 of alloc array", softmax2(alloc_matrix));

	// materialize
	multi::array<float, 2> const sofmax_copy(softmax2(alloc_matrix));
	printR2("materialized softmax", sofmax_copy);

	// materialize
	multi::array<float, 2, thrust::device_allocator<float> > const sofmax_copy_universal(softmax2(alloc_matrix));
	printR2("materialized universal softmax", sofmax_copy);

	return boost::report_errors();
}
#else
int main() {
		return boost::report_errors();
}
#endif
