// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#if __cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)
#include <boost/multi/array.hpp>
#include <boost/multi/elementwise.hpp>

#include <functional>  // for plus<>
#include <iostream>
#include <ranges>

namespace multi = boost::multi;

auto main() -> int {
	auto const A = multi::array<int, 2>{
		{0, 1, 2},
		{3, 4, 5}
	};
	auto const B = multi::array<int, 2>{
		{ 0, 10, 20},
		{30, 40, 50}
	};

	using multi::elementwise::apply;

	multi::array<int, 2> const C1 = apply(std::plus<>{}, A, B);

	BOOST_TEST( C1[1][1] == std::plus<>{}(A[1][1], B[1][1]) );  // NOLINT(build/include_what_you_use)

	using multi::elementwise::operator+;

	multi::array<int, 2> const C2 = A + B;

	BOOST_TEST( C1 == C2 );

	using multi::elementwise::operator-;

	multi::array<int, 2> const C3 = A - (-B);

	BOOST_TEST( C2 == C3 );

	{
		multi::array<double, 1> aa = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		BOOST_TEST( std::abs( exp(aa)[2] - std::exp(aa[2]) ) < 1e-12 );

		auto exp_aa_copy = +exp(aa);

		BOOST_TEST( std::abs( exp_aa_copy[2] - std::exp(aa[2]) ) < 1e-12 );
	}
	// auto softmax =
	// std::views::transform([](auto row) { return row | stdv::transform([max = maxR1(row)](auto ele) { return exp(ele - max); }); })
	// | std::views::transform([](auto nums) { return nums | stdv::transform([den = sumR1(nums)](auto num) { return num / den; }); });

	// 	{
	// 		using multi::elementwise::apply;

	// 		auto const matrix =
	// 			([](auto ii) noexcept { return static_cast<float>(ii); } ^
	// 			 multi::extents_t(6))
	// 				.partitioned(2);

	// 		constexpr auto maxR1 = []<class R, class V = std::ranges::range_value_t<R>>(R const& row, V low = std::numeric_limits<V>::lowest()) {
	// 			return std::ranges::fold_left(row, low, std::ranges::max);
	// 		};

	// 		#define BOOST_MULTI_FWD(var) std::forward<decltype(var)>(var)
	// 		using multi::elementwise::apply_front;
	// 		auto matrix_minus_row_max = apply_front([&](auto row) { return apply_front([max = maxR1(row)](auto elem) { return elem/*- max*/; }, row); }, matrix);

	// 		std::cout << std::abs( matrix_minus_row_max[0][2] ) << '\n';
	// 		BOOST_TEST( std::abs( matrix_minus_row_max[0][2] ) < 1e-12F );
	// 	}
	return boost::report_errors();
}
#else
auto main() -> int {
	return boost::report_errors();
}
#endif
