// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep  // for __cpp_lib_ranges_fold
#endif

#if (__cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)) && defined(__cpp_lib_ranges_fold) && (__cpp_lib_ranges_fold >= 202207L)
#include <boost/multi/array.hpp>
#include <boost/multi/elementwise.hpp>

#include <algorithm>   // IWYU pragma: keep  // for std::equal
#include <cmath>       // for std::abs
#include <concepts>    // for constructible_from  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep
#include <functional>  // for std::plus  // IWYU pragma: keep
#include <iostream>    // for std::cout  // NOLINT(misc-include-cleaner)
#include <iterator>    // IWYU pragma: keep
#include <limits>      // for std::numeric_limits  // IWYU pragma: keep
#include <ranges>      // IWYU pragma: keep
#include <tuple>       // for std::get  // NOLINT(misc-include-cleaner)

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

namespace lazy {

template<class A>
auto operator*(typename A::element scalar, A const& a) {
	return [scalar, &a](auto... is) { return scalar * a[is...]; } ^ a.extents();
}

namespace elementwise {

template<class A, class B>
auto operator*(A const& a, B const& b) requires(A::dimensionality == B::dimensionality) {
	return [&a, &b](auto... is) { return a[is...] * b[is...]; } ^ a.extents();
}

template<class A, class B>
auto operator+(A const& a, B const& b) requires(A::dimensionality == B::dimensionality) {
	return [&a, &b](auto... is) { return a[is...] + b[is...]; } ^ a.extents();
}

}  // namespace elementwise
}  // namespace lazy

struct re {
	void restrict() const {}
	// friend void restrict(re const&) {}
};

auto main() -> int {
	re R;
	R.restrict();
	// restrict(re);

	// test repeat
	{
		auto iota = [](multi::index i) { return i; } ^ multi::extents_t<1>(5);

		BOOST_TEST( iota.size() == 5 );
		BOOST_TEST( iota[0] == 0 );
		BOOST_TEST( iota[4] == 4 );

		auto iotax4 = iota.repeated(4);

		BOOST_TEST( iotax4.size() == 4 );

		BOOST_TEST( get<0>(iotax4.sizes()) == 4 );
		BOOST_TEST( get<1>(iotax4.sizes()) == 5 );

		BOOST_TEST( iotax4[0][0] == 0 );
		BOOST_TEST( iotax4[1][0] == 0 );
		BOOST_TEST( iotax4[2][0] == 0 );
		BOOST_TEST( iotax4[3][0] == 0 );

		BOOST_TEST( iotax4[0][1] == 1 );
		BOOST_TEST( iotax4[1][1] == 1 );
		BOOST_TEST( iotax4[2][1] == 1 );
		BOOST_TEST( iotax4[3][1] == 1 );
	}

	// subtract max
	{
		auto iota = [](multi::index i) { return i; } ^ multi::extents_t<1>(6);
		auto mat  = iota.partitioned(2);

		BOOST_TEST( mat.size() == 2 );
		BOOST_TEST( mat[0].size() == 3 );

		BOOST_TEST( mat[0][0] == 0 );
		BOOST_TEST( mat[0][1] == 1 );
		BOOST_TEST( mat[0][2] == 2 );

		BOOST_TEST( mat[1][0] == 3 );
		BOOST_TEST( mat[1][1] == 4 );
		BOOST_TEST( mat[1][2] == 5 );

		auto max_per_row = mat.transformed(maxR1);

		BOOST_TEST( max_per_row.num_elements() == 2 );
		BOOST_TEST( max_per_row.size() == 2 );

		auto max_per_row_repeat   = max_per_row.repeated(3);
		auto max_per_row_repeat_T = max_per_row_repeat.transposed();

		BOOST_TEST( max_per_row_repeat_T.size() == 2 );

		using multi::elementwise::operator-;

		auto subtract = mat - max_per_row_repeat_T;

		BOOST_TEST( subtract[0][0] == -2 );
		BOOST_TEST( subtract[0][1] == -1 );
		BOOST_TEST( subtract[0][2] ==  0 );

		BOOST_TEST( subtract[1][0] == -2 );
		BOOST_TEST( subtract[1][1] == -1 );
		BOOST_TEST( subtract[1][2] ==  0 );

		using multi::elementwise::exp;
		auto subtract_exp = exp(subtract);

		BOOST_TEST( subtract_exp.extents() == subtract.extents() );

		printR2("partial", exp(~(~mat - mat.transformed(maxR1))));

		auto&& exp_m_max = exp(mat - ~mat.transformed(maxR1).repeated(3));
		printR2("partial", exp_m_max);

		// auto rep3  = exp_m_max.transformed(sumR1).repeated(3);
		// auto final = exp_m_max / exp_m_max.transformed(sumR1).repeated(3);
		using multi::elementwise::operator/;
		using multi::elementwise::operator|;

		auto x = exp(~mat - (mat | maxR1));

		printR2("final", ~(x / ((~x) | sumR1)));

		BOOST_TEST( std::abs((~(x / ((~x) | sumR1)))[1][1] - 0.244728) < 1e-3 );
	}
	{
		auto lazy_one = [](multi::index) { return 1; } ^ multi::extents_t<1>{10};

		BOOST_TEST( lazy_one[3] == 1 );
		BOOST_TEST( decltype(lazy_one)::dimensionality == 1 );

		auto lazy_0D = []() { return 1; } ^ multi::extents_t<0>{};

		BOOST_TEST( lazy_0D[] == 1 );
		BOOST_TEST( decltype(lazy_0D)::dimensionality == 0 );

		auto rep = lazy_0D.repeated(5);
		BOOST_TEST( decltype(rep)::dimensionality == 1 );
		BOOST_TEST( rep.size() == 5 );
	}
	{
		multi::extents_t<1> xs(5);
		BOOST_TEST( xs.size() == 5 );

		auto curxs = xs.home();
		BOOST_TEST( curxs[2] == xs[2] );
	}
	{
		multi::extents_t<2> xs(3, 5);
		BOOST_TEST( xs.size() == 3 );

		auto curxs = xs.home();
		BOOST_TEST( curxs[0][0] == xs[0][0] );
		BOOST_TEST( curxs[1][1] == xs[1][1] );
	}
	{
		multi::extents_t<3> xs(3, 5, 7);
		BOOST_TEST( xs.size() == 3 );

		auto curxs = xs.home();
		BOOST_TEST( curxs[0][0][0] == xs[0][0][0] );
		BOOST_TEST( curxs[1][1][1] == xs[1][1][1] );
	}
	{
		auto v1D = [](auto ii) { return (ii * ii); } ^ multi::extents_t<1>(3);
		BOOST_TEST( v1D[2] == (2*2) );

		auto cur = v1D.home();
		BOOST_TEST( cur[2] == v1D[2] );
	}
	{
		auto v2D = [](auto ii, auto jj) { return (ii * ii) + (jj * jj); } ^ multi::extents_t<2>(3, 5);
		BOOST_TEST( v2D[2][3] == (2*2) + (3*3) );

		auto cur = v2D.home();
		BOOST_TEST( cur[2][2] == v2D[2][2] );

		// auto front = *v2D.begin();
		auto const v2D_diag = v2D.diagonal();
		BOOST_TEST( v2D_diag[2] == v2D[2][2] );
	}
	{
		std::initializer_list<int> il = {1, 2, 3};

		auto il_res = [il](auto ii) { return il.begin()[ii]; } ^ multi::extents_t<1>(static_cast<multi::ssize_t>(il.size()));
		BOOST_TEST( il_res[1] == 2 );
	}
	{
		std::initializer_list<std::initializer_list<int>> il = {
			{1, 2, 3},
			{4, 5, 6}
		};
		auto il_res = [il](auto ii, auto jj) { return il.begin()[ii].begin()[jj]; } ^ multi::extents(il);
		// multi::extents_t<2>(static_cast<multi::size_t>(il.size()), static_cast<multi::size_t>(il.begin()->size()));

		BOOST_TEST( il_res[1][1] == 5 );
	}
	{
		auto v2D = [](auto ii, auto jj) { return (ii * ii) + (jj * jj); } ^ multi::extents_t<2>(3, 5);

		auto it = v2D.begin();
		it++;
		BOOST_TEST( it == v2D.begin() + 1 );

		multi::array<multi::index, 2> v2D_copy = v2D;

		multi::array<multi::index, 2> v2D_copy2(v2D.extents());

		v2D_copy2() = v2D;
		BOOST_TEST( v2D_copy2 == v2D_copy );

		multi::array<multi::index, 2> v2D_copy3;
		v2D_copy3 = v2D;

		BOOST_TEST( v2D_copy3 == v2D_copy );

		// multi::array<multi::index, 2> v2D_copy4;

		// v2D_copy4() = v2D;  // this fails with an assert because sizes do not match
		// BOOST_TEST( v2D_copy4 == v2D_copy );

		static_assert(std::is_default_constructible_v<decltype(v2D)::iterator>);
		static_assert(std::is_default_constructible_v<decltype(v2D.elements())::iterator>);
	}
	{
		auto v2D = [](auto ii, auto jj) { return (ii * ii) + (jj * jj); } ^ multi::extents_t<2>(5, 7);
		BOOST_TEST( (v2D.begin() + 3) - v2D.begin() == 3 );
		BOOST_TEST( ((v2D.begin() + 2) + 1) - v2D.begin() == 3 );
	}
	{
		multi::iextension m(96);
		multi::iextension h(64);
		multi::iextension k(64);
		multi::iextension n(96);

		multi::array<float, 4> A = +(  // NOLINTNEXTLINE(runtime/threadsafe_fn)
			[](auto...) { return (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5f) * 100.0f; } ^ multi::extents_t<4>{m, h, k, n}
		);
	}
	{
		multi::array<double, 3> arr;
		arr = [](auto i, auto j, auto k) { return static_cast<double>(i + j + k); } ^ multi::extents_t<3>(2, 3, 4);
	}
	{
		std::vector<int> vec = {1, 2, 3};

		auto&& arr = [&vec](auto i) -> int& { return vec[static_cast<std::size_t>(i)]; } ^ multi::extents_t<1>(static_cast<multi::extents_t<1>::size_type>(vec.size()));

		arr[1] = 99;

		BOOST_TEST( vec[1] == 99 );
	}
	{
		std::vector<int> vec = {1, 2, 3};

		auto&& arr = [&vec](auto i, auto /*j*/) -> int& { return vec[static_cast<std::size_t>(i)]; } ^ multi::extents_t<2>(static_cast<multi::extents_t<1>::size_type>(vec.size()), static_cast<multi::extents_t<1>::size_type>(vec.size()));

		arr[1][1] = 99;

		BOOST_TEST( vec[1] == 99 );
	}

	{
		std::vector<std::vector<int>> vvec = {
			{1, 2, 3},
			{4, 5, 6}
		};

		auto&& arr = [&vvec](auto i, auto j) -> int& { return vvec[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]; } ^ multi::extents_t<2>(static_cast<multi::ssize_t>(vvec.size()), static_cast<multi::ssize_t>(vvec.front().size()));

		arr[1][1] = 99;

		BOOST_TEST( vvec[1][1] == 99 );
	}

	return boost::report_errors();
}
#else
auto main() -> int {
	return boost::report_errors();
}
#endif
