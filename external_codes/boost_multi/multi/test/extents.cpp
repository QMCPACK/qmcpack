// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/restriction.hpp>

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <algorithm>  // IWYU pragma: keep  // for std::equal
#include <iterator>   // IWYU pragma: keep

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#include <concepts>  // for totally_ordered
#include <ranges>    // IWYU pragma: keep
#endif

#include <tuple>        // IWYU pragma: keep
#include <type_traits>  // for std::is_same_v
// IWYU pragma: no_include <utility>  // for declval, forward, move
// IWYU pragma: no_include <variant>  // for get, iwyu bug
#include <vector>

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep  // for _GLIBCXX_RELEASE, __GLIBC...
#endif

namespace multi = boost::multi;

// Observed failure: clang-15 (as its own toolchain, not later clang versions) paired with
// GCC-11's libstdc++ (_GLIBCXX_RELEASE 11) fails to see boost::multi's iterators as satisfying
// `range`/`input_or_output_iterator` when instantiated through std::ranges::ref_view (as
// happens with std::views::reverse and other pipe adaptors), even though the iterator is
// well-formed and works with std::ranges::begin/end directly. See test/broadcast_softmax.cpp
// for the same guard.
#if defined(__clang__) && (__clang_major__ == 15) && defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE < 12)
#define BOOST_MULTI_STDRANGES_PIPE_BROKEN 1
#endif

auto main() -> int {  // NOLINT(bugprone-exception-escape,readability-function-cognitive-complexity)
	auto const A2D = multi::array<int, 2>({5, 7}, 1);

	// auto const A2Dx = A2D.extension();
	auto const A2Dx = A2D.extent();

	BOOST_TEST( &A2D() == &A2D(A2Dx) );

	auto const A2Dxs = A2D.extents();

	using std::get;

	BOOST_TEST( get<0>(A2Dxs[1][2]) == 1 );
	BOOST_TEST( get<1>(A2Dxs[1][2]) == 2 );

	BOOST_TEST( get<0>(A2Dxs) == A2Dx );
	BOOST_TEST( get<1>(A2Dxs) == A2D[0].extent() );

	BOOST_TEST( &A2D() == &A2D(get<0>(A2D.extents()), get<1>(A2D.extents())) );
	BOOST_TEST( &A2D() == &std::apply(A2D, A2Dxs) );

	BOOST_TEST( A2Dxs.size() == A2D.size() );
	BOOST_TEST( A2Dxs.sizes() == A2D.sizes() );

	auto const [ni, nj] = A2Dxs.sizes();
	for(int i = 0; i != ni; ++i) {      // NOLINT(altera-unroll-loops)
		for(int j = 0; j != nj; ++j) {  // NOLINT(altera-unroll-loops)
			auto const [first, second] = A2Dxs[i][j];
			BOOST_TEST( first == i );
			BOOST_TEST( second == j );
		}
	}

	auto const [is, js] = A2Dxs;
	for(auto const i : is) {      // NOLINT(altera-unroll-loops)
		for(auto const j : js) {  // NOLINT(altera-unroll-loops)
			auto const [first, second] = A2Dxs[i][j];
			BOOST_TEST( first == i );
			BOOST_TEST( second == j );
		}
	}

	BOOST_TEST( get<0>(A2Dxs).size() == 5 );

	// auto it2d = A2Dxs.elements().begin(); (void)it2d;

	multi::array<int, 1> const A1D({37}, 1);
	BOOST_TEST( A1D.size() == 37 );
	BOOST_TEST( A1D.num_elements() == 37 );
	BOOST_TEST( A1D.extents().num_elements() == 37 );

	BOOST_TEST( A1D.extents().elements().size() == A1D.extents().num_elements() );
	{
		auto it = A1D.extents().elements().begin();
		BOOST_TEST( get<0>(*it) == 0 );
		++it;
		BOOST_TEST( get<0>(*it) == 1 );

		it = A1D.extents().elements().end();
		--it;
		BOOST_TEST( get<0>(*it) == 36 );
	}
	{
		auto x1d = multi::extents_t<1>(3);

		BOOST_TEST( multi::extents_t<1>(3) == multi::extents_t(3) );

		auto it = x1d.elements().begin();
		BOOST_TEST( get<0>(*it) == 0 );

		++it;
		BOOST_TEST( get<0>(*it) == 1 );

		++it;
		BOOST_TEST( get<0>(*it) == 2 );

		++it;
		BOOST_TEST( it == x1d.elements().end() );

		--it;
		BOOST_TEST( get<0>(*it) == 2 );

		--it;
		BOOST_TEST( get<0>(*it) == 1 );

		--it;
		BOOST_TEST( get<0>(*it) == 0 );
		BOOST_TEST( it == x1d.elements().begin() );

		BOOST_TEST( x1d.elements().begin() != x1d.elements().end() );
		BOOST_TEST( !(x1d.elements().begin() == x1d.elements().end()) );

		BOOST_TEST( x1d.elements().begin() < x1d.elements().end() );
		BOOST_TEST( x1d.elements().begin() <= x1d.elements().end() );
		BOOST_TEST( x1d.elements().begin() <= x1d.elements().begin() );  // cppcheck-suppress [duplicateExpression];
	}
	{
		auto x1d = multi::extents_t<1>(3);

		auto it = x1d.elements().begin();
		BOOST_TEST( get<0>(*it) == 0 );

		++it;
		BOOST_TEST( get<0>(*it) == 1 );

		it += 2;

		BOOST_TEST( it == x1d.elements().end() );

		it -= 3;
		BOOST_TEST( get<0>(*it) == 0 );
		BOOST_TEST( it == x1d.elements().begin() );

		BOOST_TEST( x1d.elements().end() - x1d.elements().begin() == 3 );
	}
	{
		multi::extents_t<2> const x2d({4, 3});

		BOOST_TEST( multi::extents_t<2>(4, 3) == multi::extents_t(4, 3) );

		auto ll = [](auto xx, auto yy) {
			return xx + yy;
		};
		multi::restriction<2, decltype(ll)> const x2df({4, 2}, ll);
		(void)x2df;
		auto val = x2df[3][1];
		BOOST_TEST(val == 4);

		auto elems = x2df.elements();
		BOOST_TEST( elems[7] == 4 );
		BOOST_TEST( *(x2df.elements().begin() + 1) == 1 + 0 );

		// BOOST_TEST( *(*(x2df.begin()).begin()) == 0 )

		// multi::detail::what(x2df[1]);
		// std::cout << x2df[1][2] << std::endl;

		// auto x2d_trd = x2d.element_transformed([](auto is) { using std::get; return get<0>(is) + get<1>(is); });

		BOOST_TEST( x2d.elements().end() - x2d.elements().begin() == 12 );

		auto it = x2d.elements().begin();

		BOOST_TEST( it == x2d.elements().begin() );

		using std::get;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		++it;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		++it;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		BOOST_TEST( it - x2d.elements().begin() == 2 );

		++it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		++it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		++it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		BOOST_TEST( it - x2d.elements().begin() ==  5 );
		BOOST_TEST( x2d.elements().begin() - it == -5 );
		BOOST_TEST( x2d.elements().end() - it == 7 );
		BOOST_TEST( x2d.elements().end() - x2d.elements().begin() == 12 );

		++it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		++it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		++it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		++it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		++it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		++it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		++it;
		BOOST_TEST( it ==  x2d.elements().end() );

		--it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		--it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		--it;
		BOOST_TEST( 3 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		--it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		--it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		--it;
		BOOST_TEST( 2 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		--it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		--it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		--it;
		BOOST_TEST( 1 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		--it;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 2 == get<1>(*it) );

		--it;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 1 == get<1>(*it) );

		--it;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		BOOST_TEST( it ==  x2d.elements().begin() );
	}
	{
		multi::extents_t<2> const x2d({4, 3});

		// auto it2d = x2d.begin();

		auto it = x2d.elements().begin();

		BOOST_TEST( it == x2d.elements().begin() );

		using std::get;
		BOOST_TEST( 0 == get<0>(*it) );
		BOOST_TEST( 0 == get<1>(*it) );

		BOOST_TEST( 0 == get<0>(*(it + 2)) );
		BOOST_TEST( 2 == get<1>(*(it + 2)) );

		BOOST_TEST( 1 == get<0>(*(it + 5)) );
		BOOST_TEST( 2 == get<1>(*(it + 5)) );

		auto const it2  = it + 5;
		auto const it22 = it - (-5);
		BOOST_TEST( it2 == it22 );

		BOOST_TEST( 1 == get<0>(*(it2)) );
		BOOST_TEST( 2 == get<1>(*(it2)) );

		BOOST_TEST( 1 == get<0>(*(it2-1)) );
		BOOST_TEST( 1 == get<1>(*(it2-1)) );

		BOOST_TEST( 1 == get<0>(*(it2-2)) );
		BOOST_TEST( 0 == get<1>(*(it2-2)) );

		auto const it3  = it2 - 5;
		auto const it33 = it2 + (-5);
		BOOST_TEST( it3 == it33 );

		BOOST_TEST( it3 == it );
	}
	// {
	// 	multi::array<int, 1> const arr(10);

	// 	auto xn = decltype(arr.extent())(10);
	// 	BOOST_TEST( xn. size() == 10 );

	// 	multi::extent_t const xn2(10);
	// 	BOOST_TEST( xn2.size() == 10 );

	// 	xn = xn2;

	// 	multi::detail::extensions const xns2{xn2};
	// 	using std::get;
	// 	BOOST_TEST( get<0>(xns2) == xn2 );

	// 	multi::detail::extensions const xns2d{xn2, xn2};
	// 	auto [xns2d_a, xns2d_b] = xns2d;

	// 	BOOST_TEST( xns2d_a == xn2 );
	// 	BOOST_TEST( xns2d_b == xn2 );

	// 	multi::extents_t<2> const met2{xns2d};

	// 	multi::layout_t<2> const lyt(met2);
	// 	multi::layout_t<2> const lyt_2(xns2d);

	// 	BOOST_TEST( lyt == lyt_2 );

	// 	// multi::array<int, 1> const arr2({xn2});
	// }
	{
		auto const x2df = [](auto x, auto y) { return x + y; } ^ multi::extents_t<2>(3, 4);

		// boost::multi::restriction<2, decltype(ll)> x2df(multi::extents_t<2>(3, 4), ll);
		BOOST_TEST( x2df.elements()[0] == 0 );
		BOOST_TEST( x2df.elements()[1] == 1 );
		BOOST_TEST( x2df.elements()[2] == 2 );
		BOOST_TEST( x2df.elements()[3] == 3 );
		BOOST_TEST( x2df.elements()[4] == 1 );
		BOOST_TEST( x2df.elements()[5] == 2 );

		BOOST_TEST( x2df[2][1] == 2 + 1 );

		multi::array<multi::index, 2> const arr2df = [](auto x, auto y) { return x + y; } ^ multi::extents_t<2>(3, 4);

		BOOST_TEST( arr2df(2, 1) == 2 + 1 );
		BOOST_TEST( arr2df[2][1] == 2 + 1 );

		BOOST_TEST(std::equal(
			arr2df.elements().begin(), arr2df.elements().end(),
			([](auto x, auto y) { return x + y; } ^ multi::extents_t<2>(3, 4)).elements().begin()
		));

		using multi::experimental::operator->*;
		BOOST_TEST(std::equal(
			arr2df.elements().begin(), arr2df.elements().end(),
			(multi::extents_t<2>(3, 4)->*[](auto x, auto y) { return x + y; }).elements().begin()
		));

		BOOST_TEST(   arr2df.elements().begin() != arr2df.elements().end()  );
		BOOST_TEST( !(arr2df.elements().begin() == arr2df.elements().end()) );

		BOOST_TEST( arr2df[2][1] == ([](auto x, auto y) { return x + y; } ^ multi::extents_t<2>(3, 4))[2][1] );

		BOOST_TEST( arr2df[2][1] == ([](auto x, auto y) { return x + y; } ^ multi::extents_t(3, 4))[2][1] );
		BOOST_TEST(
			arr2df[2][1]
			== multi::extents_t<2>(3, 4).element_transformed( [](auto const& idxs) { using std::get; return get<0>(idxs) + get<1>(idxs); })[2][1]
		);
		BOOST_TEST(
			arr2df[2][1]
			== multi::extents_t<2>(3, 4).element_transformed( [](auto idxs) {auto [xx, yy] = idxs; return xx + yy; })[2][1]
		);
	}
	{
		multi::extents_t<3> const xs{3, 4, 5};

		BOOST_TEST(( multi::extents_t<3>{3, 4, 5} == multi::extents_t(3, 4, 5) ));

		BOOST_TEST( xs.sub() == multi::extents_t<2>(4, 5) );
		static_assert(std::is_same_v<decltype(xs[1][1][1]), multi::extents_t<3>::element>);
	}
	{
		multi::array<int, 2> const arr({3, 4});

		auto const& xs = arr.extents();

		using std::get;
		BOOST_TEST( get<0>(xs[0][0]) == 0 );
		BOOST_TEST( get<1>(xs[0][0]) == 0 );

		BOOST_TEST(   xs.begin() != xs.end()  );
		BOOST_TEST( !(xs.begin() == xs.end()) );

		BOOST_TEST( xs[0] == xs[0] );
		BOOST_TEST( xs[0] != xs[1] );

		BOOST_TEST( xs[0] == *xs.begin() );
		BOOST_TEST( xs[1] == *(xs.begin() + 1) );

		auto it = xs.begin();
		++it;
		BOOST_TEST( *it == xs[1] );

		auto const& values = [](auto ii, auto jj) { return ii + jj; } ^ arr.extents();

		BOOST_TEST( values.dimensionality == 2 );
		BOOST_TEST( values.extents() == arr.extents() );
		BOOST_TEST( *values.elements().begin() == 0 );
		BOOST_TEST( values.elements().begin() < values.elements().end() );
		BOOST_TEST( values.elements().begin() != values.elements().end() );
		BOOST_TEST( values[0][0] == 0 );
		BOOST_TEST( values.begin() != values.end() );

		{
			auto arr2 = multi::array<boost::multi::index, 2>(arr.extents());

			arr2.elements() = values.elements();
			BOOST_TEST( std::equal(arr2.elements().begin(), arr2.elements().end(), values.elements().begin(), values.elements().end()) );
		}
		{
			auto arr2 = multi::array<boost::multi::index, 2>(arr.extents());

			arr2() = values;
			BOOST_TEST( std::equal(arr2.elements().begin(), arr2.elements().end(), values.elements().begin(), values.elements().end()) );
		}
		{
			auto arr2 = multi::array<boost::multi::index, 2>(arr.extents());

			arr2 = values;
			BOOST_TEST( std::equal(arr2.elements().begin(), arr2.elements().end(), values.elements().begin(), values.elements().end()) );
		}

#ifdef __cpp_deduction_guides
		{
			multi::array<multi::index, 2> const arr_gold = values;
			multi::array const                  arr2     = values;
			BOOST_TEST( arr_gold == arr2 );
		}
#endif
	}
	{
		auto xs1D = multi::extents_t(10);
		BOOST_TEST( xs1D.size() == 10 );
		using std::get;
		BOOST_TEST( get<0>(xs1D[3]) == 3 );

		BOOST_TEST( xs1D.begin() != xs1D.end() );
		BOOST_TEST( !(xs1D.begin() == xs1D.end()) );
		BOOST_TEST( xs1D.begin() + 10 == xs1D.end() );
		BOOST_TEST( xs1D.begin() == xs1D.end() - 10 );

		BOOST_TEST( *(xs1D.begin() + 3) == xs1D[3] );

#ifdef __NVCC__  // nvcc gets confused with inline lambdas
		auto fun = [](auto ii) noexcept { return ii * ii; };
		auto v1D = fun ^ multi::extents_t(10);
#else
		auto v1D = [](auto ii) noexcept { return ii * ii; } ^ multi::extents_t(10);
#endif

		BOOST_TEST( v1D.size() == 10 );
		BOOST_TEST( v1D[4] == 16 );

		BOOST_TEST( v1D.elements().size() == 10 );
		BOOST_TEST( v1D.elements()[4] == v1D[4] );

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<1>::iterator>);
		static_assert(std::random_access_iterator<multi::extents_t<1>::iterator>);
		static_assert(std::ranges::random_access_range<decltype(xs1D)>);

		BOOST_TEST( xs1D.begin() == std::ranges::begin(xs1D) );
		BOOST_TEST( xs1D.end()   == std::ranges::end(xs1D)   );

#ifndef BOOST_MULTI_STDRANGES_PIPE_BROKEN
		auto xs1Dr = xs1D | std::ranges::views::reverse;

		BOOST_TEST( *xs1Dr.begin() == 9 );
		BOOST_TEST( *(xs1Dr.end() - 1) == 0 );

		BOOST_TEST( xs1Dr[9] == xs1D[0]	);
		BOOST_TEST( xs1Dr[0] == xs1D[9]	);
#endif

		// auto xs1D_elements = xs1D.elements();
		BOOST_TEST( xs1D.elements().begin() == std::ranges::begin(xs1D.elements()) );

		static_assert(std::input_or_output_iterator<decltype(v1D)::iterator>);

		BOOST_TEST( std::ranges::begin(v1D) == v1D.begin() );
		BOOST_TEST( std::ranges::end(v1D) == v1D.end() );

		static_assert(std::totally_ordered<decltype(v1D)::iterator>);
		static_assert(std::random_access_iterator<decltype(v1D)::iterator>);

#ifndef BOOST_MULTI_STDRANGES_PIPE_BROKEN
		auto v1Dr = v1D | std::views::reverse;
		BOOST_TEST( v1Dr[0] == v1D[9] );
		BOOST_TEST( v1Dr[9] == v1D[0] );
#endif
#endif
	}
	{
		auto xs2D = multi::extents_t<2>(5, 7);
		BOOST_TEST( xs2D.size() == 5 );

		using std::get;
		BOOST_TEST( get<0>(xs2D[3][2]) == 3 );
		BOOST_TEST( get<1>(xs2D[3][2]) == 2 );

		BOOST_TEST( xs2D.begin() != xs2D.end() );
		BOOST_TEST( !(xs2D.begin() == xs2D.end()) );
		BOOST_TEST( xs2D.begin() + xs2D.size() == xs2D.end() );
		BOOST_TEST( xs2D.begin() == xs2D.end() - xs2D.size() );

		BOOST_TEST( *(xs2D.begin() + 3) == xs2D[3] );

		// auto it = xs2D.begin();
		// multi::detail::what(*it);

		auto xs3D = multi::extents_t<3>(5, 7, 21);
		BOOST_TEST( xs3D.size() == 5 );
		// multi::detail::what(*xs3D.begin());

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
		using xs2D_iterator = multi::extents_t<2>::iterator;

		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<0>>);
		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<2>>);

		static_assert(std::is_trivially_default_constructible_v<multi::range<multi::index, multi::index>>);
		static_assert(std::is_trivially_default_constructible_v<multi::extent_t<multi::index, multi::index>>);

		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<1>::base_>);
		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<1>>);
		static_assert(std::is_trivially_default_constructible_v<multi::extents_t<2>::iterator>);

		static_assert(std::is_constructible_v<xs2D_iterator>);
		static_assert(std::constructible_from<xs2D_iterator, xs2D_iterator>);
		static_assert(std::default_initializable<multi::extents_t<2>::iterator>);
		static_assert(std::semiregular<multi::extents_t<2>::iterator>);
		static_assert(std::regular<multi::extents_t<2>::iterator>);
		static_assert(std::incrementable<multi::extents_t<2>::iterator>);

		static_assert(std::weakly_incrementable<multi::extents_t<2>::iterator>);
		static_assert(std::input_iterator<multi::extents_t<2>::iterator>);
		static_assert(std::forward_iterator<multi::extents_t<2>::iterator>);
		static_assert(std::bidirectional_iterator<multi::extents_t<2>::iterator>);
		static_assert(std::random_access_iterator<multi::extents_t<2>::iterator>);
		static_assert(std::ranges::random_access_range<multi::extents_t<2>>);

		BOOST_TEST( xs2D.begin() == std::ranges::begin(xs2D) );
		BOOST_TEST( xs2D.end()   == std::ranges::end(xs2D)   );

#ifndef BOOST_MULTI_STDRANGES_PIPE_BROKEN
		auto xs2Dr = xs2D | std::ranges::views::reverse;

		BOOST_TEST( *xs2Dr.begin() == *(xs2D.end() - 1) );
		BOOST_TEST( *(xs2Dr.end() - 1) == *(xs2D.begin()) );

		BOOST_TEST( xs2Dr[xs2D.size() - 1] == xs2D[0] );
		BOOST_TEST( xs2Dr[0] == xs2D[xs2D.size() - 1] );
#endif
#endif
	}
	{
		auto v2D = [](auto ii, auto jj) { return (ii * ii) + (jj * jj); } ^ multi::extents_t<2>(3, 5);
		BOOST_TEST( v2D[2][3] == (2*2) + (3*3) );
		// auto front = *v2D.begin();

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
		BOOST_TEST( v2D.begin() == std::ranges::begin(v2D) );
		BOOST_TEST( v2D.end()   == std::ranges::end(v2D)   );

#ifndef BOOST_MULTI_STDRANGES_PIPE_BROKEN
		auto v2Dr = v2D | std::ranges::views::reverse;

		BOOST_TEST( (*v2Dr.begin())[4] == (*(v2D.end() - 1))[4] );
		BOOST_TEST( (*(v2Dr.end() - 1))[4] == (*(v2D.begin()))[4] );

		BOOST_TEST( v2Dr[v2D.size() - 1][5] == v2D[0][5] );
		BOOST_TEST( v2Dr[0][5] == v2D[v2D.size() - 1][5] );
#endif

		// auto const v2DT = v2D.transposed() | std::views::reverse;  // TODO(correaa)
		// BOOST_TEST( v2DT[1][5] == v2D[2][1] );
		{
			auto matrix =
				([](auto ii) noexcept { return static_cast<float>(ii); } ^
				 multi::extents_t(6))
					.partitioned(2);

			auto [matrix_is, matrix_js] = matrix.extents();
			BOOST_TEST( matrix_is.size() == 2 );
			BOOST_TEST( matrix_js.size() == 3 );

#ifndef __NVCC__
			static_assert(std::movable<std::decay_t<decltype(matrix[0])>::iterator>);
#endif
			BOOST_TEST( std::ranges::begin(matrix[0]) == matrix[0].begin() );
			BOOST_TEST( std::ranges::end(matrix[0]) == matrix[0].end() );
		}
#endif
	}
	{
		multi::extents_t<2> const x2D(6, 5);
		multi::extents_t<3> const p3D = multi::layout_t<2>(x2D).partition(2).extents();

		using std::get;
		BOOST_TEST( get<0>(p3D).size() == 2 );
		BOOST_TEST( get<1>(p3D).size() == 3 );
		BOOST_TEST( get<2>(p3D).size() == 5 );
	}
	{
		auto exts = multi::extents_t<2>(3, 4);

		// auto something = exts[-1];
		// (void)something;

		BOOST_TEST( exts.extents() == exts );

		static_assert(std::is_default_constructible_v<decltype(exts.elements())::iterator>);

		auto it = exts.elements().begin() + 8;
		BOOST_TEST( it + 0 == it );  // cppcheck-suppress knownConditionTrueFalse

		it += -6;
		BOOST_TEST( it == exts.elements().begin() + 2 );
	}
	{
		auto exts = multi::extents_t<2>(3, 4);

		{
			auto [ext1, ext2] = exts;

			BOOST_TEST( ext1.size() == 3 );
			BOOST_TEST( ext2.size() == 4 );
		}
		{
			auto [ext1, ext2] = exts.transpose();

			BOOST_TEST( ext1.size() == 4 );
			BOOST_TEST( ext2.size() == 3 );
		}

		BOOST_TEST( exts == exts.transpose().transpose() );
	}
	{
		auto exts = multi::extents_t<1>(10);
		BOOST_TEST( exts.size() == 10 );
		BOOST_TEST( (exts.end() - 1) - (exts.begin() + 1) == exts.size() - 2 );

		auto it = exts.begin();
		it += 2;
		BOOST_TEST( it == exts.begin() + 2 );

		BOOST_TEST( exts.elements().size() == 10 );
		BOOST_TEST( (exts.elements().end() - 1) - (exts.elements().begin() + 1) == exts.elements().size() - 2 );
		BOOST_TEST( (exts.elements().end() - 1) - (exts.elements().begin() + 2) == exts.elements().size() - 3 );

		auto it2 = exts.begin();
		auto it3 = it2 + 3;
		auto it4 = it3 + 2;

		BOOST_TEST( it4 == exts.begin() + 5 );
		BOOST_TEST( *it4 == 5 );

		auto it5 = it4 - 3;
		BOOST_TEST( it5 == exts.begin() + 2 );
		BOOST_TEST( *it5 == 2 );
	}
	{
		auto exts = multi::extents_t<1>(10);

		// BOOST_TEST( exts[-1] != decltype(exts[-1]){} );  gives an out-of-bounds assert

		BOOST_TEST( exts.size() == 10 );
		BOOST_TEST( (exts.end() - 1) - (exts.begin() + 1) == exts.size() - 2 );
	}
	{
		std::vector<int> const vec(20);

		multi::array<int, 2> const arr(multi::extents_t<2>(vec.size(), vec.size()), 99);
	}
	{
		std::vector<int> const vec(20);

		multi::array<int, 2> const arr(multi::extents_t(vec.size(), vec.size()), 99);
	}
	{
		std::vector<int> const vec(20);

		multi::array<int, 2> const arr({static_cast<multi::ssize_t>(vec.size()), static_cast<multi::ssize_t>(vec.size())}, 99);  // needs to be explicit, good
	}
	// {
	// 	std::vector<int> v(20);
	// 	multi::array<int, 2> A({v.size(), v.size()}, 99);  // error, needs to be explicit, good
	// }

	return boost::report_errors();
}
