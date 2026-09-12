// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/elementwise.hpp>
#include <boost/multi/restriction.hpp>  // for restriction, operator!=

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <algorithm>   // IWYU pragma: keep  // for std::equal
#include <cmath>       // for std::abs
#include <functional>  // for std::plus  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep
#include <iostream>
#include <iterator>  // IWYU pragma: keep
#include <limits>    // for std::numeric_limits  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep
#include <numeric>
// IWYU pragma: no_include <tuple>    // for apply
#include <utility>  // for forward  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep  // for _GLIBCXX_RELEASE, __GLIBC...
#endif

namespace multi = boost::multi;

// NOLINTBEGIN(readability-identifier-length)
auto main() -> int {  // NOLINT(readability-function-cognitive-complexit,bugprone-exception-escape,readability-function-cognitive-complexity)
	{
		// multi::array const a = {1.0, 2.0, 3.0};
		multi::array<double, 1> const a = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		auto c = exp(a);

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		multi::array<double, 1> const a = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		auto c = exp(a);

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		multi::array<double, 1> const a = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		auto c = exp(a);

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		multi::array<double, 1> const a = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		auto c = exp(a);

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		multi::array<double, 1> a = {1.0, 2.0, 3.0};

		using multi::elementwise::exp;
		auto c = exp(a());

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		auto r = [](auto i) constexpr { return static_cast<double>(i + 1); } ^ multi::extents_t<1>{3};
		using multi::elementwise::exp;
		auto c = exp(r);

		BOOST_TEST( c.extents() == r.extents() );

		BOOST_TEST( std::abs(c[0] - std::exp(1.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[1] - std::exp(2.0)) < 1e-4 );
		BOOST_TEST( std::abs(c[2] - std::exp(3.0)) < 1e-4 );
	}
	{
		multi::array<int, 1> a = {1, -2, -3};

		using multi::elementwise::abs;
		auto const& c = abs(a);

		BOOST_TEST( c.extents() == a.extents() );

		BOOST_TEST( c[0] == std::abs(a[0]) );
		BOOST_TEST( c[1] == std::abs(a[1]) );
		BOOST_TEST( c[2] == std::abs(a[2]) );

		multi::array<int, 1> const c_copy1 = c;

		BOOST_TEST( c_copy1.extents() == c.extents() );

		multi::array const c_copy2 = c;
		BOOST_TEST( c_copy2 == c_copy1 );

		auto const c_copy3 = +c;
		BOOST_TEST( c_copy3 == c_copy1 );
		BOOST_TEST( c_copy3.base() != nullptr );
	}
	// {
	// 	using multi::elementwise::exp;
	// 	auto c = exp(
	// 		{{1.0, 2.0, 3.0},
	// 		{4.0, 5.0, 6.0}}
	// 	);

	// 	BOOST_TEST( std::abs(c[0][0] - std::exp(1.0)) < 1e-4 );
	// 	BOOST_TEST( std::abs(c[0][1] - std::exp(2.0)) < 1e-4 );
	// 	BOOST_TEST( std::abs(c[0][2] - std::exp(3.0)) < 1e-4 );

	// 	BOOST_TEST( std::abs(c[1][0] - std::exp(4.0)) < 1e-4 );
	// 	BOOST_TEST( std::abs(c[1][1] - std::exp(5.0)) < 1e-4 );
	// 	BOOST_TEST( std::abs(c[1][2] - std::exp(6.0)) < 1e-4 );
	// }
	{
		multi::array<int, 1> const a = {-1, -2, 3};

		using multi::elementwise::abs;
		auto c = abs(a);

		BOOST_TEST( c[0] == 1 );
		BOOST_TEST( c[1] == 2 );
		BOOST_TEST( c[2] == 3 );
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};
		multi::array<int, 1> const b = {4, 5, 6};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto c = a + b;

		BOOST_TEST(( c == multi::array<int, 1>{5, 7, 9} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};
		multi::array<int, 1> const b = {4, 5, 6};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto const& c = a + b;

		BOOST_TEST(( c == multi::array<int, 1>{5, 7, 9} ));
	}
	{
		auto const A = multi::array<int, 2>{
			{0, 1, 2},
			{3, 4, 5}
		};
		auto const B = multi::array<int, 2>{
			{ 0, 10, 20},
			{30, 40, 50}
		};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const C = A + B;

		BOOST_TEST( C[1][1] == A[1][1] + B[1][1] );
	}
	{
		auto const A = multi::array<int, 1>{0, 1, 2};
		auto const B = multi::array<int, 1>{0, 10, 20};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];
		using multi::elementwise::operator*;  // cppcheck-suppress [constStatement];

		multi::array<int, 1> const C = A + (2 * B);

		BOOST_TEST( C[1] == A[1] + (2*B[1]) );
	}
	{
		auto const A = multi::array<int, 2>{
			{0, 1},
			{2, 3}
		};
		auto const B = multi::array<int, 2>{
			{ 0, 10},
			{20, 30}
		};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];
		using multi::elementwise::operator*;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const C = A + (2 * B);

		BOOST_TEST( C[1][1] == A[1][1] + (2 * B[1][1]) );
	}
	{
		auto const A = multi::array<int, 2>{
			{0, 1, 2},
			{3, 4, 5}
		};
		auto const B = multi::array<int, 2>{
			{ 0, 10, 20},
			{30, 40, 55}
		};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];
		using multi::elementwise::operator*;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const C = A + (2 * B);

		BOOST_TEST( C[1][1] == A[1][1] + (2 * B[1][1]) );
	}
	{
		auto const A = multi::array<int, 2>{
			{0, 1, 2},
			{3, 4, 5}
		};
		auto const B = multi::array<int, 2>{
			{ 0, 10, 20},
			{30, 40, 50}
		};
		auto const C = multi::array<int, 2>{
			{  0, 100, 200},
			{300, 400, 500}
		};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const D = A + B + C;

		BOOST_TEST( D[1][1] == A[1][1] + B[1][1] + C[1][1] );
	}
	{
		auto const A = multi::array<int, 2>{
			{0, 1, 2},
			{3, 4, 5}
		};
		auto const B = multi::array<int, 2>{
			{ 0, 10, 20},
			{30, 40, 50}
		};

		using multi::elementwise::operator*;  // cppcheck-suppress [constStatement];
		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const C = A + (A * B);

		BOOST_TEST( C[1][1] == A[1][1] + (A[1][1] * B[1][1]) );
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		auto f1d = [](auto) { return 1; } ^ multi::extents_t<1>{3};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto const& c = a + f1d;

		BOOST_TEST(( multi::array<int, 1>{2, 3, 4} == c ));
		BOOST_TEST(( c == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		auto f = []() { return 1; } ^ multi::extents_t<0>{};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto const& c = a + f;

		BOOST_TEST(( multi::array<int, 1>{2, 3, 4} == c ));
		BOOST_TEST(( c == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto const& c = a + ([]() { return 1; } ^ multi::extents_t<0>{});

		BOOST_TEST(( multi::array<int, 1>{2, 3, 4} == c ));
		BOOST_TEST(( c == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		std::cout << (a + 1)[1] << '\n';
		BOOST_TEST(( multi::array<int, 1>{2, 3, 4} == a + 1 ));
		BOOST_TEST(( a + 1 == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

		auto c = a + 1;

		BOOST_TEST(( multi::array<int, 1>{2, 3, 4} == c ));
		BOOST_TEST(( c == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 2> const A = {
			{1, 2, 3},
			{4, 5, 6}
		};

		multi::array<int, 1> const b = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const C = A + b;

		BOOST_TEST((
			C ==
			multi::array<int, 2>{
				{2, 4, 6},
				{5, 7, 9}
			}
		));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];

		multi::array<int, 1> const b = a + 1;

		BOOST_TEST(( b == multi::array<int, 1>{2, 3, 4} ));
	}
	{
		multi::array<int, 1> const a = {1, 2, 3};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];
		BOOST_TEST(( a + 1 == multi::array<int, 1>{2, 3, 4} ));
	}
#if (!defined(__GNUC__) || (__GNUC__ > 8)) || defined(__clang__)
	{
		multi::array<int, 2> const A = {
			{0, 1, 2},
			{3, 4, 5}
		};
		multi::array<int, 2> const B = {
			{0, 1, 2},
			{3, 4, 5}
		};
		multi::array<int, 2> const C = {
			{0, 1, 2},
			{3, 4, 5}
		};

		using multi::elementwise::operator+;  // cppcheck-suppress [constStatement];
		using multi::elementwise::operator*;  // cppcheck-suppress [constStatement];

		multi::array<int, 2> const D = A + A * B + 2 * C;

		BOOST_TEST( D[1][1] == A[1][1] + (A[1][1] * B[1][1]) + (2 * C[1][1]) );

		auto const& r = (A + A * B + 2 * C).diagonal();

		auto trace_D = std::accumulate(r.begin(), r.end(), 0);  // NOLINT(misc-include-cleaner) std::reduce unavailable in libstdc++ < 9 (e.g. clang-8 CI)

		BOOST_TEST(trace_D == std::accumulate(D.diagonal().begin(), D.diagonal().end(), 0) );
	}
#endif
	{
		using multi::elementwise::eye;
		auto arr = +eye(5);
		BOOST_TEST( arr.size() == 5 );
		BOOST_TEST( arr[2][2] == 1 );
	}
	{
		using multi::elementwise::zeros;
		auto arr = +zeros<2>({3, 4});
		BOOST_TEST( arr.size() == 3 );
		BOOST_TEST( arr[2][2] == 0 );
	}
	// clang-format off
	{
		using multi::elementwise::operator+;  // cppcheck-suppress constStatement ;

		using T = multi::array<int, 2>;
		multi::array<T, 2> XT({2, 3}, T({{1, 0}, {0, 1}}));
		multi::array<T, 2> YT({2, 3}, T({{1, 0}, {0, 1}}));

		auto xy = XT[0][0] + YT[0][0];  // ok
		auto XY = XT + YT;  // error
	}
	// clang-format on

	// {
	// 	multi::array<int, 1> const a = {1, 2, 3};

	// 	using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

	// 	auto const& c = a + ([]() { return 1; } ^ multi::extents_t<0>{});

	// 	BOOST_TEST(( multi::array{2, 3, 4} == c ));
	// 	BOOST_TEST(( c == multi::array{2, 3, 4} ));
	// }
	// {
	// 	multi::array<int, 1> const a = {1, 2, 3};

	// 	using multi::elementwise::operator+;  // cppcheck-suppress constStatement;

	// 	auto const& c = a + 1;

	// 	BOOST_TEST(( multi::array{2, 3, 4} == c ));
	// 	BOOST_TEST(( c == multi::array{2, 3, 4} ));
	// }

	// {
	// 	multi::array const A = {0, 1, 2};        // NOLINT(llvm-header-guard)
	// 	multi::array const B = {10, 11, 12};     // NOLINT(llvm-header-guard)
	// 	multi::array const C = {100, 111, 222};  // NOLINT(llvm-header-guard)

	// 	// multi::array<int, 1> B = {10, 11, 12};

	// 	using multi::elementwise::operator+;     // cppcheck-suppress constStatement;
	// 	using multi::elementwise::operator*;     // cppcheck-suppress constStatement;
	// 	multi::array const D = A + B + 2 * C;  // NOLINT(llvm-header-guard)

	// 	BOOST_TEST( D[2] == A[2] + B[2] + (2 * C[2]) );
	// }
	// np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + np.array([10, 20, 30])
	// array([[11, 22, 33],
	//        [14, 25, 36],
	//        [17, 28, 39]])

	return boost::report_errors();
}
// NOLINTEND(readability-identifier-length)
