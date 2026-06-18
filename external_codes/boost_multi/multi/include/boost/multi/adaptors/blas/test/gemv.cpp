// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/axpy.hpp>  // for operator-
#include <boost/multi/adaptors/blas/core.hpp>  // for gemv, context, dot, nrm2
#include <boost/multi/adaptors/blas/dot.hpp>   // for dot, dot_ref
#include <boost/multi/adaptors/blas/gemv.hpp>  // for gemv_range, gemv, oper...
#include <boost/multi/adaptors/blas/nrm2.hpp>  // for operator^
#include <boost/multi/array.hpp>               // for array, layout_t, array...
#include <boost/multi/elementwise.hpp>         // for operations
#include <boost/multi/restriction.hpp>         // for restriction

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for generate, transform
#include <cmath>      // for abs  // IWYU pragma: keep
#include <complex>    // for complex, operator*
// #include <cstdlib>    // IWYU pragma: keep
// IWYU pragma: no_include <cstdlib>  // for abs
#include <iterator>  // for size, begin
// IWYU pragma: no_include <memory>       // for allocator
#include <numeric>      // for inner_product
#include <random>       // for normal_distribution
#include <type_traits>  // for is_same_v
#include <utility>      // for move, forward
// IWYU pragma: no_include <stdlib.h>       // for abs

namespace multi = boost::multi;
namespace blas  = multi::blas;

// using fp_types = boost::mpl::list<double, float>;  // old versions of Boost.Test need MPL Type lists explicitly

namespace {
template<class M, class VI, class VO>
auto MV(M const& a, VI const& x, VO&& y) -> VO&& {  // NOLINT(readability-identifier-naming,readability-identifier-length) BLAS naming
	std::transform(
		a.begin(), a.end(), y.begin(),
		[&x](auto const& row) { return std::inner_product(row.begin(), row.end(), x.begin(), typename VI::value_type{}); }
	);
	return std::forward<VO>(y);
}
}  // end unnamed namespace

namespace {
void gemv_broadcast() {
	// NOLINTNEXTLINE(readability-identifier-length)
	multi::array<double, 2> const a = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0}
	};

	{
		multi::array<double, 1> const ones({3}, 1.0);

		BOOST_TEST( ones.stride() == 1 );

		BOOST_TEST( ones[0] == 1.0 );
		BOOST_TEST( ones[1] == 1.0 );
		BOOST_TEST( ones[2] == 1.0 );

		multi::array<double, 1> sum_by_rows({2}, 0.0);
		blas::gemv_n(1.0, a.begin(), 2, ones.begin(), 0.0, sum_by_rows.begin());

		BOOST_TEST( std::abs( sum_by_rows[0] - (1.0 + 2.0 + 3.0)) < 1.0e-8 );
		BOOST_TEST( std::abs( sum_by_rows[1] - (4.0 + 5.0 + 6.0)) < 1.0e-8 );
	}
	// BLAS GEMV doesn't work with stride zero
	// {
	//  multi::array<double, 0> const one(1.0);
	//  auto const& ones = one.broadcasted();

	//  BOOST_TEST( ones.stride() == 0 );

	//  BOOST_TEST( ones[0] == 1.0 );
	//  BOOST_TEST( ones[1] == 1.0 );
	//  BOOST_TEST( ones[2] == 1.0 );

	//  multi::array<double, 1> sum_by_rows({2}, 0.0);
	//  blas::gemv_n(1.0, a.begin(), 2, ones.begin(), 0.0, sum_by_rows.begin());

	//  std::cout << sum_by_rows[0] << " " << sum_by_rows[1] << "\n";
	//  BOOST_TEST( std::abs( sum_by_rows[0] - (1.0 + 2.0 + 3.0)) < 1.0e-8 );
	//  BOOST_TEST( std::abs( sum_by_rows[1] - (4.0 + 5.0 + 6.0)) < 1.0e-8 );
	// }
}

}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_double)
	{
		using T = double;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<T, 2> const a = {
			{ 9.0, 24.0, 30.0, 9.0},
			{ 4.0, 10.0, 12.0, 7.0},
			{14.0, 16.0, 36.0, 1.0},
		};
		multi::array<T, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemv_n(1.0, a.begin(), a.size(), x.begin(), 0.0, y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3 ) < 0.0001);
			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( y[2] -  +blas::dot(a[2], x)) < 0.0001);
			}
		}
		{
			multi::array<T, 1>       y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 2> const aT{~a};
			blas::gemv_n(1.0, (~aT).begin(), (~aT).size(), x.begin(), 0.0, y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3 ) < 0.0001);

			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( y[2] - +blas::dot(a[2], x)) < 0.0001);
			}
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			auto               mv = blas::gemv(1.0, a, x);
			copy_n(mv.begin(), mv.size(), y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001);

			multi::array<T, 1> w2(multi::extents_t<1>{multi::iextension{size(a)}});
			MV(a, x, w2);
			BOOST_TEST( std::abs(w2[0] - y[0]) < 0.00001);
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			y = blas::gemv(1.0, a, x);
			BOOST_TEST( std::abs(y[1] - 91.3) < 0.00001);
		}
		{
			multi::array<T, 1> y = blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( std::abs(y[1] - 91.3) < 0.00001);
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}}, 0.);  // NOLINT(readability-identifier-length) BLAS naming
			y += blas::gemv(1.0, a, x);
			BOOST_TEST( std::abs(y[1] - 91.3) < 0.00001);
		}
		{
			multi::array<T, 1> y = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemv(1.1, a, x, 1.0, y);           // y = a*M*x + b*y
			BOOST_TEST( std::abs(y[1] - 105.43) < 0.00001);
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_float)
	{
		using T = float;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<T, 2> const a = {
			{ 9.0, 24.0, 30.0, 9.0},
			{ 4.0, 10.0, 12.0, 7.0},
			{14.0, 16.0, 36.0, 1.0},
		};
		multi::array<T, 1> const x = {1.1F, 2.1F, 3.1F, 4.1F};  // NOLINT(readability-identifier-length) BLAS naming
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemv_n(1.0, a.begin(), a.size(), x.begin(), 0.0, y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3F ) < 0.0001F );
			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( y[2] -  +blas::dot(a[2], x)) < 0.0001F );
			}
		}
		{
			multi::array<T, 1>       y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 2> const aT{~a};
			blas::gemv_n(1.0, (~aT).begin(), (~aT).size(), x.begin(), 0.0, y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3F ) < 0.0001F );

			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( y[2] - +blas::dot(a[2], x)) < 0.0001F );
			}
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			auto               mv = blas::gemv(1.0, a, x);
			copy_n(mv.begin(), mv.size(), y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3F) < 0.00001F);

			multi::array<T, 1> w2(multi::extents_t<1>{multi::iextension{size(a)}});
			MV(a, x, w2);
			BOOST_TEST( std::abs(w2[0] - y[0]) < 0.0001F );
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) BLAS naming
			y = blas::gemv(1.0, a, x);
			BOOST_TEST( std::abs(y[1] - 91.3F) < 0.00001F );
		}
		{
			multi::array<T, 1> y = blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( std::abs(y[1] - 91.3F) < 0.00001F);
		}
		{
			multi::array<T, 1> y(multi::extents_t<1>{multi::iextension{size(a)}}, 0.);  // NOLINT(readability-identifier-length) BLAS naming
			y += blas::gemv(1.F, a, x);
			BOOST_TEST( std::abs(y[1] - 91.3F) < 0.00001F );
		}
		{
			multi::array<T, 1> y = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemv(1.1F, a, x, 1.0F, y);         // y = a*M*x + b*y
			BOOST_TEST( std::abs( y[1] - 105.43F) < 0.00001F );
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_double)
	{
		using T        = double;
		namespace blas = multi::blas;

		using std::abs;  // NOLINT(misc-include-cleaner) bug in clang-tidy 21.1.2
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<T, 2> const a = {
			{ 9.0, 24.0, 30.0, 9.0},
			{ 4.0, 10.0, 12.0, 7.0},
			{14.0, 16.0, 36.0, 1.0},
		};
		multi::array<T, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming
		{
			multi::array<T, 1> y     = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			T const            alpha = 1.1;
			T const            beta  = 1.2;
			blas::gemv(alpha, a, x, beta, y);  // y = a*M*x + b*y

			multi::array<T, 1> const y3 = {214.02, 106.43, 188.37};
			BOOST_TEST( std::abs(y[1] - y3[1]) < 2e-14 );
		}
		if constexpr(!std::is_same_v<T, float>) {
			auto Y = +blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( std::abs( Y[0] - +blas::dot(a[0], x)) < 0.00001 );
			BOOST_TEST( std::abs( Y[1] - +blas::dot(a[1], x)) < 0.00001);
			BOOST_TEST( std::abs( Y[2] - +blas::dot(a[2], x)) < 0.00001);
		}
		{
			multi::array<T, 1> const x_shadow = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 1> const y        = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 1> const dot      = blas::gemv(1.0, multi::array<T, 2>({x_shadow}), y);
			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( dot[0] - blas::dot(x_shadow, y) ) < 1e-10 );
			}
		}
		{
			using blas::operators::operator%;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			using blas::operators::operator-;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			using blas::operators::operator^;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			BOOST_TEST((((~+~a) % x - a % x) ^ 2) < 1e-9);
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_float)
	{
		using T        = float;
		namespace blas = multi::blas;

		using std::abs;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<T, 2> const a = {
			{ 9.0, 24.0, 30.0, 9.0},
			{ 4.0, 10.0, 12.0, 7.0},
			{14.0, 16.0, 36.0, 1.0},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<T, 1> const x = {static_cast<T>(1.1), static_cast<T>(2.1), static_cast<T>(3.1), static_cast<T>(4.1)};
		{
			multi::array<T, 1> y     = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			auto const         alpha = static_cast<T>(1.1);
			auto const         beta  = static_cast<T>(1.2);
			blas::gemv(alpha, a, x, beta, y);  // y = a*M*x + b*y

			multi::array<T, 1> const y3 = {static_cast<T>(214.02), static_cast<T>(106.43F), static_cast<T>(188.37)};
			BOOST_TEST( std::abs(y[1] - y3[1]) < 2e-14F );
		}
		if constexpr(!std::is_same_v<T, float>) {
			auto Y = +blas::gemv(1.0, a, x);  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( std::abs( Y[0] - +blas::dot(a[0], x)) < 0.0001F );
			BOOST_TEST( std::abs( Y[1] - +blas::dot(a[1], x)) < 0.0001F );
			BOOST_TEST( std::abs( Y[2] - +blas::dot(a[2], x)) < 0.0001F );
		}
		{
			multi::array<T, 1> const x_shadow = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 1> const y        = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) BLAS naming
			multi::array<T, 1> const dot      = blas::gemv(1.0, multi::array<T, 2>({x_shadow}), y);
			if(!std::is_same_v<T, float>) {  // workaround Apple Accelerate BLAS bug in dot
				BOOST_TEST( std::abs( dot[0] - +blas::dot(x_shadow, y) ) < 1e-10F );
			}
		}
		{
			using blas::operators::operator%;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			using blas::operators::operator-;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			using blas::operators::operator^;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			BOOST_TEST((((~+~a) % x - a % x) ^ 2) < 1e-9);
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_complex_double)
	{
		using T        = double;
		namespace blas = multi::blas;
		using complex  = std::complex<T>;
		using std::abs;

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const M = {
			{ {9.0, 0.0}, {24.0, 0.0}, {30.0, 0.0}, {9.0, 0.0}},
			{ {4.0, 0.0}, {10.0, 0.0}, {12.0, 0.0}, {7.0, 0.0}},
			{{14.0, 0.0}, {16.0, 0.0}, {36.0, 0.0}, {1.0, 0.0}},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> const X = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0},
			{4.1, 0.0},
		};
		{
			// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
			multi::array<complex, 1> Y = {
				{4.0, 0.0},
				{5.0, 0.0},
				{6.0, 0.0},
			};

			auto const alpha = T{1.1};
			auto const beta  = T{1.2};

			blas::gemv(alpha, M, X, beta, Y);  // y = a*M*x + b*y

			multi::array<complex, 1> const Y3 = {
				{214.02, 0.0},
				{106.43, 0.0},
				{188.37, 0.0},
			};

			using blas::operators::operator-;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			T const n2{blas::nrm2(Y - Y3)};
			BOOST_TEST(n2 < T{1.0e-4});
		}
	}

	// BOOST_AUTO_TEST_CASE(blas_gemv_complex_float_mimic_cdot)
	{
		multi::array<std::complex<float>, 2> const v1 = {
			{std::complex<float>{1.0F, 2.0F}, std::complex<float>{3.0F, 4.0F}, std::complex<float>{5.0F, 6.0F}}
		};
		BOOST_TEST( v1.size() == 1 );
		BOOST_TEST( v1.num_elements() == 3 );

		multi::array<std::complex<float>, 1> const v2 = {
			std::complex<float>{ 7.0F,  8.0F},
			std::complex<float>{ 9.0F, 10.0F},
			std::complex<float>{11.0F, 12.0F}
		};
		BOOST_TEST( v2.size() == 3 );

		multi::array<std::complex<float>, 1> res({1}, std::complex<float>{});
		BOOST_TEST( res.size() == 1 );

		blas::gemv(1.F, v1, v2, 0.F, res);

		BOOST_TEST( std::abs(res[0] - (v1[0][0]*v2[0] + v1[0][1]*v2[1] + v1[0][2]*v2[2])) < 1e-8F );

		std::complex<float> res_dot;

		blas::dot(v1[0], v2, res_dot);

		BOOST_TEST( std::abs(res[0] - res_dot) < 1e-8F );
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_complex_float)
	{
		using T        = float;
		namespace blas = multi::blas;
		using complex  = std::complex<T>;
		using std::abs;

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const M = {
			{ {9.0, 0.0}, {24.0, 0.0}, {30.0, 0.0}, {9.0, 0.0}},
			{ {4.0, 0.0}, {10.0, 0.0}, {12.0, 0.0}, {7.0, 0.0}},
			{{14.0, 0.0}, {16.0, 0.0}, {36.0, 0.0}, {1.0, 0.0}},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> const X = {
			{1.1F, 0.0F},
			{2.1F, 0.0F},
			{3.1F, 0.0F},
			{4.1F, 0.0F},
		};
		{
			// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
			multi::array<complex, 1> Y = {
				{4.0, 0.0},
				{5.0, 0.0},
				{6.0, 0.0},
			};

			auto const alpha = 1.1F;
			auto const beta  = 1.2F;

			blas::gemv(alpha, M, X, beta, Y);  // y = a*M*x + b*y

			multi::array<complex, 1> const Y3 = {
				{214.02F, 0.0F},
				{106.43F, 0.0F},
				{188.37F, 0.0F},
			};

			using blas::operators::operator-;  // cppcheck-suppress constStatement ; bug in cppcheck 2.18
			T const n2{blas::nrm2(Y - Y3)};
			BOOST_TEST( std::abs(n2) < 1e-4F );
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex)
	{
		namespace blas = multi::blas;
		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		using std::abs;

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const a = {
			{2.0 + 3.0 * I, 2.0 + 1.0 * I, 1.0 + 2.0 * I},
			{4.0 + 2.0 * I, 2.0 + 4.0 * I, 3.0 + 1.0 * I},
			{7.0 + 1.0 * I, 1.0 + 5.0 * I, 0.0 + 3.0 * I},
		};
		multi::array<complex, 1> const x = {1.0 + 2.0 * I, 2.0 + 1.0 * I, 9.0 + 2.0 * I};  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST(( +blas::gemv(1., a, x) == multi::array<complex, 1>{4.0 + 31.*I, 25.0 + 35.0*I, -4.0 + 53.0*I} ));

		auto aT = +~a;
		BOOST_TEST(( +blas::gemv(1., ~aT, x) == multi::array<complex, 1>{4.0 + 31.0*I, 25.0 + 35.0*I, -4.0 + 53.0*I} ));

		BOOST_TEST( +blas::gemv(1., ~a, x) == (multi::array<complex, 1>{63.0 + 38.0*I, -1.0 + 62.0*I, -4.0 + 36.0*I}) );
		BOOST_TEST( +blas::gemv(1., ~a, x) == + blas::gemv(1.0, aT, x) );
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_temporary)
	{
		using complex = std::complex<double>;

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const A = {
			{{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
			{{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}},
			{{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		auto const B = [](auto array) {
			// NOLINTNEXTLINE(cert-msc32-c,cert-msc51-cpp,bugprone-random-generator-seed) test purposes
			auto rand = [gauss = std::normal_distribution<>{}, gen = std::mt19937{}]() mutable {
				return complex{gauss(gen), gauss(gen)};
			};
			std::generate(array.elements().begin(), array.elements().end(), rand);
			return array;
		}(multi::array<complex, 2>({3, 3}));

		// using blas::operators::operator*;
		// using blas::operators::operator-;
		// using blas::operators::operator^;
		// BOOST_TEST( (((+(A*B))[0] - B[0])^2) == 0.0 );
		// BOOST_TEST( (((+(A*B))[1] - B[1])^2) == 0.0 );
		// BOOST_TEST( (((+(A*B))[2] - B[2])^2) == 0.0 );
	}

	// BOOST_AUTO_TEST_CASE(multi_blas_gemv_context)
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const a = {
			{ 9.0, 24.0, 30.0, 9.0},
			{ 4.0, 10.0, 12.0, 7.0},
			{14.0, 16.0, 36.0, 1.0},
		};
		multi::array<double, 1> const x = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) conventional name in BLAS

		blas::context ctxt;
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			blas::gemv_n(&ctxt, 1.0, begin(a), size(a), begin(x), 0.0, begin(y));
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.0001);
			BOOST_TEST( std::abs( y[2] - +blas::dot(a[2], x)) < 0.0001);
		}
		{
			multi::array<double, 1>       y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			multi::array<double, 2> const aT{~a};
			blas::gemv_n(&ctxt, 1.0, begin(~aT), size(~aT), begin(x), 0.0, begin(y));
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001);
			BOOST_TEST( std::abs( y[2] - +blas::dot(a[2], x)) < 0.00001);
		}
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			auto&&                  mv = blas::gemv(&ctxt, 1.0, a, x);
			copy_n(mv.begin(), mv.size(), y.begin());
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			y = blas::gemv(&ctxt, 1.0, a, x);
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			y = blas::gemv(1.0, a, x);
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}});  // NOLINT(readability-identifier-length) conventional name in BLAS
			y() = blas::gemv(1.0, a, x);
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y = blas::gemv(&ctxt, 1.0, a, x);  // NOLINT(readability-identifier-length) conventional name in BLAS
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y(multi::extents_t<1>{multi::iextension{size(a)}}, 0.0);  // NOLINT(readability-identifier-length) conventional name in BLAS
			y += blas::gemv(&ctxt, 1.0, a, x);
			BOOST_TEST( std::abs( y[1] - 91.3) < 0.00001 );
		}
		{
			multi::array<double, 1> y = {4.0, 5.0, 6.0};  // NOLINT(readability-identifier-length) conventional name in BLAS
			y += blas::gemv(&ctxt, 1.1, a, x);
			BOOST_TEST( std::abs( y[1] - 105.43) < 0.00001 );
		}
	}

	gemv_broadcast();

#ifndef __NVCC__
	{
		multi::array<double, 2> const arr = {
			{1.0, 2.0},
			{3.0, 4.0}
		};
		multi::array<double, 1> const vec = {1.0, 2.0};

		multi::array<double, 1> vec2 = {0.0, 0.0};

		vec2 = multi::blas::gemv(5.0, arr, vec);

		BOOST_TEST( vec2[1] - 55.0 < 1e-7 );

		using multi::elementwise::operator+;  // cppcheck-suppress constStatement; bug in v2.19.0
		using multi::blas::gemv;
		using multi::elementwise::exp;
		using multi::elementwise::log;

		auto ret = log(+gemv(5.0, arr, vec) + exp(vec));

		BOOST_TEST( std::abs( ret[1] - 4.13339 ) < 1e-4 );
	}
#endif

	return boost::report_errors();
}
