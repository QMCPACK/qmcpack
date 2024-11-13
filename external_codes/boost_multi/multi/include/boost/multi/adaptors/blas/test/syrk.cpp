// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/filling.hpp>     // for filling
#include <boost/multi/adaptors/blas/operations.hpp>  // for transposed, T
#include <boost/multi/adaptors/blas/syrk.hpp>        // for syrk
// IWYU pragma: no_include "boost/multi/adaptors/blas/traits.hpp"  // for blas, multi

#include <boost/multi/array.hpp>                     // for array, layout_t

#include <complex>  // for operator*, complex

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_blas_syrk_real) {
		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0},
		};
		{
			multi::array<double, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)
			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::lower, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[2][1] ==   19.0 );
			BOOST_TEST( c[1][2] == 9999.0 );
		}
		{
			multi::array<double, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)
			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::upper, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[1][2] ==   19.0 );
			BOOST_TEST( c[2][1] == 9999.0 );
		}
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)
			namespace blas = multi::blas;

			using blas::filling;
			using blas::syrk;

			syrk(filling::lower, 1.0, a, 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[1][0] ==   34.0 );
			BOOST_TEST( c[0][1] == 9999.0 );
		}
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::upper, 1.0, a, 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, a⸆a, `c` in lower triangular

			BOOST_TEST( c[0][1] ==   34.0 );
			BOOST_TEST( c[1][0] == 9999.0 );
		}
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::upper, 1.0, a, 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, a⸆a, `c` in lower triangular

			BOOST_TEST( c[0][1] ==   34.0 );
			BOOST_TEST( c[1][0] == 9999.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_real_special_case) {
		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
		};
		{
			multi::array<double, 2> c({1, 1}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;
			using blas::filling;

			syrk(filling::lower, 1.0, a, 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[0][0] == 1.0*1.0 + 3.0*3.0 + 4.0*4.0 );
		}
		{
			multi::array<double, 2> c({1, 1}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;
			using blas::filling;

			syrk(filling::upper, 1.0, a, 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[0][0] == 1.0*1.0 + 3.0*3.0 + 4.0*4.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_complex_real_case) {
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<complex, 2> const a = {
			{1.0 + I * 0.0, 3.0 + I * 0.0, 4.0 + I * 0.0},
			{9.0 + I * 0.0, 7.0 + I * 0.0, 1.0 + I * 0.0},
		};
		{
			multi::array<complex, 2> c({3, 3}, 9999.0 + I * 0.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::lower, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( real(c[2][1]) ==   19.0 );
			BOOST_TEST( real(c[1][2]) == 9999.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_complex) {
		using complex = std::complex<double>;

		auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<complex, 2> const a = {
			{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
			{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I},
		};
		{
			multi::array<complex, 2> c({3, 3}, 9999.0 + I * 0.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			syrk(blas::filling::lower, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( real(c[2][1]) == - 3.0 );
			BOOST_TEST( imag(c[2][1]) == -34.0 );
		}
		{
			multi::array<complex, 2> c({2, 2}, 9999.0 + I * 0.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			syrk(blas::filling::lower, 1.0, a, 0.0, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( c[1][0] == complex(18.0, -21.0) );
			BOOST_TEST( c[0][1] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({2, 2}, 9999.0 + I * 0.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			syrk(blas::filling::upper, 1.0, a, 0.0, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( c[0][1] == complex(18.0, -21.0) );
			BOOST_TEST( c[1][0] == 9999.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_operation_complex) {
		using complex = std::complex<double>;

		auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<complex, 2> const a = {
			{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
			{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I},
		};

		{
			multi::array<complex, 2> c({2, 2}, 9999.0 + I * 0.0);  // NOLINT(readability-identifier-length)

			using multi::blas::filling;

			syrk(filling::lower, 1.0, a, 0.0, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( c[1][0] == complex(18.0, -21.0) );
			BOOST_TEST( c[0][1] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)  // NOLINT(fuchsia-default-arguments-calls)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::lower, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(aa⸆)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( c[2][1] == complex(-3.0, -34.0) );
			BOOST_TEST( c[1][2] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)  // NOLINT(fuchsia-default-arguments-calls)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::lower, 1.0, a.rotated(), 0.0, c);  // c⸆=c=a⸆a=(aa⸆)⸆, `c` in lower triangular  // NOLINT(fuchsia-default-arguments-calls)

			BOOST_TEST( c[2][1] == complex(-3.0, -34.0) );
			BOOST_TEST( c[1][2] == 9999.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_operation_real) {
		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0},
		};
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			using multi::blas::filling;

			syrk(filling::lower, 1.0, a, 0.0, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular

			BOOST_TEST( c[1][0] == 34.0 );
			BOOST_TEST( c[0][1] == 9999.0 );
		}
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			using multi::blas::filling;

			syrk(filling::upper, 1.0, a, 0.0, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular

			BOOST_TEST( c[0][1] == 34.0 );
			BOOST_TEST( c[1][0] == 9999.0 );
		}
		{
			multi::array<double, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)

			using multi::blas::filling;

			syrk(filling::lower, 1.0, a.rotated(), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[2][1] == 19.0 );
			BOOST_TEST( c[1][2] == 9999.0 );
		}
		{
			multi::array<double, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::lower, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular

			BOOST_TEST( c[2][1] == 19.0 );
			BOOST_TEST( c[1][2] == 9999.0 );
		}
		{
			multi::array<double, 2> c({3, 3}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			using blas::filling;

			syrk(filling::upper, 1.0, blas::T(a), 0.0, c);  // c⸆=c=a⸆a=(a⸆a)⸆, `c` in upper triangular

			BOOST_TEST( c[1][2] == 19.0 );
			BOOST_TEST( c[2][1] == 9999.0 );
		}
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			namespace blas = multi::blas;

			syrk(blas::filling::upper, 1.0, a, 0.0, blas::T(c));  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular

			BOOST_TEST( c[0][1] == 9999.0 );
			BOOST_TEST( c[1][0] == 34.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_implicit_zero) {
		// NOLINTNEXTLINE(readability-identifier-length)
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0},
		};
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length)

			using multi::blas::filling;

			syrk(filling::lower, 1.0, a, c);  // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular

			BOOST_TEST( c[1][0] == 34.0 );
			BOOST_TEST( c[0][1] == 9999.0 );
		}
	}
	return boost::report_errors();
}
