// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/blas/core.hpp>        // for context
#include <boost/multi/adaptors/blas/gemm.hpp>        // for gemm, gemm_range
// IWYU pragma: no_include "boost/multi/adaptors/blas/numeric.hpp"     // for involuter, conju...
#include <boost/multi/adaptors/blas/operations.hpp>  // for H, T, (anonymous)

#include <boost/multi/array.hpp>                     // for layout_t, array

#include <cmath>  // for abs  // IWYU pragma: keep
#include <complex>   // for complex, operator*
#include <iterator>  // for begin, size
// IWYU pragma: no_include <memory>

namespace multi = boost::multi;
namespace blas  = multi::blas;

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/
#define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )
// #define BOOST_REQUIRE_SMALL(X, ToL) BOOST_TEST( std::abs( X ) < (ToL) )

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_T_sub) {
		namespace blas = multi::blas;

		multi::array<double, 2> A({100, 4}, 1.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> B({4, 4}, 1.0);    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({100, 1}, 0.0);  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(1.0, A({0, 100}, {1, 2}), blas::T(B)({0, 1}, {0, 1}), 0.0, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_H_sub) {
		multi::array<double, 2> A({100, 4}, 1.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> B({4, 4}, 1.0);    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({100, 1}, 0.0);  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(1., A({0, 100}, {1, 2}), blas::H(B)({0, 1}, {0, 1}), 0.0, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_H_sub_6) {
		multi::array<double, 2> A({100, 4}, 2.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> B({4, 4}, 3.0);    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({100, 1}, 0.0);  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(1., A({0, 100}, {1, 2}), blas::H(B)({0, 1}, {0, 1}), 0.0, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 6.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_H_copy) {
		multi::array<double, 2> A({100, 4}, 1.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> B({4, 4}, 1.0);    // NOLINT(readability-identifier-length) BLAS naming

		auto C = +blas::gemm(1., A({0, 100}, {1, 2}), blas::H(B)({2, 3}, {2, 3}));  // c=ab, c⸆=b⸆a⸆  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_complex_100x1_1x1) {
		using complex = std::complex<double>;
		multi::array<complex, 2> const A({100, 1}, {1.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const B({1, 1}, {1.0, 0.0});    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<complex, 2> C({100, 1}, {0.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm({1.0, 0.0}, A, B, {0.0, 0.0}, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_complex_100x1_1x1_T) {
		using complex = std::complex<double>;
		multi::array<complex, 2> const A({100, 1}, complex{1.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const B({1, 1}, complex{1.0, 0.0});    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<complex, 2> C({100, 1}, complex{0.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(complex{1.0, 0.0}, A, blas::T(B), complex{0.0, 0.0}, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( C[99][0] == 1.0 );
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_complex_100x1_1x1_H) {
		using complex = std::complex<double>;                    // complex const I{0, 1};
		multi::array<complex, 2> const A({100, 1}, {1.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const B({1, 1}, {1.0, 0.0});    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<complex, 2> C({100, 1}, {0.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm({1.0, 0.0}, A, blas::H(B), {0.0, 0.0}, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST( C[99][0] == 1.0 );
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1) {
		using complex = std::complex<double>;
		multi::array<complex, 2> const A({100, 1}, {1.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const B({1, 1}, {1.0, 0.0});    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<complex, 2> C({100, 1}, {0.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm({1.0, 0.0}, A, B, {0.0, 0.0}, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_T) {
		multi::array<double, 2> const A({100, 1}, 1.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const B({1, 1}, 1.0);    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({100, 1}, 0.0);  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(1.0, A, blas::T(B), 0.0, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(adaptor_blas_double_100x1_1x1_H) {
		multi::array<double, 2> const A({100, 1}, 1.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const B({1, 1}, 1.0);    // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({100, 1}, 0.0);  // NOLINT(readability-identifier-length) BLAS naming

		blas::gemm(1.0, A, blas::H(B), 0.0, C);  // c=ab, c⸆=b⸆a⸆
		BOOST_TEST(C[99][0] == 1.0);
	}

	BOOST_AUTO_TEST_CASE(multi_blas_gemm_square_real) {
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0},
			{1.0, 2.0, 3.0},
		};
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const b = {
			{11.0, 12.0, 4.0},
			{ 7.0, 19.0, 1.0},
			{11.0, 12.0, 4.0},
		};
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, b, 0.0, c);
			BOOST_TEST( c[2][1] == 86.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			BOOST_TEST( size( a) == size( c) );
			BOOST_TEST( size(~b) == size(~c) );
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 86.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, blas::T(b), 0.0, c);
			BOOST_TEST( c[2][1] == 48.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1., a.begin(), a.size(), blas::T(b).begin(), 0.0, c.begin());
			BOOST_TEST( c[2][1] == 48.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, blas::T(a), b, 0.0, c);
			BOOST_TEST( c[2][1] == 103.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(blas::T(a)), size(blas::T(a)), begin(b), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 103.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, blas::T(a), blas::T(b), 0.0, c);
			BOOST_TEST( c[2][1] == 50.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 50.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, blas::T(b), 0.0, c);
			BOOST_TEST( c[2][1] == 48.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(a), size(a), begin(blas::T(b)), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 48.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, blas::T(a), b, 0.0, c);
			BOOST_TEST( c[2][1] == 103.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)}, 9999.0);  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(blas::T(a)), size(blas::T(a)), begin(b), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 103.0 );
		}
		{
			multi::array<double, 2> c({a.size(), b.rotated().size()}, 9999.0);  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemm(2.0, blas::H(a), blas::H(b), 0.0, c);
			BOOST_TEST( c[2][1] == 100.0 );
		}
		{
			multi::array<double, 2> c = blas::gemm(2.0, blas::H(a), blas::H(b));  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( c[2][1] == 100.0 );
		}
		{
			multi::array<double, 2> const c = blas::gemm(2.0, blas::H(a), blas::H(b));  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( c[2][1] == 100.0 );
		}
		{
			multi::array<double, 2> c({a.size(), b.rotated().size()}, 9999.0);  // NOLINT(readability-identifier-length) BLAS naming
			c = blas::gemm(2.0, blas::H(a), blas::H(b));
			BOOST_TEST( c[2][1] == 100.0 );
		}
		{
			multi::array<double, 2> c;  // NOLINT(readability-identifier-length) BLAS naming
			c = blas::gemm(2.0, blas::H(a), blas::H(b));
			BOOST_TEST( c[2][1] == 100.0 );
		}
		{
			multi::array<double, 2> c({a.size(), b.rotated().size()}, 9999.0);  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemm_n(2.0, begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0.0, begin(c));
			BOOST_TEST( c[2][1] == 100.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square) {
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const a = {
			{1.0, 3.0},
			{9.0, 7.0},
		};
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
		};
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, b, 0.0, c);      // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 148.0 );
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming

			blas::context const ctxt;
			blas::gemm_n(&ctxt, 1.0, begin(a), size(a), begin(b), 0.0, begin(c));
			BOOST_TEST( c[1][0] == 148.0 );
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, ~a, b, 0.0, c);     // c=a⸆b, c⸆=b⸆a
			BOOST_TEST(( c[1][1] == 169.0 && c[1][0] == 82.0 ));
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming

			blas::context const ctxt;
			blas::gemm_n(&ctxt, 1.0, begin(~a), size(~a), begin(b), 0.0, begin(c));
			BOOST_TEST(( c[1][1] == 169 && c[1][0] == 82 ));
		}
		{
			multi::array<double, 2> const c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming

			blas::context const ctxt;
			blas::gemm_n(&ctxt, 1.0, begin(~a), size(~a), begin(b), 0.0, begin(~c));
			BOOST_TEST( (~c)[1][1] == 169 );
			BOOST_TEST( (~c)[1][0] ==  82 );
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, ~b, 0.0, c);     // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == 183.0 );
		}
		{
			// TODO(correaa) fix sfinae of const c
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming

			blas::context const ctxt;
			blas::gemm_n(&ctxt, 1.0, begin(a), size(a), begin(~b), 0.0, begin(c));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == 183.0 );
		}
		{
			// NOLINTNEXTLINE(misc-const-correctness) TODO(correaa) fix sfinae of const c
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, ~b, 0.0, ~c);    // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( (~c)[1][0] == 183.0 );
		}
		{
			// NOLINTNEXTLINE(misc-const-correctness) TODO(correaa) fix sfinae of const c
			multi::array<double, 2> c({2, 2});                                // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(a), size(a), begin(~b), 0.0, begin(~c));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( (~c)[1][0] == 183.0 );
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, ~a, ~b, 0.0, c);    // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == 117.0 );
		}
		{
			multi::array<double, 2> c({2, 2});                                 // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(~a), size(~a), begin(~b), 0.0, begin(c));  // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == 117.0 );
		}
		{
			multi::array<double, 2> c({2, 2});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, ~a, ~b, 0.0, ~c);   // c⸆=a⸆b⸆, c=ba
			BOOST_TEST( c[0][1] == 117.0 );
		}
		{
			multi::array<double, 2> c({2, 2});                                  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(~a), size(~a), begin(~b), 0.0, begin(~c));  // c⸆=a⸆b⸆, c=ba
			BOOST_TEST( c[0][1] == 117.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare) {
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const a = {
			{1.0, 3.0, 1.0},
			{9.0, 7.0, 1.0},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const b = {
			{11.0, 12.0, 1.0},
			{ 7.0, 19.0, 1.0},
			{ 1.0,  1.0, 1.0},
		};
		{
			multi::array<double, 2> c({2, 3});  // NOLINT(readability-identifier-length) BLAS naming
			blas::gemm(1.0, a, b, 0.0, c);      // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 17 );
		}
		{
			multi::array<double, 2> c({2, 3});                              // NOLINT(readability-identifier-length) BLAS naming
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 17.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic) {
		namespace blas = multi::blas;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const a = {
			{1.0, 3.0, 1.0},
			{9.0, 7.0, 1.0},
		};
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<double, 2> const b = {
			{11.0, 12.0, 4.0, 8.0},
			{ 7.0, 19.0, 2.0, 7.0},
			{ 5.0,  3.0, 3.0, 1.0},
		};
		{
			multi::array<double, 2> c({size(a), size(~b)});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(1.0, a, b, 0.0, c);                   // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 53.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)});                 // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 53.0 );
		}
		{
			multi::array<double, 2> c({2, 4});  // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm(0.1, a, b, 0.0, c);      // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
		{
			multi::array<double, 2> c({2, 4});                             // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n(0.1, begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
		{
			auto c = +blas::gemm(0.1, a, b);  // c=ab, c⸆=b⸆a⸆  // NOLINT(readability-identifier-length) conventional BLAS naming
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
		{
			multi::array c = blas::gemm(0.1, a, b);  // NOLINT(readability-identifier-length) conventional BLAS naming
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const a = {
			{1.0 - 2.0 * I, 9.0 - 1.0 * I},
			{2.0 + 3.0 * I, 1.0 - 2.0 * I},
		};
		{
			auto c = +blas::gemm(1.0, a, blas::H(a));  // c=aa†, c†=aa†  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array c = blas::gemm(1.0, a, blas::H(a));  // c=aa†, c†=aa†  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( c[1][0] == 7.-10.*I );
			BOOST_TEST( c[0][1] == 7.+10.*I );
		}
		{
			multi::array<complex, 2> c = blas::gemm(1.0, a, blas::H(a));  // c=aa†, c†=aa†  // NOLINT(readability-identifier-length) conventional BLAS naming
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) conventional BLAS naming
			c = blas::gemm(1.0, a, blas::H(a));                 // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) conventional BLAS naming
			c() = blas::gemm(1.0, a, blas::H(a));               // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.*I );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});     // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm({1.0, 0.0}, a, blas::H(a), {0.0, 0.0}, c);  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});                                     // NOLINT(readability-identifier-length) conventional BLAS naming
			blas::gemm_n({1.0, 0.0}, begin(a), size(a), begin(blas::H(a)), {0.0, 0.0}, begin(c));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7. - 10.*I );
			BOOST_TEST( c[0][1] == 7. + 10.*I );
		}
	}

#if defined(CUDA_FOUND)
	#include <thrust/complex.h>
	BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh_thrust) {
		using complex = thrust::complex<double>;
		complex const                  I{0.0, 1.0};
		multi::array<complex, 2> const a = {
			{1.0 - 2.0 * I, 9.0 - 1.0 * I},
			{2.0 + 3.0 * I, 1.0 - 2.0 * I}
		};
		{
			auto c = +blas::gemm(1.0, a, blas::hermitized(a));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array c = blas::gemm(1.0, a, blas::hermitized(a));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c = blas::gemm(1.0, a, blas::hermitized(a));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			c = blas::gemm(1.0, a, blas::hermitized(a));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, blas::hermitized(a), 0.0, c);  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1.0, begin(a), size(a), begin(blas::H(a)), 0.0, begin(c));  // c=aa†, c†=aa†
			BOOST_TEST( c[1][0] == 7.0 - 10.0*I );
			BOOST_TEST( c[0][1] == 7.0 + 10.0*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated) {
		using complex = std::complex<double>;
		complex const                  I{0.0, 1.0};
		multi::array<complex, 2> const a = {
			{1.0 - 2.0 * I, 9.0 - 1.0 * I}
		};
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm(1.0, a, blas::H(a), 0.0, c);  // c=aa†, c†=aa†
			BOOST_TEST( c[0][0] == 87.0 + 0.0*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1.0, begin(a), size(a), begin(blas::H(a)), 0.0, begin(c));  // c=aa†, c†=aa†
			BOOST_TEST( c[0][0] == 87.0 + 0.0*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bisbis) {
		using complex = std::complex<double>;
		complex const                  I{0.0, 1.0};
		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I},
			{9.0 - 1.0 * I},
			{1.0 + 1.0 * I}
		};
		multi::array<complex, 2> const b = {
			{11.0 - 2.0 * I, 7.0 - 3.0 * I, 8.0 - 1.0 * I}
		};
		{
			multi::array<complex, 2> c({1, 1});

			BOOST_TEST( size(blas::H(a)) == 1 );
			BOOST_TEST( size(blas::H(b)[0]) == 1 );

			blas::gemm(1.0, blas::H(a), blas::H(b), 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 84.0 + 7.0*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1.0, begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 84.0 + 7.0*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_empty) {
		multi::array<double, 2> const a({0, 5});
		BOOST_TEST( size( a) == 0 );
		BOOST_TEST( size(~a) == 5 );
		BOOST_TEST( a.is_empty() );

		multi::array<double, 2> const b({5, 0});
		BOOST_TEST( size( b) == 0 );
		BOOST_TEST( size(~b) == 0 );
		BOOST_TEST( b.is_empty() );
		{
			multi::array<double, 2> c;
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
		}
		{
			multi::array<double, 2> c;
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare2) {
		multi::array<double, 2> const a = {
			{1.0, 3.0},
			{9.0, 7.0},
			{1.0, 1.0},
		};
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
		};
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[2][1] == 31.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[2][1] == 31.0 );
		}
		{
			multi::array<double, 2> c({size(~b), size(a)});
			blas::gemm(1.0, a, b, 0.0, ~c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 31.0 );
		}
		{
			multi::array<double, 2> c({size(~b), size(a)});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(~c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 31.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({3, 2});
			blas::gemm(1.0, ~ar, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[2][1] == 31.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({3, 2});
			blas::gemm_n(1.0, begin(~ar), size(~ar), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[2][1] == 31.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({2, 3});
			blas::gemm(1.0, ~ar, b, 0.0, ~c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 31.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({2, 3});
			blas::gemm_n(1.0, begin(~ar), size(~ar), begin(b), 0.0, begin(~c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 31.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x2_2x2) {
		multi::array<double, 2> const a = {
			{1.0, 3.0},
			{9.0, 4.0},
		};
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
		};
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1.0, ~a, b, 0.0, c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 61.0 );

			blas::gemm(1.0, ~a, b, 0.0, ~c);  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[0][1] == 61.0 );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm_n(1.0, begin(~a), size(~a), begin(b), 0.0, begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 61.0 );

			blas::gemm_n(1.0, begin(~a), size(~a), begin(b), 0.0, begin(~c));  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[0][1] == 61.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x2) {
		multi::array<double, 2> const a = {
			{1.0, 3.0},
			{9.0, 4.0},
			{1.0, 5.0},
		};
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
			{ 8.0,  1.0},
		};
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1.0, ~a, b, 0.0, c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 101.0 );

			blas::gemm(1., ~a, b, 0., ~c);  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[0][1] == 101 );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm_n(1.0, begin(~a), size(~a), begin(b), 0.0, begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 101.0 );

			blas::gemm_n(1.0, begin(~a), size(~a), begin(b), 0.0, begin(~c));  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[0][1] == 101.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x2) {
		multi::array<double, 2> const a = {
			{1.0, 9.0, 1.0}
		};
		BOOST_TEST( stride(~a) == 1 );
		BOOST_TEST( stride( a) == 3 );
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
			{ 8.0,  1.0},
		};
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({size(~b), size(~ar)});
			blas::gemm(1.0, ~ar, b, 0.0, ~c);  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
		{
			auto ar = +~a;
			BOOST_TEST( size(~ar) == 1 );
			BOOST_TEST( begin(~ar).stride() == 1 );
			BOOST_TEST( begin(~ar)->stride() == 1 );
			BOOST_TEST( begin( ar)->stride() == 1 );

			multi::array<double, 2> c({size(~b), size(~ar)});
			BOOST_TEST( begin( c).stride() == 1 );
			BOOST_TEST( begin(~c).stride() == 1 );
			BOOST_TEST( begin(c)->stride() == 1 );

			BOOST_TEST( begin(b) );
			blas::gemm_n(1.0, begin(~ar), size(~ar), begin(b), 0.0, begin(~c));  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_3x2) {
		using complex                    = std::complex<double>;
		multi::array<complex, 2> const a = {
			{1.0, 9.0, 1.0}
		};
		BOOST_TEST( stride(~a) == 1 );
		BOOST_TEST( stride( a) == 3 );
		multi::array<complex, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
			{ 8.0,  1.0},
		};
		{
			multi::array<complex, 2> c({size(a), size(~b)});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			multi::array<complex, 2> c({size(a), size(~b)});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({size(~b), size(~ar)});
			blas::gemm(1.0, ~ar, b, 0.0, ~c);  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({size(~b), size(~ar)});
			blas::gemm_n(1.0, begin(~ar), size(~ar), begin(b), 0.0, begin(~c));  // c⸆=a⸆b, c=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_part_3x2) {
		multi::array<double, 2> const a = {
			{1.0, 9.0, 1.0},
			{3.0, 3.0, 3.0},
		};
		BOOST_TEST( stride(~a) == 1 );
		BOOST_TEST( stride( a) == 3 );
		multi::array<double, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
			{ 8.0,  1.0},
		};
		{
			multi::array<double, 2> c({size(a({0, 1})), size(~b)});
			blas::gemm(1.0, a({0, 1}), b, 0.0, c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			multi::array<double, 2> c({size(a({0, 1})), size(~b)});
			blas::gemm_n(1.0, begin(a({0, 1})), size(a({0, 1})), begin(b), 0.0, begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm(1.0, ~(ar(extension(ar), {0, 1})), b, 0.0, ~c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm_n(1.0, begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0., begin(~c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[1][0] == 184.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_part_3x2) {
		using complex                    = std::complex<double>;
		multi::array<complex, 2> const a = {
			{1.0, 9.0, 1.0},
			{3.0, 3.0, 3.0},
		};
		BOOST_TEST( stride(~a) == 1 );
		BOOST_TEST( stride( a) == 3 );
		multi::array<complex, 2> const b = {
			{11.0, 12.0},
			{ 7.0, 19.0},
			{ 8.0,  1.0}
		};
		{
			multi::array<complex, 2> c({size(a({0, 1})), size(~b)});
			blas::gemm(1.0, a({0, 1}), b, 0.0, c);
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			multi::array<complex, 2> c({size(a({0, 1})), size(~b)});
			blas::gemm_n(1.0, begin(a({0, 1})), size(a({0, 1})), begin(b), 0.0, begin(c));
			BOOST_TEST( c[0][1] == 184.0 );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm(1.0, ~(ar(extension(ar), {0, 1})), b, 0.0, ~c);
			BOOST_TEST( c[1][0] == 184.0 );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm_n(1.0, begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0.0, begin(~c));
			BOOST_TEST( c[1][0] == 184.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1) {
		multi::array<double, 2> const a = {
			{1.0, 9.0, 1.0},
			{3.0, 3.0, 3.0},
		};
		BOOST_TEST( stride(~a) == 1 );
		BOOST_TEST( stride( a) == 3 );
		multi::array<double, 2> const b = {
			{11.0},
			{7.0},
			{8.0}
		};
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 82.0 );
			BOOST_TEST( c[1][0] == 78.0 );
		}
		{
			multi::array<double, 2> c({size(a), size(~b)});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[0][0] == 82.0 );
			BOOST_TEST( c[1][0] == 78.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm(1.0, ~(ar(extension(ar), {0, 1})), b, 0.0, ~c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			auto                    ar = +~a;
			multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
			blas::gemm_n(1., begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0., begin(~c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST( c[0][0] == 82.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1_bis) {
		multi::array<double, 2> const a = {
			{1.0, 9.0, 1.0},
			{3.0, 4.0, 5.0},
		};
		multi::array<double, 2> const b = {
			{11.0},
			{7.0},
			{8.0}
		};

		{
			multi::array<double, 2> c({1, 2});
			blas::gemm(1.0, a, b, 0.0, ~c);  // c⸆=ab, c=b⸆a⸆
			BOOST_TEST( (~c)[0][0] ==  82.0 );
			BOOST_TEST( (~c)[1][0] == 101.0 );
		}
		{
			multi::array<double, 2> c({1, 2});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(~c));  // c⸆=ab, c=b⸆a⸆
			BOOST_TEST( (~c)[0][0] ==  82.0 );
			BOOST_TEST( (~c)[1][0] == 101.0 );
		}
		{
			multi::array<double, 2> c({2, 1});
			blas::gemm(1.0, a, b, 0.0, c);  // c⸆=ab, c=b⸆a⸆
			BOOST_TEST( (~c)[0][1] == 101.0 );
			BOOST_TEST(    c[1][0] == 101.0 );
		}
		{
			multi::array<double, 2> c({2, 1});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c⸆=ab, c=b⸆a⸆
			BOOST_TEST( (~c)[0][1] == 101.0 );
			BOOST_TEST(    c[1][0] == 101.0 );
		}
		{
			multi::array<double, 2> c({1, 2});
			auto                    ar = +~a;
			blas::gemm(1., ~ar, b, 0., ~c);  // c⸆=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 101.0 );
		}
		{
			multi::array<double, 2> c({1, 2});
			auto                    ar = +~a;
			blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(~c));  // c⸆=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 101.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x1) {
		multi::array<double, 2> const a = {
			{1.0, 9.0, 1.0},
		};
		multi::array<double, 2> const b = {
			{11.0},
			{7.0},
			{8.0}
		};
		{
			multi::array<double, 2> c({1, 1});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			auto                    ar = +~a;
			blas::gemm(1.0, ~ar, b, 0.0, c);
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			auto                    ar = +~a;
			blas::gemm_n(1.0, begin(~ar), size(~ar), begin(b), 0.0, begin(c));
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			auto                    br = +~b;
			blas::gemm(1.0, a, ~br, 0.0, c);
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			BOOST_TEST( begin(c). stride() == 1 );
			BOOST_TEST( begin(c)->stride() == 1 );

			auto br = +~b;
			//  BOOST_TEST( begin(br). stride() == 1 );
			BOOST_TEST( begin( br)->stride() == 1 );

			BOOST_TEST(begin(a)->stride() == 1);
			BOOST_TEST( begin(~br). stride() == 1 );
			//  BOOST_TEST( begin(~br)->stride() == 1 );
			BOOST_TEST(begin(c)->stride() == 1);
			BOOST_TEST(begin(c).stride() == 1);
			BOOST_TEST(size(a) == 1);

			blas::gemm_n(1.0, begin(a), size(a), begin(~br), 0.0, begin(c));
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			auto                    br = +~b;
			blas::gemm(1.0, a, blas::H(br), 0.0, c);
			BOOST_TEST( c[0][0] == 82.0 );
		}
		{
			multi::array<double, 2> c({1, 1});
			auto                    br = +~b;
			blas::gemm_n(1.0, begin(a), size(a), begin(blas::H(br)), 0.0, begin(c));
			BOOST_TEST( c[0][0] == 82.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square) {
		using complex = std::complex<double>;
		constexpr complex              I{0.0, 1.0};
		multi::array<complex, 2> const a = {
			{1.0 + 3.0 * I, 3.0 + 2.0 * I},
			{9.0 + 1.0 * I, 7.0 + 1.0 * I},
		};
		multi::array<complex, 2> const b = {
			{11.0 + 2.0 * I, 12.0 + 4.0 * I},
			{ 7.0 + 1.0 * I, 19.0 - 9.0 * I},
		};
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1.0, a, b, 0.0, c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 145.0 + 43.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 145. + 43.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., ~a, b, 0., c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST(( c[1][1] == 170.-8.*I && c[1][0] == 77.+42.*I ));
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST(( c[1][1] == 170.-8.*I && c[1][0] == 77.+42.*I ));
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, ~b, 0., c);  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == 177.+69.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(~b), 0., begin(c));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == 177.+69.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::T(a), blas::T(b), 0., c);  // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == 109. + 68.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c));  // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == 109. + 68.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::T(a), blas::T(b), 0., blas::T(c));  // c⸆=a⸆b⸆, c=ba
			BOOST_TEST( c[0][1] == 109.+68.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(blas::T(c)));  // c⸆=a⸆b⸆, c=ba
			BOOST_TEST( c[0][1] == 109.+68.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x1) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I, 9. - 1. * I, 1. + 1. * I},
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I},
			{7. - 3. * I},
			{8. - 1. * I}
		};
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     ar = +~a;
			blas::gemm(1., ~ar, b, 0., c);  // c=ab, c⸆=ba
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     ar = +~a;
			blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(c));  // c=ab, c⸆=ba
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     br = +~b;
			blas::gemm(1., a, ~br, 0., c);
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     br = +~b;
			blas::context            ctxt;
			blas::gemm_n(ctxt, 1., begin(a), size(a), begin(~br), 0., begin(c));
			BOOST_TEST( c[0][0] == 84.-7.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     br = +~b;
			blas::gemm(1., a, blas::H(br), 0., ~c);
			BOOST_TEST( c[0][0] == 80. + 53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     br = +~b;
			blas::gemm_n(1., begin(a), size(a), begin(blas::H(br)), 0., begin(~c));
			BOOST_TEST( c[0][0] == 80. + 53.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_hermitized_square) {
		using complex = std::complex<double>;
		constexpr complex              I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 3. * I, 3. + 2. * I},
			{9. + 1. * I, 7. + 1. * I},
		};
		multi::array<complex, 2> const b = {
			{11. + 2. * I, 12. + 4. * I},
			{ 7. + 1. * I, 19. - 9. * I},
		};
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, b, 0., c);  // c=ab, c†=b†a†
			BOOST_TEST( c[1][0] == 145. + 43.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c†=b†a†
			BOOST_TEST( c[1][0] == 145. + 43.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), blas::H(b), 0., c);  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == 109.0 - 68.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c));  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == 109.0 - 68.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), blas::H(b), 0., blas::H(c));  // c†=a†b†, c=ba
			BOOST_TEST( c[1][0] == 184.0 - 40.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=a†b, c†=b†a
			BOOST_TEST( c[1][0] == 87.0 - 16.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=a†b, c†=b†a
			BOOST_TEST( c[1][0] == 87.0 - 16.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, blas::H(b), 0., c);  // c=ab†, c†=ba†
			BOOST_TEST( c[1][0] == 189.0 - 23.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			c = blas::gemm(1., a, blas::H(b));  // c=ab†, c†=ba†
			BOOST_TEST( c[1][0] == 189.0 - 23.0*I );
		}
		{
			multi::array<complex, 2> c = blas::gemm(1., a, blas::H(b));  // c=ab†, c†=ba†
			BOOST_TEST( size(c) == 2 );
			BOOST_TEST( c[1][0] == 189.0 - 23.0*I );
		}
		{
			auto c = multi::array<complex, 2>(blas::gemm(1., a, blas::H(b)));  // c=ab†, c†=ba†
			BOOST_TEST( size(c) == 2 );
			BOOST_TEST( c[1][0] == 189.0 - 23.0*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));  // c=ab†, c†=ba†
			BOOST_TEST( c[1][0] == 189. - 23.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), blas::H(b), 0., c);  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == 109. - 68.*I);
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c));  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == 109. - 68.*I);
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I},
			{9. - 1. * I},
			{1. + 1. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I},
			{7. - 3. * I},
			{8. - 1. * I}
		};
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 80.-53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 80.-53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=a†b, c†=b†a
			BOOST_TEST( c[0][0] == 80.-53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=a†b, c†=b†a
			BOOST_TEST( c[0][0] == 80.-53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     ha = +blas::hermitized(a);
			blas::gemm(1., ha, b, 0., c);
			BOOST_TEST( c[0][0] == 80.-53.*I );

			blas::gemm(1., blas::H(b), a, 0., c);
			BOOST_TEST( c[0][0] == 80.+53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			auto                     ha = +blas::hermitized(a);
			blas::gemm_n(1., begin(ha), size(ha), begin(b), 0., begin(c));
			BOOST_TEST( c[0][0] == 80.-53.*I );

			blas::gemm_n(1., begin(blas::H(b)), size(blas::H(b)), begin(a), 0., begin(c));
			BOOST_TEST( c[0][0] == 80.+53.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x2) {
		using complex = std::complex<double>;
		constexpr complex              I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I, 9. - 1. * I, 1. + 1. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I, 5. + 2. * I},
			{ 7. - 3. * I, 2. + 1. * I},
			{ 8. - 1. * I, 1. + 1. * I}
		};
		{
			multi::array<complex, 2> c({1, 2});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 20.+21.*I );
		}
		{
			multi::array<complex, 2> c({1, 2});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 20.+21.*I );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({1, 2});
			blas::gemm(1., blas::H(ar), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 28.+3.*I );
		}
		{
			auto                     ar = +~a;
			multi::array<complex, 2> c({1, 2});
			blas::gemm_n(1., begin(blas::H(ar)), size(blas::H(ar)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 28.+3.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x2) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I},
			{9. - 1. * I},
			{1. + 1. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I, 5. + 2. * I},
			{ 7. - 3. * I, 2. + 1. * I},
			{ 8. - 1. * I, 1. + 1. * I}
		};
		{
			multi::array<complex, 2> c({1, 2});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 28.+3.*I );
		}
		{
			multi::array<complex, 2> c({1, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][1] == 28.+3.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x2) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I, 5. + 2. * I},
			{9. - 1. * I, 9. + 1. * I},
			{1. + 1. * I, 2. + 2. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I, 5. + 2. * I},
			{ 7. - 3. * I, 2. + 1. * I},
			{ 8. - 1. * I, 1. + 1. * I}
		};
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 125.-84.*I );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 125.-84.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x1) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I, 5. + 2. * I},
			{9. - 1. * I, 9. + 1. * I},
			{1. + 1. * I, 2. + 2. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I},
			{7. - 3. * I},
			{8. - 1. * I}
		};
		{
			multi::array<complex, 2> c({2, 1});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 125.-84.*I );
		}
		{
			multi::array<complex, 2> c({2, 1});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 125.-84.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bis) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I},
			{9. - 1. * I},
			{1. + 1. * I}
		};
		multi::array<complex, 2> const b = {
			{11. - 2. * I},
			{7. - 3. * I},
			{8. - 1. * I}
		};
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 80. - 53.*I );
		}
		{
			multi::array<complex, 2> c({1, 1});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[0][0] == 80. - 53.*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square_automatic) {
		multi::array<double, 2> const a = {
			{1., 3.},
			{9., 7.},
		};
		multi::array<double, 2> const b = {
			{11., 12.},
			{ 7., 19.},
		};
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 148 && c[1][1] == 241 );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == 148 && c[1][1] == 241 );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1., a, blas::T(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][1] == 196. );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1., blas::T(a), b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][1] == 169. );
			BOOST_TEST( c[1][0] ==  82. );
		}
		{
			multi::array<double, 2> c({2, 2});
			blas::gemm(1., blas::T(a), blas::T(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][1] == 154. );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a = {
			{1. + 2. * I, 3. - 3. * I},
			{9. + 1. * I, 7. + 4. * I},
		};
		multi::array<complex, 2> const b = {
			{11. + 1. * I, 12. + 1. * I},
			{ 7. + 8. * I, 19. - 2. * I},
		};
		namespace blas = multi::blas;
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == complex(115, 104) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == complex(115, 104) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, blas::T(b), 0., c);  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == complex(178, 75) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(a), size(a), begin(blas::T(b)), 0., begin(c));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][0] == complex(178.0, 75.0) );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic_part2) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};

		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I, 3.0 - 3.0 * I},
			{9.0 + 1.0 * I, 7.0 + 4.0 * I},
		};
		multi::array<complex, 2> const b = {
			{11.0 + 1.0 * I, 12.0 + 1.0 * I},
			{ 7.0 + 8.0 * I, 19.0 - 2.0 * I},
		};
		namespace blas = multi::blas;
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::T(a), b, 0., c);  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST(( c[1][1] == complex(180, 29) && c[1][0] == complex(53, 54) ));
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(b), 0., begin(c));  // c=a⸆b, c⸆=b⸆a
			BOOST_TEST(( c[1][1] == complex(180, 29) && c[1][0] == complex(53, 54) ));
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::T(a), blas::T(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][1] == complex(186.0, 65.0) );
			BOOST_TEST( c[1][0] == complex(116, 25) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST(( c[1][1] == complex(186, 65) && c[1][0] == complex(116.0, 25.0) ));
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == complex(115, 104) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1.0, begin(a), size(a), begin(b), 0.0, begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][0] == complex(115, 104) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), b, 0., c);  // c=a†b, c†=b†a
			BOOST_TEST( c[1][0] == complex(111, 64) && c[1][1] == complex(158.0, -51.0) );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic_part3) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};

		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I, 3.0 - 3.0 * I},
			{9.0 + 1.0 * I, 7.0 + 4.0 * I},
		};
		multi::array<complex, 2> const b = {
			{11.0 + 1.0 * I, 12.0 + 1.0 * I},
			{ 7.0 + 8.0 * I, 19.0 - 2.0 * I},
		};
		namespace blas = multi::blas;
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c));  // c=a†b, c†=b†a
			BOOST_TEST( c[1][0] == complex(111, 64) && c[1][1] == complex(158, -51) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., a, blas::H(b), 0.0, c);  // c=ab†, c†=ba†
			BOOST_TEST( c[1][0] == complex(188, 43) && c[1][1] == complex(196, 25) );
			auto c2 = +blas::gemm(1., a, blas::H(b));
			BOOST_TEST( c2 == c );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(0.1, a, blas::H(b), 0.0, c);  // c=ab†, c†=ba†
			auto const c2 = +blas::gemm(0.1, a, blas::H(b));
			BOOST_TEST( c2 == c );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::H(a), blas::H(b), 0.0, c);  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == complex(116, -25) && c[1][1] == complex(186, -65) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c));  // c=a†b†, c†=ba
			BOOST_TEST( c[1][0] == complex(116, -25) && c[1][1] == complex(186, -65) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1., blas::T(a), blas::H(b), 0.0, c);  // c=a⸆b†, c†=ba⸆†
			BOOST_TEST( c[1][0] == complex(118, 5) && c[1][1] == complex(122, 45) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::H(b)), 0., begin(c));  // c=a⸆b†, c†=ba⸆†
			BOOST_TEST( c[1][0] == complex(118, 5) && c[1][1] == complex(122, 45) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm(1.0, blas::T(a), blas::T(b), 0.0, c);  // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == complex(116, 25) && c[1][1] == complex(186, 65) );
		}
		{
			multi::array<complex, 2> c({2, 2});
			blas::gemm_n(1.0, begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c));  // c=a⸆b⸆, c⸆=ba
			BOOST_TEST( c[1][0] == complex(116, 25) && c[1][1] == complex(186, 65) );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};

		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I, 3.0 - 3.0 * I, 1.0 - 9.0 * I},
			{9.0 + 1.0 * I, 7.0 + 4.0 * I, 1.0 - 8.0 * I},
		};
		multi::array<complex, 2> const b = {
			{11.0 + 1.0 * I, 12.0 + 1.0 * I, 4.0 + 1.0 * I, 8.0 - 2.0 * I},
			{ 7.0 + 8.0 * I, 19.0 - 2.0 * I, 2.0 + 1.0 * I, 7.0 + 1.0 * I},
			{ 5.0 + 1.0 * I,  3.0 - 1.0 * I, 3.0 + 8.0 * I, 1.0 + 1.0 * I}
		};
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(112, 12) );
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(112, 12) );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_realcomplex_complex_nonsquare_automatic) {
		using complex = std::complex<double>;

		complex const I{0.0, 1.0};

		multi::array<complex, 2> const a = {
			{1.0, 3.0, 1.0},
			{9.0, 7.0, 1.0},
		};
		multi::array<complex, 2> const b = {
			{11. + 1. * I, 12. + 1. * I, 4. + 1. * I, 8. - 2. * I},
			{ 7. + 8. * I, 19. - 2. * I, 2. + 1. * I, 7. + 1. * I},
			{ 5. + 1. * I,  3. - 1. * I, 3. + 8. * I, 1. + 1. * I}
		};
		{
			multi::array<complex, 2> c = blas::gemm(1., a, b);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(53, 24) );
		}
		{
			multi::array<complex, 2> c({2, 4});
			c = blas::gemm(1., a, b);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(53, 24) );
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm(1., a, b, 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(53, 24) );
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == complex(53, 24) );
		}
		{
			multi::array<double, 2> const a_real = {
				{1., 3., 1.},
				{9., 7., 1.},
			};
			multi::array<complex, 2> c({2, 4});
			blas::real_doubled(c) = blas::gemm(1., a_real, blas::real_doubled(b));

			BOOST_TEST( c[1][2] == complex(53, 24) );
		}
	}

	BOOST_AUTO_TEST_CASE(submatrix_result_issue_97) {
		using complex = std::complex<double>;

		constexpr complex I{0.0, 1.0};

		multi::array<complex, 2> M = {
			{2.0 + 3.0 * I, 2.0 + 1.0 * I, 1.0 + 2.0 * I},
			{4.0 + 2.0 * I, 2.0 + 4.0 * I, 3.0 + 1.0 * I},
			{7.0 + 1.0 * I, 1.0 + 5.0 * I, 0.0 + 3.0 * I},
		};
		auto M2 = +M({0, 3}, {0, 1});
		BOOST_TEST( M2 == M({0, 3}, {0, 1}) );
	}

	BOOST_AUTO_TEST_CASE(blas_context_gemm) {
		using complex = std::complex<double>;

		static constexpr complex I{0, 1};

		auto rand = [d = std::normal_distribution<>{}, g = std::mt19937{}]() mutable {  // NOLINT(cert-msc32-c, cert-msc51-cpp) for test purposes
			return d(g) + d(g) * I;
		};

		multi::array<complex, 2> A({30, 40});
		multi::array<complex, 2> B({40, 50});

		std::generate(A.elements().begin(), A.elements().end(), rand);
		std::generate(B.elements().begin(), B.elements().end(), rand);
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_gemm_range) {
		multi::array<double, 2> const a({2, 3}, 0.);
		multi::array<double, 2> const b({4, 3}, 0.);
		{
			multi::array<double, 2> c({2, 4});
			c() = blas::gemm(0.1, a, blas::H(b));
			BOOST_REQUIRE_CLOSE(c[1][2], 0., 0.00001);
		}
		{
			multi::array<double, 2> c = blas::gemm(0.1, a, blas::H(b));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][2] == 0. );
		}
		{
			multi::array<double, 2> const a = {
				{1, 3, 1},
				{9, 7, 1},
			};
			(void)a;
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_hermitized_second_gemm_range) {
		using complex = std::complex<double>;

		multi::array<complex, 2> const a({2, 3}, 0.);
		multi::array<complex, 2> const b({4, 3}, 0.);
		{
			multi::array<complex, 2> c({2, 4}, 999.);
			blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
			BOOST_TEST( c[1][2] != 999. );
		}
		{
			multi::array<complex, 2> c = blas::gemm(1., a, blas::H(b));  // c=ab⸆, c⸆=ba⸆
			BOOST_TEST( c[1][2] == 0. );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second) {
		namespace blas = multi::blas;

		multi::array<double, 2> const a = {
			{1, 3, 1},
			{9, 7, 1},
		};
		multi::array<double, 2> const b = {
			{11,  7, 5},
			{12, 19, 3},
			{ 4,  2, 3},
			{ 8,  7, 1}
		};
		{
			multi::array<double, 2> c({2, 4});
			blas::gemm(1., a, blas::H(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 53. );
		}
		{
			multi::array<double, 2> c({2, 4});
			blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 53. );
		}
		{
			multi::array<double, 2> c({2, 4});
			blas::gemm(0.1, a, blas::H(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
		{
			multi::array<double, 2> c({2, 4});
			blas::gemm_n(0.1, begin(a), size(a), begin(blas::H(b)), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
		{
			multi::array<double, 2> c({2, 4});
			c() = blas::gemm(0.1, a, blas::H(b));
		}
		{
			multi::array<double, 2> c = blas::gemm(0.1, a, blas::H(b));  // c=ab⸆, c⸆=ba⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 5.3, 0.00001);
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_real_nonsquare_hermitized_second) {
		namespace blas = multi::blas;
		using complex  = std::complex<double>;

		multi::array<complex, 2> const a = {
			{1., 3., 1.},
			{9., 7., 1.},
		};
		multi::array<complex, 2> const b = {
			{11.,  7., 5.},
			{12., 19., 3.},
			{ 4.,  2., 3.},
			{ 8.,  7., 1.}
		};
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm(1., a, blas::H(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_TEST( c[1][2] == 53. );
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(c[1][2], 53.0, 1E-6);
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm(0.1, a, blas::H(b), 0., c);  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(real(c[1][2]), 5.3, 0.00001);
		}
		{
			multi::array<complex, 2> c({2, 4});
			blas::gemm_n(0.1, begin(a), size(a), begin(blas::H(b)), 0., begin(c));  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE_CLOSE(real(c[1][2]), 5.3, 0.00001);
		}
		{
			multi::array<complex, 2> c({2, 4});
			c() = blas::gemm(0.1, a, blas::H(b));
		}
		{
			multi::array<complex, 2> c = blas::gemm(0.1, a, blas::H(b));  // c=ab⸆, c⸆=ba⸆
			BOOST_REQUIRE_CLOSE(real(c[1][2]), 5.3, 0.00001);
		}
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_1xn_complex) {
		using complex = std::complex<double>;
		multi::array<complex, 2> const a({1, 100}, 1.);
		multi::array<complex, 2> const b({1, 100}, 1.);

		multi::array<complex, 2> c({1, 1}, 999.);
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
		BOOST_TEST( c[0][0] == 100. );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_complex_inq_hydrogen_case) {
		using complex = std::complex<double>;
		complex const I{0, 1};

		multi::array<complex, 2> const a({3, 1}, 2. + 1. * I);
		multi::array<complex, 2> const b({1, 1}, 3. + 4. * I);

		multi::array<complex, 2> c({3, 1}, 999.);
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
		BOOST_TEST_REQUIRE( c[0][0] == (2. + 1.*I)*std::conj(3. + 4.*I) );
		BOOST_TEST_REQUIRE( c[1][0] == (2. + 1.*I)*std::conj(3. + 4.*I) );
		BOOST_TEST_REQUIRE( c[2][0] == (2. + 1.*I)*std::conj(3. + 4.*I) );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_complex_inq_hydrogen_case_no_n_interface) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a({3, 1}, 2. + 1. * I);
		multi::array<complex, 2> const b({10, 1}, 3. + 4. * I);

		multi::array<complex, 2> c({3, 10}, 999.);
		blas::gemm(1., a, blas::H(b), 0., c);
		BOOST_TEST_REQUIRE( c[0][0] == (2. + 1.*I)*std::conj(3. + 4.*I)  );
		BOOST_TEST_REQUIRE( c[1][0] == (2. + 1.*I)*std::conj(3. + 4.*I)  );
		BOOST_TEST_REQUIRE( c[0][1] == (2. + 1.*I)*std::conj(3. + 4.*I)  );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_1x1_complex_inq_hydrogen_case_complex_value_hermitized) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a({1, 1}, 2. + 1. * I);
		multi::array<complex, 2> const b({1, 1}, 3. + 4. * I);

		multi::array<complex, 2> c({1, 1}, 999.);
		c = blas::gemm(1., a, blas::H(b));
		BOOST_TEST( c[0][0] == (2. + 1.*I)*std::conj(3. + 4.*I) );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_1x1_complex_inq_hydrogen_case_complex_value) {
		using complex = std::complex<double>;
		complex const                  I{0, 1};
		multi::array<complex, 2> const a({1, 1}, 2. + 1. * I);
		multi::array<complex, 2> const b({1, 1}, 3. + 4. * I);

		multi::array<complex, 2> c({1, 1}, 999.);
		c = blas::gemm(1., a, b);
		BOOST_TEST( c[0][0] == (2. + 1.*I)*(3. + 4.*I) );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_1x1_complex_inq_hydrogen_case) {
		using complex = std::complex<double>;
		multi::array<complex, 2> const a({1, 1}, 2.);
		multi::array<complex, 2> const b({1, 1}, 3.);

		multi::array<complex, 2> c({1, 1}, 999.);
		c = blas::gemm(1., a, b);
		BOOST_TEST( c[0][0] == 6. );
	}

	BOOST_AUTO_TEST_CASE(blas_gemm_inq_case) {  // https://gitlab.com/correaa/boost-multi/-/issues/97
		using complex = std::complex<double>;
		complex const I{0, 1};

		multi::array<complex, 2> mat({10, 2}, 1.0 + 3. * I);
		multi::array<complex, 2> vec({10, 1}, -2.0 + 4. * I);

		mat({0, 10}, {1, 2}) = vec;

		namespace blas = multi::blas;

		{
			auto olap1 = +blas::gemm(1., blas::H(mat), vec);
			auto olap2 = +blas::gemm(1., blas::H(mat({0, 10}, {0, 1})), vec);

			BOOST_TEST( blas::H(mat)[1].size() == (~vec)[0].size() );
			BOOST_TEST( blas::dot(blas::H(mat)[0], (~vec)[0]) == olap1[0][0] );
			BOOST_TEST( std::inner_product(blas::H(mat)[0].begin(), blas::H(mat)[0].end(), (~vec)[0].begin(), complex{0}) == olap1[0][0] );

			multi::array<complex, 2> mat2  = mat({0, 10}, {0, 1});
			auto                     olap3 = +blas::gemm(1., blas::H(mat2), vec);

			BOOST_TEST(olap1[0][0] == olap2[0][0]);
			BOOST_TEST(olap3[0][0] == olap2[0][0]);
		}
		{
			multi::array<complex, 2> mat2  = mat({0, 3}, {0, 1});
			auto                     olap3 = +blas::gemm(1., blas::H(mat({0, 3}, {0, 1})), vec);
			BOOST_TEST( (+blas::gemm(1., blas::H(mat2), vec))[0][0] == (+blas::gemm(1., blas::H(mat({0, 3}, {0, 1})), vec))[0][0] );
		}
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109_part2) {
		multi::array<double, 2> const A({3, 4}, 5.);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const B({2, 3}, 7.);  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({2, 4}, 999.0);  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm(1.0, ~A, ~B, 0.0, ~C);

		BOOST_TEST_REQUIRE( C[0][0] == 105.0 );
		BOOST_TEST_REQUIRE( C[0][1] == 105.0 );
		BOOST_TEST_REQUIRE( C[1][0] == 105.0 );
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109) {
		multi::array<double, 2> const A({3, 4}, 5.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const B({2, 3}, 7.0);  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({4, 2}, 999.0);  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm(1., ~A, ~B, 0., C);

		BOOST_TEST_REQUIRE( C[0][0] == 105.0 );
		BOOST_TEST_REQUIRE( C[0][1] == 105.0 );
		BOOST_TEST_REQUIRE( C[1][0] == 105.0 );
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109_part2_complex) {
		multi::array<std::complex<double>, 2> const A({3, 4}, {5.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<std::complex<double>, 2> const B({2, 3}, {7.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<std::complex<double>, 2> C({2, 4}, {999.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm({1.0, 0.0}, ~A, ~B, {0.0, 0.0}, ~C);

		BOOST_TEST_REQUIRE( C[0][0] == 105.0 );
		BOOST_TEST_REQUIRE( C[0][1] == 105.0 );
		BOOST_TEST_REQUIRE( C[1][0] == 105.0 );
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109_complex) {
		multi::array<std::complex<double>, 2> const A({3, 4}, {5.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<std::complex<double>, 2> const B({2, 3}, {7.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<std::complex<double>, 2> C({4, 2}, {999.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm({1.0, 0.0}, ~A, ~B, {0.0, 0.0}, C);

		BOOST_TEST( C[0][0] == 105.0 );
		BOOST_TEST( C[0][1] == 105.0 );
		BOOST_TEST( C[1][0] == 105.0 );
	}
#endif

	BOOST_AUTO_TEST_CASE(blas_issue_109_complex_mx2) {
		multi::array<std::complex<double>, 2> const A({3, 4}, {5.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<std::complex<double>, 2> const B({2, 3}, {7.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<std::complex<double>, 2> C({4, 2}, {999.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm({1.0, 0.0}, ~A, ~B, {0.0, 0.0}, C);

		BOOST_TEST( C[0][0] == 105.0 );
		BOOST_TEST( C[1][0] == 105.0 );
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109_complex_mx1) {
		multi::array<std::complex<double>, 2> const A({3, 4}, {5.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<std::complex<double>, 2> const B({1, 3}, {7.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<std::complex<double>, 2> C({4, 1}, {999.0, 0.0});  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm(std::complex<double>{1.0, 0.0}, ~A, ~B, std::complex<double>{0.0, 0.0}, C);

		BOOST_TEST( C[0][0] == 105.0 );
		BOOST_TEST( C[1][0] == 105.0 );
	}

	BOOST_AUTO_TEST_CASE(blas_issue_109_double_mx1) {
		multi::array<double, 2> const A({3, 4}, 5.0);  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 2> const B({1, 3}, 7.0);  // NOLINT(readability-identifier-length) BLAS naming

		multi::array<double, 2> C({4, 1}, 999.0);  // NOLINT(readability-identifier-length) BLAS naming
		blas::gemm(1.0, ~A, ~B, 0.0, C);

		BOOST_TEST( C[0][0] == 105.0 );
		BOOST_TEST( C[1][0] == 105.0 );
	}

	return boost::report_errors();
}  // NOLINT(readability/fn_size)
