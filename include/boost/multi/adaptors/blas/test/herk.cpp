// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/blas/filling.hpp>     // for filling
#include <boost/multi/adaptors/blas/gemm.hpp>        // for gemm, gemm_range
#include <boost/multi/adaptors/blas/herk.hpp>        // for herk
#include <boost/multi/adaptors/blas/nrm2.hpp>        // for nrm2_ref, nrm2
#include <boost/multi/adaptors/blas/numeric.hpp>     // for involuted, under...
#include <boost/multi/adaptors/blas/operations.hpp>  // for H, T, (anonymous)
// IWYU pragma: no_include "boost/multi/adaptors/blas/traits.hpp"      // for blas  // needed for iwyu-clang-macos

#include <boost/multi/array.hpp>                     // for array, layout_t

#include <cmath>  // for sqrt
// IWYU pragma: no_include <cstdlib>
#include <complex>      // for operator*, opera...
#include <iostream>     // for operator<<, basi...
#include <iterator>     // for size
#include <limits>       // for numeric_limits
#include <string>       // for char_traits, bas...
#include <type_traits>  // for is_same

namespace multi = boost::multi;

// NOLINTNEXTLINE(fuchsia-default-arguments-declarations,fuchsia-default-arguments-calls)
template<class M> auto print(M const& mat, std::string const& msg = "") -> decltype(auto) {
	using multi::size;
	using std::cout;
	cout << msg << "\n"
		 << '{';
	for(int i = 0; i != size(mat); ++i) {
		cout << '{';
		for(auto j : mat[i].extension()) {  // NOLINT(altera-unroll-loops)
			cout << mat[i][j];
			if(j + 1 != size(mat[i])) {
				cout << ", ";
			}
		}
		cout << '}' << '\n';
		if(i + 1 != size(mat)) {
			cout << ", ";
		}
	}
	return cout << '}' << '\n';
}

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_blas_herk) {
		namespace blas = multi::blas;
		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<complex, 2> const a = {
			{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
			{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I},
		};
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) conventional name in BLAS
			blas::herk(a, c);
			BOOST_TEST( c[1][0] == complex(50.0, -49.0) );
			BOOST_TEST( c[0][1] == complex(50.0, +49.0) );

			multi::array<complex, 2> const c_copy = blas::herk(1.0, a);
			BOOST_TEST( c == c_copy );

			BOOST_TEST( +blas::gemm(1.0, a, blas::H(a)) == blas::herk(a) );
		}
	}

	BOOST_AUTO_TEST_CASE(inq_case) {
		namespace blas = multi::blas;
		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const a = {
			{0.0,  1.0,  2.0},
			{3.0,  4.0,  5.0},
			{6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0},
		};

		BOOST_TEST( (+blas::gemm(1.0, a, blas::T(a)))[1][2] == 86.0 );
		{
			multi::array<double, 2> c({4, 4});  // NOLINT(readability-identifier-length) conventional name in BLAS
			blas::herk(1.0, a, c);
			BOOST_TEST( c[1][2] == (+blas::gemm(1.0, a, blas::T(a)))[1][2] );
			//  BOOST_TEST( c[2][1] == (+blas::gemm(1., a, blas::T(a)))[2][1] );
		}
		{
			multi::array<double, 2> const c = blas::herk(1.0, a);  // NOLINT(readability-identifier-length) conventional name in BLAS
			BOOST_TEST( c == +blas::gemm(1., a, blas::T(a)) );
			BOOST_TEST( blas::herk(a) == +blas::gemm(1.0, a, blas::T(a)) );
			BOOST_TEST( blas::herk(2.0, a) == +blas::gemm(2.0, a, blas::T(a)) );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk_real) {
		namespace blas = multi::blas;
		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const a = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0},
		};
		{
			multi::array<double, 2> c({2, 2}, 9999.0);  // NOLINT(readability-identifier-length) BLAS naming
			blas::herk(1.0, a, c);
			BOOST_TEST( c[0][1] == 34.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case) {
		namespace blas = multi::blas;
		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const a = {
			{1.0, 2.0, 3.0},
		};
		multi::array<double, 2> b = blas::herk(a);  // NOLINT(readability-identifier-length) BLAS naming

		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( b[0][0] == 1.0*1.0 + 2.0*2.0 + 3.0*3.0 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case_scale) {
		namespace blas = multi::blas;
		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<double, 2> const a = {
			{1.0, 2.0, 3.0},
		};

		multi::array<double, 2> b = blas::herk(0.1, a);  // NOLINT(readability-identifier-length) BLAS naming

		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( std::abs( b[0][0] - (1.0*1.0 + 2.0*2.0 + 3.0*3.0)*0.1 ) < 1E-6 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case) {
		namespace blas = multi::blas;

		using complex = std::complex<double>;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const a = {
			{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}},
		};
		multi::array<complex, 2> b = blas::herk(1.0, a);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( b[0][0] == 1.0*1.0 + 2.0*2.0 + 3.0*3.0 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case_scale) {
		namespace blas = multi::blas;

		using complex = std::complex<double>;
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 2> const a = {
			{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}},
		};
		multi::array<complex, 2> b = blas::herk(0.1, a);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( std::abs( real( b[0][0]/0.1 ) - (1.0*1.0 + 2.0*2.0 + 3.0*3.0) ) < 1E-6 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case) {
		namespace blas = multi::blas;

		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I, 2.0 + 3.0 * I, 3.0 + 4.0 * I},
		};
		multi::array<complex, 2> b = blas::herk(a);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( b[0][0] == std::norm(1.0 + 2.0*I) + std::norm(2.0 + 3.0*I) + std::norm(3.0 + 4.0*I) );

		BOOST_TEST( std::sqrt(real(blas::herk(a)[0][0])) == blas::nrm2(a[0]) );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_out_param) {
		namespace blas = multi::blas;
		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const a = {{1.0 + 2.0 * I}, {2.0 + 3.0 * I}, {3.0 + 4.0 * I}};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<complex, 2>       b({1, 1});                                                // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST( size(b) == 1 );

		blas::herk(blas::filling::upper, 1.0, blas::H(a), 0.0, b);

		BOOST_TEST( b[0][0] == std::norm(1.0 + 2.0*I) + std::norm(2.0 + 3.0*I) + std::norm(3.0 + 4.0*I) );

		//  BOOST_TEST( std::sqrt(real(b[0][0])) == blas::nrm2(blas::T(a)[0])() );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized) {
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) conventional name in BLAS

		// NOLINTNEXTLINE(readability-identifier-length) conventional name in BLAS
		multi::array<complex, 2> const a = {
			{1.0 + 2.0 * I},
			{2.0 + 3.0 * I},
			{3.0 + 4.0 * I},
		};

		namespace blas = multi::blas;

		multi::array<complex, 2> b = blas::herk(blas::H(a));  // NOLINT(readability-identifier-length) BLAS naming

		BOOST_TEST( size(b) == 1 );
		BOOST_TEST( b[0][0] == std::norm(1.0 + 2.0*I) + std::norm(2.0 + 3.0*I) + std::norm(3.0 + 4.0*I) );

		BOOST_TEST( std::sqrt(real(blas::herk(blas::H(a))[0][0])) == blas::nrm2(a.rotated()[0]) );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_auto) {
		namespace blas = multi::blas;
		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const arr = {
			{1.0 + 2.0 * I},
			{2.0 + 3.0 * I},
			{3.0 + 4.0 * I},
		};
		auto arr2 = blas::herk(1.0, blas::hermitized(arr));
		static_assert(std::is_same<decltype(arr2), multi::array<complex, 2>>{});
		BOOST_TEST( size(arr2) == 1 );
		BOOST_TEST( arr2[0][0] == std::norm(1.0 + 2.0*I) + std::norm(2.0 + 3.0*I) + std::norm(3.0 + 4.0*I) );

		BOOST_TEST( std::sqrt(real(blas::herk(blas::H(arr))[0][0])) == blas::nrm2(arr.rotated()[0]) );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_identity) {
		namespace blas = multi::blas;
		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const arr = {
			{1.0 + 3.0 * I, 3.0 - 2.0 * I, 4.0 + 1.0 * I},
			{9.0 + 1.0 * I, 7.0 - 8.0 * I, 1.0 - 3.0 * I},
		};

		{
			multi::array<complex, 2> arr2({2, 2}, {9999.0, 0.0});   // NOLINT(readability-identifier-length) conventional one-letter operation BLAS
			blas::herk(blas::filling::lower, 1.0, arr, 0.0, arr2);  // c^dagger = c = a a^dagger = (a a^dagger)^dagger, `c` is lower triangular
			BOOST_TEST(( arr2[1][0] == complex{50.0, -49.0} ));
			BOOST_TEST( arr2[0][1] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) conventional one-letter operation BLAS
			static_assert(blas::is_conjugated<decltype(blas::H(c))>::value);

			blas::herk(blas::filling::lower, 1.0, arr, 0.0, blas::H(c));  // c^dagger = c = a a^dagger = (aa^dagger)^daggerr, `c` in upper triangular

			BOOST_TEST(( blas::H(c)[1][0] == complex{50.0, -49.0} ));
			BOOST_TEST( blas::H(c)[0][1] == 9999.0 );
		}
		{
			// NOLINTNEXTLINE(readability-identifier-length) : conventional one-letter operation BLASs
			multi::array<complex, 2> c({3, 3}, {9999.0, 0.0});
			herk(blas::filling::lower, 1.0, blas::T(arr), 0.0, blas::T(c));  // c†=c=aT(aT)† not supported
			BOOST_TEST(( c.transposed()[1][0] == complex{52.0, -90.0} ));
			BOOST_TEST(  c.transposed()[0][1] == 9999.0                );
		}
		{
			// NOLINTNEXTLINE(readability-identifier-length) : conventional one-letter operation BLASs
			multi::array<complex, 2> c({3, 3}, {9999.0, 0.0});
			blas::herk(blas::filling::lower, 1.0, blas::T(arr), 0.0, blas::H(blas::T(c)));  // c†=c=aT(aT)† not supported
			BOOST_TEST(( blas::H(blas::T(c))[1][0] == complex{52.0, -90.0} ));
			BOOST_TEST( blas::H(blas::T(c))[0][1] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});   // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
			blas::herk(blas::filling::upper, 1.0, arr, 0.0, c);  // c†=c=aa†=(aa†)†, `c` in upper triangular
			BOOST_TEST(( c[0][1] == complex{50.0, +49.0} ));
			BOOST_TEST( c[1][0] == 9999.0 );
		}
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
			blas::herk_nm(blas::filling::upper, 1.0, arr.home(), arr.size(), (~arr).size(), 0.0, c.home());
			BOOST_TEST(( c[0][1] == complex{50.0, +49.0} ));
			BOOST_TEST( c[1][0] == 9999.0 );
		}
		// {
		//  multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
		//  c() = blas::herk(blas::filling::upper, 1.0, arr);
		//  BOOST_TEST(( c[0][1] == complex{50.0, +49.0} ));
		//  // BOOST_TEST( c[1][0] == 9999.0 );
		// }
		{
			multi::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
			blas::herk(1.0, arr, c);                             // c†=c=aa†=(aa†)†
			BOOST_TEST(( c[0][1] == complex{50.0, +49.0} ));
			BOOST_TEST(( c[1][0] == complex{50.0, -49.0} ));
		}
		{
			multi::array<complex, 2> c({3, 3}, {9999.0, 0.0});            // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
			blas::herk(blas::filling::lower, 1.0, blas::H(arr), 0.0, c);  // c†=c=aa†=(aa†)†, `c` in lower triangular
			BOOST_TEST(( c[1][0] == complex{52.0, 90.0} ));
			BOOST_TEST( c[0][1] == 9999.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_square) {
		namespace blas = multi::blas;

		using complex  = std::complex<double>;
		auto const I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
		auto const nan = std::numeric_limits<double>::quiet_NaN();

		// NOLINTNEXTLINE(readability-identifier-length) lapack conventional name
		multi::array<complex, 2> const A = {
			{12.9388 + I * 0.0, 9.80028 + I * -0.00011091, 9.66966 + I * -0.0114817},
			{    nan + I * nan,         8.44604 + I * 0.0,  3.78646 + I * 0.0170734},
			{    nan + I * nan,             nan + I * nan,        7.70655 + I * 0.0},
		};

		// NOLINTNEXTLINE(readability-identifier-length) lapack conventional name
		multi::array<complex, 2> C({3, 3}, complex{0.0, 0.0});

		blas::herk(boost::multi::blas::filling::upper, complex{1.0, 0.0}, A, complex{0.0, 0.0}, C);
	}

	return boost::report_errors();
}
