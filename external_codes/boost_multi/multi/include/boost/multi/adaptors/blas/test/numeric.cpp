// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/numeric.hpp>
#include <boost/multi/adaptors/blas/operations.hpp>
// IWYU pragma: no_include "boost/multi/adaptors/blas/complex_traits.hpp"  // for blas  // needed by iwyu-clang-macos

#include <boost/multi/array.hpp>

#include <cmath>
#include <complex>

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag) {
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		namespace blas                       = multi::blas;
		multi::array<complex, 1> const array = {1.0 + 2.0 * I, 3.0 + 5.0 * I, 9.0 + 2.0 * I};
		BOOST_TEST( blas::imag(array)[2] == 2.0 );
		BOOST_TEST( blas::real(array)[2] == 9.0 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_conjugated) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> array = {
			{1.0 - 3.0 * I, 6.0 + 2.0 * I},
			{8.0 + 2.0 * I, 2.0 + 4.0 * I},
			{2.0 - 1.0 * I, 1.0 + 1.0 * I}
		};
		BOOST_TEST( array[0][0] == 1.0 - 3.0*I );

		multi::array<complex, 2> const carray = {
			{1.0 - 3.0 * I, 6.0 + 2.0 * I},
			{8.0 + 2.0 * I, 2.0 + 4.0 * I},
			{2.0 - 1.0 * I, 1.0 + 1.0 * I}
		};
		BOOST_TEST( carray[0][0] == 1.0 - 3.0*I );

		namespace blas = multi::blas;
		auto conjr     = blas::make_conjugater(array.data_elements());

		decltype(blas::make_conjugater(carray.data_elements())) ppp;  // = BdataC;
		ppp = conjr;

		BOOST_TEST( *ppp == 1.0 + 3.0*I );

		//  static_assert(    multi::blas::is_complex_array<multi::array<thrust::complex<double>, 2>>{}, "!");
		static_assert(blas::is_complex_array<decltype(array)>{});
		static_assert(!blas::is_conjugated<decltype(array)>{});

		auto&& conjd_array = blas::conj(array);
		static_assert(blas::is_conjugated<decltype(conjd_array)>{});

		BOOST_TEST( conjd_array[0][0] == 1.0 + 3.0*I );
		BOOST_TEST( std::abs( imag(*conjd_array.base()) - +3.0 ) < 1e-10 );

		//  BOOST_TEST_REQUIRE( base(Bconj)->imag() == +3 );
		BOOST_TEST( conjd_array[0][1] == conjd_array.rotated()[1][0] );
		BOOST_TEST( conjd_array.rotated()[1][0] == conjd_array[0][1] );

		//  BOOST_TEST( base(Bconj) == -3.0*I );
		static_assert(blas::is_complex_array<decltype(conjd_array)>{});

		BOOST_TEST( blas::conj(conjd_array) == array );

		BOOST_TEST( blas::conj(array)[1][0] == std::conj(array[1][0]) );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_numeric_decay) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> arr = {
			{1.0 - 3.0 * I, 6.0 + 2.0 * I, 9.0 + 3.0 * I},
			{8.0 + 2.0 * I, 2.0 + 4.0 * I, 9.0 + 3.0 * I},
			{2.0 - 1.0 * I, 1.0 + 1.0 * I, 9.0 + 3.0 * I},
			{9.0 + 3.0 * I, 9.0 + 3.0 * I, 9.0 + 3.0 * I}
		};

		namespace blas = multi::blas;
		multi::array<complex, 2>       conj_arr{blas::conj(arr)};
		multi::array<complex, 2> const conj_arr2{blas::conj(arr)};

		BOOST_TEST( conj_arr[2][1] == std::conj(arr[2][1]) );
		BOOST_TEST( blas::conj(arr)[2][1] == std::conj(arr[2][1]) );

		BOOST_TEST( blas::T(arr)[1][2] == arr[2][1] );
		BOOST_TEST( blas::T(arr) == ~arr );

		BOOST_TEST( blas::conj(arr)[1][2]    == blas::hermitized(arr)[2][1] );
		BOOST_TEST( blas::conj(blas::T(arr)) == blas::hermitized(arr) );

		BOOST_TEST( blas::hermitized(arr)[2][1] == blas::conj(arr)[1][2] );
		BOOST_TEST( blas::hermitized(arr)       == blas::conj(blas::T(arr)) );

		BOOST_TEST( std::abs( blas::real(arr)[2][1] - std::real(arr[2][1]) ) < 1e-10 );
		BOOST_TEST( std::abs( blas::imag(arr)[2][1] - std::imag(arr[2][1]) ) < 1e-10 );

		multi::array<double, 2> const B_real_doubled = {
			{1.0, -3.0, 6.0, 2.0, 9.0, 3.0},
			{8.0,  2.0, 2.0, 4.0, 9.0, 3.0},
			{2.0, -1.0, 1.0, 1.0, 9.0, 3.0},
			{9.0,  3.0, 9.0, 3.0, 9.0, 3.0}
		};
		BOOST_TEST( blas::real_doubled(arr).sizes() == B_real_doubled.sizes() );
		BOOST_TEST( blas::real_doubled(arr)         == B_real_doubled         );
	}

#if defined(CUDA_FOUND) and CUDA_FOUND
	#include <thrust/complex.h>

	BOOST_AUTO_TEST_CASE(multi_blas_numeric_decay_thrust) {
		using complex = thrust::complex<double>;
		complex const I{0.0, 1.0};

		multi::array<complex, 2> B = {
			{1.0 - 3.0 * I, 6.0 + 2.0 * I},
			{8.0 + 2.0 * I, 2.0 + 4.0 * I},
			{2.0 - 1.0 * I, 1.0 + 1.0 * I},
		};

		namespace blas                 = multi::blas;
		multi::array<complex, 2> conjB = blas::conj(B);
		BOOST_TEST( conjB[1][2] == conj(B[1][2]) );
	}
#endif

	BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_imag_part) {
		using complex = std::complex<double>;
		complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<double, 2> arr = {
			{1.0, 3.0, 4.0},
			{9.0, 7.0, 1.0}
		};
		multi::array<complex, 2> complex_arr = arr;
		BOOST_TEST( complex_arr[1][1] == arr[1][1] );

		multi::array<complex, 2> arr2 = {
			{1.0 - 3.0 * I, 6.0 + 2.0 * I},
			{8.0 + 2.0 * I, 2.0 + 4.0 * I},
			{2.0 - 1.0 * I, 1.0 + 1.0 * I}
		};

		multi::array<double, 2> const arr2_real = {
			{1.0, 6.0},
			{8.0, 2.0},
			{2.0, 1.0},
		};
		multi::array<double, 2> const arr2_imag = {
			{-3.0, +2.0},
			{+2.0, +4.0},
			{-1.0, +1.0},
		};

		using multi::blas::imag;
		using multi::blas::real;

		BOOST_TEST( arr2_real == real(arr2) );
		BOOST_TEST( real(arr2) == arr2_real );
		BOOST_TEST( imag(arr2) == arr2_imag );

		BOOST_TEST( arr2[1][0] == 8.0 + 2.0*I );
		BOOST_TEST( arr2[1][0].imag() == 2.0 );

		namespace blas = multi::blas;
		BOOST_TEST( blas::hermitized(arr2)[1][2] == std::conj( arr2[2][1] ) );

		blas::hermitized(arr2)[1][2] = 20.0 + 30.0 * I;
		BOOST_TEST( arr2[2][1] == 20.0 - 30.0*I );
	}
	return boost::report_errors();
}
