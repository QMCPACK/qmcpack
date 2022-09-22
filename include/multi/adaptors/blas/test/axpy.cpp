// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2022

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS axpy"
#define BOOST_TEST_DYN_LINK 
#include<boost/test/unit_test.hpp>

#include "config.hpp"

#include "multi/adaptors/blas/axpy.hpp"
#include "multi/adaptors/blas/operations.hpp"
#include "multi/array.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_axpy_real) {
	multi::array<double, 2> a = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = a;
	multi::array<double, 1> const b = a[2];  // NOLINT(readability-identifier-length) BLAS naming

	blas::axpy(2., b, a[1]); // daxpy
	BOOST_REQUIRE( a[1][2] == 2.*b[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_double) {
	multi::array<double, 2> const const_arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	multi::array<double, 2> arr = const_arr;
	multi::array<double, 1> const b = const_arr[2];  // NOLINT(readability-identifier-length) conventional name in BLAS

	blas::axpy(2., b, arr[1]); // A[1] = 2*b + A[1], A[1]+= a*A[1]
	BOOST_REQUIRE( arr[1][2] == 2.*b[2] + const_arr[1][2] );

	using complex = std::complex<double>; complex const I = {0, 1};  // NOLINT(readability-identifier-length) imaginary unit
	multi::array<complex, 1> AC = {1. + 2.*I, 3. + 4.*I, 4. - 8.*I};
	multi::array<complex, 1> BC(extensions(AC), complex{0.});

	blas::axpy(+1., blas::real(AC), blas::real(BC));
	blas::axpy(-1., blas::imag(AC), blas::imag(BC));

	BOOST_REQUIRE( BC[2] == std::conj(AC[2]) );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex) {
	using complex = std::complex<double>;
	multi::array<complex, 2> arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const const_arr = arr;
	multi::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
	blas::axpy(2., x, arr[1]);  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( arr[1][2] == 2.*x[2] + const_arr[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_plus_equal) {
	using complex = std::complex<double>;
	multi::array<complex, 2> arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const carr = arr;
	multi::array<complex, 1> const y = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
	arr[1] += blas::axpy(2., y);  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( arr[1][2] == 2.*y[2] + carr[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_minus_equal) {
	using complex = std::complex<double>;
	multi::array<complex, 2> arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = arr;
	multi::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
	arr[1] -= blas::axpy(2., x);  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( arr[1][2] == -2.*x[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_context) {
	using complex = std::complex<double>;
	multi::array<complex, 2> arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const arr_copy = arr;
	multi::array<complex, 1> const arr2 = arr[2];
	blas::axpy(blas::context{}, 2., arr2, arr[1]); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( arr[1][2] == 2.*arr2[2] + arr_copy[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_operator_minus) {
	using complex = std::complex<double>;
	multi::array<complex, 1> x = {10., 11., 12., 13.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> y = x;  // NOLINT(readability-identifier-length) BLAS naming

	using blas::operators::operator-;

	BOOST_REQUIRE( (x - y)[0] == 0. );
	BOOST_REQUIRE( (y - x)[0] == 0. );

	using blas::operators::operator+;
	
	BOOST_REQUIRE( (x - (y+y))[0] == -x[0] );
	BOOST_REQUIRE( ((x+x) - y)[0] == +x[0] );

	multi::array<complex, 2> arr = {{1., 2.}, {3., 4.}};
	multi::array<complex, 1> arr2 = {1., 2.};
	BOOST_REQUIRE( (arr[0] - arr2)[0] == 0. );
	BOOST_REQUIRE( (arr[0] - arr2)[1] == 0. );

	multi::array<complex, 1> X = {10., 11., 12., 13.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> Y = {10., 11., 12., 13.};  // NOLINT(readability-identifier-length) BLAS naming

	using blas::operators::operator-=;
	X -= Y;
	BOOST_REQUIRE( X[0] == 0. );
}

#if CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_thrust) {
	using complex = thrust::complex<double>;
	multi::array<complex, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = A;
	multi::array<complex, 1> const B = A[2];
	blas::axpy(2., B, A[1]); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
}
#endif

