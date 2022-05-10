// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS axpy"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "config.hpp"

#include "../../../array.hpp"
#include       "../../blas/axpy.hpp"
#include       "../../blas/operations.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_axpy_real) {
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = A;
	multi::array<double, 1> const B = A[2];

	blas::axpy(2., B, A[1]); // daxpy
	BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_double) {
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	multi::array<double, 2> A = cA;
	multi::array<double, 1> const b = cA[2];

	blas::axpy(2., b, A[1]); // A[1] = 2*b + A[1], A[1]+= a*A[1]
	BOOST_REQUIRE( A[1][2] == 2.*b[2] + cA[1][2] );

	using complex = std::complex<double>; complex const I = {0, 1};
	multi::array<complex, 1> AC = {1. + 2.*I, 3. + 4.*I, 4. - 8.*I};
	multi::array<complex, 1> BC(extensions(AC), complex{0.});

	blas::axpy(+1., blas::real(AC), blas::real(BC));
	blas::axpy(-1., blas::imag(AC), blas::imag(BC));

	BOOST_REQUIRE( BC[2] == std::conj(AC[2]) );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex) {
	using complex = std::complex<double>;
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

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_plus_equal) {
	using complex = std::complex<double>;
	multi::array<complex, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = A;
	multi::array<complex, 1> const B = A[2];
	A[1] += blas::axpy(2., B); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_minus_equal) {
	using complex = std::complex<double>;
	multi::array<complex, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = A;
	multi::array<complex, 1> const B = A[2];
	A[1] -= blas::axpy(2., B); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( A[1][2] == -2.*B[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_context) {
	using complex = std::complex<double>;
	multi::array<complex, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	auto const AC = A;
	multi::array<complex, 1> const B = A[2];
	blas::axpy(blas::context{}, 2., B, A[1]); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
	BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
}

BOOST_AUTO_TEST_CASE(multi_blas_axpy_operator_minus) {
	using complex = std::complex<double>;
	multi::array<complex, 1> x = {10., 11., 12., 13.};
	multi::array<complex, 1> y = x;

	using blas::operators::operator-;
	using blas::operators::operator+;
	using blas::operators::operator-=;

	BOOST_REQUIRE( (x - y)[0] == 0. );
	BOOST_REQUIRE( (y - x)[0] == 0. );

	BOOST_REQUIRE( (x - (y+y))[0] == -x[0] );
	BOOST_REQUIRE( ((x+x) - y)[0] == +x[0] );

	multi::array<complex, 2> A = {{1., 2.}, {3., 4.}};
	multi::array<complex, 1> B = {1., 2.};
	BOOST_REQUIRE( (A[0] - B)[0] == 0. );
	BOOST_REQUIRE( (A[0] - B)[1] == 0. );

	multi::array<complex, 1> X = {10., 11., 12., 13.};
	multi::array<complex, 1> Y = {10., 11., 12., 13.};
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

