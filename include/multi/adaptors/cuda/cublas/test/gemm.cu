// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS gemm"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/cuda/cublas.hpp>

#include <multi/adaptors/blas/gemm.hpp>
#include <multi/adaptors/blas/axpy.hpp>
#include <multi/adaptors/blas/nrm2.hpp>
#include <multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 2> const M_gpu = {
		{ { 9.0, 0.0}, {24.0, 0.0}, {30.0, 0.0}, {9.0, 0.0} },
		{ { 4.0, 0.0}, {10.0, 0.0}, {12.0, 0.0}, {7.0, 0.0} },
		{ {14.0, 0.0}, {16.0, 0.0}, {36.0, 0.0}, {1.0, 0.0} },
	};

	multi::thrust::cuda::array<complex, 1> const X_gpu = { {1.1, 0.0}, {2.1, 0.0}, {3.1, 0.0}, {4.1, 0.0} };

	multi::thrust::cuda::array<complex, 1> Y_gpu = { {4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0} };

	blas::gemv(/*alpha*/ 1.1, M_gpu, X_gpu, /*beta*/ 1.2, Y_gpu);  // y = a*M*x + b*y

	multi::array<complex, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_REQUIRE_SMALL( +blas::nrm2(Y_copy - multi::array<complex, 1>{ {214.02, 0.0}, {106.43, 0.0}, {188.37, 0.0} }) , 1e-13);
}

BOOST_AUTO_TEST_CASE(cublas_gemv_real) {
	namespace blas = multi::blas;
	using T = double;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<T, 2> const M_gpu = {
		{  9.0, 24.0, 30.0, 9.0 },
		{  4.0, 10.0, 12.0, 7.0 },
		{ 14.0, 16.0, 36.0, 1.0 },
	};

	multi::thrust::cuda::array<T, 1> const X_gpu = { 1.1, 2.1, 3.1, 4.1 };

	multi::thrust::cuda::array<T, 1> Y_gpu = { 4.0, 5.0, 6.0 };

	blas::gemv(/*alpha*/ 1.1, M_gpu, X_gpu, /*beta*/ 1.2, Y_gpu);  // y = a*M*x + b*y

	multi::array<T, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_REQUIRE_SMALL( +blas::nrm2(Y_copy - multi::array{ 214.02, 106.43, 188.37 }) , 1e-13);
}

BOOST_AUTO_TEST_CASE(cublas_gemm_nh) {
	namespace blas  = multi::blas;

	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 2> const a = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	{
		multi::thrust::cuda::array<complex, 2> c({2, 2}, {9999.0, 0.0});  // NOLINT(readability-identifier-length) conventional BLAS naming
		blas::gemm({1.0, 0.0}, a, a, {0.0, 0.0}, c);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_REQUIRE( c_copy[1][0] == 16.0 -  2.0*I );
		BOOST_REQUIRE( c_copy[0][1] == 14.0 - 38.0*I );
	}
	{
		auto const c = +blas::gemm(complex{1.0, 0.0}, a, a);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_REQUIRE( c_copy[1][0] == 16.0 -  2.0*I );
		BOOST_REQUIRE( c_copy[0][1] == 14.0 - 38.0*I );
	}
	{
		multi::thrust::cuda::array<complex, 2> c({2, 2}, {0.0, 0.0});
		c += blas::gemm(complex{1.0, 0.0}, a, a);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_REQUIRE( c_copy[1][0] == 16.0 -  2.0*I );
		BOOST_REQUIRE( c_copy[0][1] == 14.0 - 38.0*I );
	}
}