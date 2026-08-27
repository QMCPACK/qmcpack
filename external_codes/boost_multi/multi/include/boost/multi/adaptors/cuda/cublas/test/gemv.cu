// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023-2024 Alfredo A. Correa

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/gemv.hpp>
#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

#define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )

int main() {
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
	BOOST_REQUIRE_CLOSE( +blas::nrm2(Y_copy - multi::array<complex, 1>{ {214.02, 0.0}, {106.43, 0.0}, {188.37, 0.0} }) , 0.0, 1e-13);
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex_value) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 2> const M_gpu = {
		{ { 9.0, 0.0}, {24.0, 0.0}, {30.0, 0.0}, {9.0, 0.0} },
		{ { 4.0, 0.0}, {10.0, 0.0}, {12.0, 0.0}, {7.0, 0.0} },
		{ {14.0, 0.0}, {16.0, 0.0}, {36.0, 0.0}, {1.0, 0.0} },
	};

	multi::thrust::cuda::array<complex, 1> const X_gpu = { {1.1, 0.0}, {2.1, 0.0}, {3.1, 0.0}, {4.1, 0.0} };

	auto const Y_gpu = +blas::gemv(/*alpha*/ 1.1, M_gpu, X_gpu);  // y = a*M*x

	multi::array<complex, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_REQUIRE_CLOSE( +blas::nrm2(Y_copy - multi::array<complex, 1>{ {209.22, 0.0}, {100.43, 0.0}, {181.17, 0.0} }) , 0.0, 1e-13);
}

BOOST_AUTO_TEST_CASE(cublas_gemv_real) {
	namespace blas = multi::blas;
	using T = double;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<T, 2> const M_gpu =  {
		{  9.0, 24.0, 30.0, 9.0 },
		{  4.0, 10.0, 12.0, 7.0 },
		{ 14.0, 16.0, 36.0, 1.0 },
	};

	multi::thrust::cuda::array<T, 1> const X_gpu = { 1.1, 2.1, 3.1, 4.1 };

	multi::thrust::cuda::array<T, 1> Y_gpu = { 4.0, 5.0, 6.0 };

	blas::gemv(/*alpha*/ 1.1, M_gpu, X_gpu, /*beta*/ 1.2, Y_gpu);  // y = a*M*x + b*y

	multi::array<T, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_REQUIRE_CLOSE( +blas::nrm2(Y_copy - multi::array<T, 1>{ 214.02, 106.43, 188.37 }) , 0.0, 1e-13);
}

return boost::report_errors();}
