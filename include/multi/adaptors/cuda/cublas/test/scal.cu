// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS gemv"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/cuda/cublas.hpp>

#include <multi/adaptors/blas/axpy.hpp>
#include <multi/adaptors/blas/gemm.hpp>
#include <multi/adaptors/blas/nrm2.hpp>
#include <multi/adaptors/blas/scal.hpp>

#include <multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(cublas_scal_complex_column) {
	namespace blas = multi::blas;

	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};
	multi::thrust::cuda::array<complex, 2> arr = {
		{1.0 + I*0.0,  2.0 + I*0.0,  3.0 + I*0.0,  4.0 + I*0.0},
		{5.0 + I*0.0,  6.0 + I*0.0,  7.0 + I*0.0,  8.0 + I*0.0},
		{9.0 + I*0.0, 10.0 + I*0.0, 11.0 + I*0.0, 12.0 + I*0.0},
	};

	blas::scal(2.0, (~arr)[1]);

	multi::array<complex, 2> arr_copy = arr;

	BOOST_REQUIRE(( (~arr_copy)[1][2] == complex{20.0, 0.0} ));
}

BOOST_AUTO_TEST_CASE(cublas_scal_complex) {
	namespace blas = multi::blas;

	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};
	multi::array<complex, 1> const x_copy = { {1.0 + I*1.0}, {2.0 + I*2.0}, {3.0 + I*3.0} };
	auto x = x_copy;

	auto const alpha = complex{2.0, 3.0};
	blas::scal(alpha, x);  // x <- alpha*x

	BOOST_REQUIRE(( x[1] == alpha*x_copy[1] ));
}