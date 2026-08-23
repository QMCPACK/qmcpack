// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS gemv"
// #include<boost/test/unit_test.hpp>

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>
#include <boost/multi/adaptors/blas/scal.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

// #define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )

int main() {
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

	BOOST_TEST(( (~arr_copy)[1][2] == complex{20.0, 0.0} ));
}

BOOST_AUTO_TEST_CASE(cublas_scal_complex) {
	namespace blas = multi::blas;

	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};
	multi::array<complex, 1> const x_copy = { {1.0 + I*1.0}, {2.0 + I*2.0}, {3.0 + I*3.0} };
	auto x = x_copy;

	auto const alpha = complex{2.0, 3.0};
	blas::scal(alpha, x);  // x <- alpha*x

	BOOST_TEST(( x[1] == alpha*x_copy[1] ));
}

return boost::report_errors();}
