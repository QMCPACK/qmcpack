// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS gemv"
// #include<boost/test/unit_test.hpp>

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

#define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )

int main() {
BOOST_AUTO_TEST_CASE(blas_axpy_complex) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;

	{
		multi::thrust::cuda::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const const_arr = arr;

		multi::thrust::cuda::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		blas::axpy(complex{2.0, 0.0}, x, arr[1]);  // arr can't be const
		
		multi::array<complex, 2> arr_copy = arr;
		BOOST_TEST(( arr_copy[1][0] == complex{23.0, 0.0} ));
	}
	{
		multi::thrust::cuda::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const const_arr = arr;

		multi::thrust::cuda::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		arr[1] += blas::axpy(complex{2.0, 0.0}, x);

		multi::array<complex, 2> arr_copy = arr;
		BOOST_TEST(( arr_copy[1][0] == complex{23.0, 0.0} ));
	}
	{
		multi::thrust::cuda::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const const_arr = arr;

		multi::thrust::cuda::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		arr[1] = blas::axpy(2.0, x);  // blas::axpy(complex{2.0, 0.0}, x);

		multi::array<complex, 2> arr_copy = arr;
		std::cout << arr_copy[1][0] << std::endl;
		BOOST_TEST(( arr_copy[1][0] == complex{23.0, 0.0} ));
	}
}

return boost::report_errors();}
