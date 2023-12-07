// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS gemv"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/cuda/cublas.hpp>

#include <multi/adaptors/blas/gemm.hpp>
#include <multi/adaptors/blas/axpy.hpp>
#include <multi/adaptors/blas/nrm2.hpp>
#include <multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

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
		BOOST_REQUIRE(( arr_copy[1][0] == complex{23.0, 0.0} ));
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
		BOOST_REQUIRE(( arr_copy[1][0] == complex{23.0, 0.0} ));
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
		BOOST_REQUIRE(( arr_copy[1][0] == complex{23.0, 0.0} ));
	}
}
