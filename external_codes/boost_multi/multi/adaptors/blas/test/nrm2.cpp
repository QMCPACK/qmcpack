#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS nrm2"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas.hpp"
#include "../../../array.hpp"
#include "../../../adaptors/cuda.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include<complex>

namespace multi = boost::multi;

using complex = std::complex<double>; constexpr complex I{0,1};

BOOST_AUTO_TEST_CASE(multi_blas_nrm2){
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( blas::nrm2(A[1]) == std::sqrt(blas::dot(A[1], A[1])) );

	{
		multi::array<complex, 1> A = {1.+I, 3.+2.*I, 3.+4.*I};
		BOOST_REQUIRE( blas::dot(A, A)() == (1.+I)*(1.+I) + (3.+2.*I)*(3.+2.*I) + (3.+4.*I)*(3.+4.*I) );
	}
	{
		multi::cuda::array<double, 2> const Agpu = A;
		multi::cuda::static_array<double, 0> n = 1.2;
		blas::nrm2(Agpu[1], n);
	}
	{
		multi::cuda::array<double, 2> Agpu = A;
		double n = 99.;
		blas::nrm2(Agpu[1], n); // cuda supports putting scalar results in CPU
		double n2{blas::nrm2(Agpu[1])};
		BOOST_REQUIRE( n == n2 );
	}
}

