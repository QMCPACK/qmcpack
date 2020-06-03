#ifdef COMPILATION_INSTRUCTIONS
$CXX -Wall -Wextra -Wpedantic $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS axpy"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include       "../../blas.hpp"
#include "../../../array.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_axpy){
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		auto const AC = A;
		multi::array<double, 1> const B = A[2];
		using multi::blas::axpy;
		axpy(2., B, A[1]); // daxpy
		BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
	}
	{
		using Z = std::complex<double>;
		multi::array<Z, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		auto const AC = A;
		multi::array<Z, 1> const B = A[2];
		using multi::blas::axpy;
		axpy(2., B, A[1]); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		BOOST_REQUIRE( A[1][2] == 2.*B[2] + AC[1][2] );
	}
}

