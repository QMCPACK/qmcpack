#ifdef COMPILATION_INSTRUCTIONS
$CXX -Wall -Wextra -Wpedantic $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas.hpp"
#include "../../../array.hpp"
#include "../../../adaptors/cuda.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include<complex>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

template<class T> void what(T&&) = delete;

using complex = std::complex<double>; complex const I{0,1};

BOOST_AUTO_TEST_CASE(multi_blas_nrm2){
	multi::array<double, 2> const A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	{
		using multi::blas::dot;
		using multi::blas::nrm2;
	//	BOOST_REQUIRE( nrm2(A[1]) == std::sqrt(dot(A[1], A[1])) );
	}
	{

		multi::array<complex, 1> A = {1.+I, 3.+2.*I, 3.+4.*I};
		using blas::dot;
		using blas::hermitized;
		using blas::nrm2;
		auto n = nrm2(A);
		using std::real;
		BOOST_TEST( dot(A, A)() == (1.+I)*(1.+I) + (3.+2.*I)*(3.+2.*I) + (3.+4.*I)*(3.+4.*I) );
	}
	{
	//	multi::cuda::array<double, 2> Agpu = A;
	//	multi::cuda::static_array<double, 0> n = 1.2;
	//	blas::nrm2(Agpu[1], n);
	}
	{
		multi::cuda::array<double, 2> Agpu = A;
		double n = 99.;
		blas::nrm2(Agpu[1], n); // cuda supports putting scalar results in CPU
	}
}

