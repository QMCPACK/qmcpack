#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS/cuBLAS iamax"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/iamax.hpp"

#include "../../../array.hpp"
#include "../../../adaptors/cuda.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include<complex>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax){
	multi::array<complex, 2> const A = {
		{1. + 2.*I,  2.,  3.,  4.},
		{5.,  6. + 3.*I,  7.,  8.},
		{9., 10., 11.+ 4.*I, 12.}
	};
	using blas::iamax;
	auto chess = [](auto const& a, auto const& b){
		using std::abs; 
		return abs(real(a))+abs(imag(a)) < abs(real(b))+abs(imag(b));
	};
	BOOST_REQUIRE(iamax(A[1])==std::max_element(begin(A[1]), end(A[1]), chess)-begin(A[1]));
	BOOST_REQUIRE(A[1][iamax(A[1])]==*std::max_element(begin(A[1]), end(A[1]), chess));
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax_cuda){
	multi::cuda::array<complex, 2> const A = {
		{1. + 2.*I,  2.       ,  3.      ,  4.},
		{5.       ,  6. + 3.*I,  7.      ,  8.},
		{9.       , 10.       , 11.+ 4.*I, 12.}
	};
	using blas::iamax;
	BOOST_REQUIRE(iamax(A[1])==1);
}


