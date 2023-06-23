#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS asum"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/asum.hpp"
#include "../../blas/cuda.hpp"
#include "../../../array.hpp"
#include "../../../adaptors/cuda.hpp"

#include<complex>
#include<numeric>

using std::cout;
namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_asum_double){
	multi::array<double, 2> const A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	using multi::blas::asum;
	BOOST_REQUIRE(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a + std::abs(b);}));
}

BOOST_AUTO_TEST_CASE(multi_blas_asum_complex){
	using Z = std::complex<double>; Z const I{0, 1};
	multi::array<Z, 2> const A = {
		{1. + 2.*I,  2.,  3.,  4.},
		{5.,  6. + 3.*I,  7.,  8.},
		{9., 10., 11.+ 4.*I, 12.}
	};
	using multi::blas::asum;
	BOOST_REQUIRE(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a + std::abs(real(b)) + std::abs(imag(b));}));
}

BOOST_AUTO_TEST_CASE(multi_blas_asum_double_cuda){
	multi::cuda::array<double, 2> const A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	using multi::blas::asum;
	BOOST_REQUIRE(asum(A[1]) == 26 );
}

using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(multi_blas_asum_complex_cuda){
	namespace blas = multi::blas;
	multi::cuda::array<complex, 2> const A = {
		{1. + 2.*I,  2.,  3.,  4.},
		{5.,  6. + 3.*I,  7.,  8.},
		{9., 10., 11.+ 4.*I, 12.}
	};

	BOOST_REQUIRE( blas::asum(A[1]) == 29. );
	BOOST_REQUIRE( blas::asum(A[1]({0, 4})) == 29. );
}

BOOST_AUTO_TEST_CASE(multi_blas_asum_complex_cuda_mutable){
	using Z = std::complex<double>; Z const I{0, 1};
	multi::cuda::array<Z, 2> A = {
		{1. + 2.*I,  2.,  3.,  4.},
		{5.,  6. + 3.*I,  7.,  8.},
		{9., 10., 11.+ 4.*I, 12.}
	};
	using multi::blas::asum;
	BOOST_REQUIRE( asum(A[1]) == Z{29.} );
	BOOST_REQUIRE( asum(A[1]({0, 4})) == Z{29.} );
}


