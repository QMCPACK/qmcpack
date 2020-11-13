#ifdef COMPILATION_INSTRUCTIONS
$CXXX $CXXFLAGS $0 -o $0x  `pkg-config --libs blas` -Wno-deprecated-declarations `pkg-config --cflags --libs cudart-11.0 cublas-11.0 blas` -lboost_unit_test_framework&&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS iamax"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../../adaptors/blas.hpp"
#include "../../../../adaptors/cuda.hpp"
#include "../../../../adaptors/blas/cuda.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_iamax){
	using complex = std::complex<double>; complex const I{0, 1};
	{
		multi::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::managed::array<complex, 1> const A = {1. + 2.*I, 2., 3. + 3.*I, 4.};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
}

