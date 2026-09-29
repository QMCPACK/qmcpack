// © Alfredo A. Correa 2019-2024

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS iamax"
#define BOOST_TEST_DYN_LINK

// #include<boost/test/unit_test.hpp>

#include "../../../../adaptors/blas.hpp"
#include "../../../../adaptors/cuda.hpp"
#include "../../../../adaptors/blas/cuda.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(const multi_adaptors_blas_cuda_iamax){
	using complex = std::complex<double>; complex const I{0.0, 1.0};
	{
		multi::array<complex, 1> const A = {1.0 + 2.0*I, 2.0, 3.0 + 3.0*I, 4.0};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::array<complex, 1> const A = {1.0 + 2.0*I, 2.0, 3.0 + 3.0*I, 4.0};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
	{
		multi::cuda::managed::array<complex, 1> const A = {1.0 + 2.0*I, 2.0, 3.0 + 3.0*I, 4.0};
		using multi::blas::iamax;
		BOOST_REQUIRE( iamax(A) == 2 );
	}
}
