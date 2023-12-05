#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../../../adaptors/cuda.hpp" // multi::cuda ns
#include "../../../../adaptors/blas/gemm.hpp"
#include "../../../../adaptors/cuda/cublas.hpp"

#include "../../../../adaptors/thrust.hpp"
// #include "../../../complex.hpp"

#include<thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_cublas_gemm_double){
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
//	multi::thrust::cuda::array<double, 2> const a_gpu = a;
}

BOOST_AUTO_TEST_CASE(multi_cublas_gemm_complex){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
//	multi::thrust::cuda::array<complex, 2> const a_gpu = a;
}

//BOOST_AUTO_TEST_CASE(multi_cublas_gemm_thrust_complex){
//	using complex = thrust::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const a = {
//		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//	};
////	multi::thrust::cuda::array<complex, 2> const a_gpu = a;
//}

BOOST_AUTO_TEST_CASE(multi_cublas_gemm_complex2){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const a = {
//		{1. + 2.*I, 5. + 2.*I}, 
//		{9. - 1.*I, 9. + 1.*I}, 
//		{1. + 1.*I, 2. + 2.*I}
//	};
//	multi::array<complex, 2> const b = {
//		{ 11. - 2.*I, 5. + 2.*I},
//		{  7. - 3.*I, 2. + 1.*I},
//		{  8. - 1.*I, 1. + 1.*I}
//	};
////	multi::thrust::cuda::array<complex, 2> const a_gpu = a;
////	multi::thrust::cuda::array<complex, 2> const b_gpu = b;
//	namespace blas = multi::blas;
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//	//	blas::gemm(1., a, blas::H(b), 0., c);

//	//	multi::thrust::cuda::array<complex, 2> const c_gpu;
//	//	blas::gemm(1., a_gpu, b_gpu, c_gpu);
//	//	BOOST_REQUIRE( c == c_gpu );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk(1., blas::H(a), c);
//		BOOST_REQUIRE( c[2][1] == complex(41, +2) );
//		BOOST_REQUIRE( c[1][2] == complex(41, -2) );

//		multi::array<complex, 2> const c_copy = blas::herk(1., blas::H(a));
//		BOOST_REQUIRE( c_copy == c );
//	}
}

