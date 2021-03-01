// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS scal"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/scal.hpp"

#include "../../../array.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_n){
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11. );

	blas::scal_n(2., A[2].begin(), A[2].size());
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11.*2. );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_it){
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11. );

	blas::scal(2., A[2].begin(), A[2].end());
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11.*2. );
}

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_real){
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11. );

	BOOST_REQUIRE(  blas::scal(1., A[2]) ==  A[2] );
	BOOST_REQUIRE( &blas::scal(1., A[2]) == &A[2] );
	BOOST_REQUIRE( +blas::scal(1., A[2]) ==  A[2] );

	blas::scal(2., A[2]);
	BOOST_REQUIRE( A[0][2] == 3. and A[2][2] == 11.*2. );

	BOOST_REQUIRE( &blas::scal(1., A[2]) == &A[2] );

}

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_complex_real_case){
//	using complex = std::complex<double>;
//	multi::array<complex, 2> A = {
//		{1.,  2.,  3.,  4.},
//		{5.,  6.,  7.,  8.},
//		{9., 10., 11., 12.}
//	};
//	BOOST_TEST( A[0][2] == 3. );
//	BOOST_TEST( A[2][2] == 11. );

//	blas::scal(2., A[2]); // zscal (2. is promoted to complex later)
//	BOOST_TEST( A[0][2] == 3. );
//	BOOST_REQUIRE( A[2][2] == 11.*2. );

//	blas::scal(1./2, A[2]); // zdscal
//	BOOST_TEST( A[0][2] == 3. );
//	BOOST_TEST( A[2][1] == 10. );
//	BOOST_TEST( A[2][2] == 11. );

//}

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_complex){
//	multi::array<complex, 2> A = {
//		{1. + 2.*I, 2. + 3.*I, 3. + 4.*I, 4. + 5.*I},
//		{5. + 2.*I, 6. + 3.*I, 7. + 4.*I, 8. + 5.*I},
//		{1. + 1.*I, 2. + 2.*I, 3. + 3.*I, 4. + 4.*I}
//	};
//	blas::scal(2., A[1]); // zscal (2. is promoted to complex later)
//	BOOST_TEST( A[1][2] == 14. + 8.*I );

//	blas::scal(3.*I, A[0]);
//	BOOST_TEST( A[0][1] == (2. + 3.*I)*3.*I );

//	blas::scal(2., blas::imag(A[2]));
//	assert( A[2][1] == 2. + 4.*I );
//}

////BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_cuda_noconst){
////	namespace cuda = multi::cuda;
////	cuda::array<complex, 2> A = {
////		{1. + 2.*I, 2. + 3.*I, 3. + 4.*I, 4. + 5.*I},
////		{5. + 2.*I, 6. + 3.*I, 7. + 4.*I, 8. + 5.*I},
////		{1. + 1.*I, 2. + 2.*I, 3. + 3.*I, 4. + 4.*I}
////	};
////	blas::scal(2., A[1]); // zscal (2. is promoted to complex later)
////	BOOST_REQUIRE( A[1][2] == 14. + 8.*I );

////	cuda::array<complex, 1> a = {1. + 10.*I, 2. + 20.*I, 3. + 30.*I};
////	blas::scal(2., a);
////	BOOST_REQUIRE(( a[1] == complex{4, 40} ));

//////	blas::scal(3., blas::imag(a)); // gives internal compilation error in gcc
//////	BOOST_REQUIRE(( a[1] == complex{4, 120} ));
////}

////BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_cuda_const){
////	namespace cuda = multi::cuda;
////	cuda::array<complex, 2> const A = {
////		{1. + 2.*I, 2. + 3.*I, 3. + 4.*I, 4. + 5.*I},
////		{5. + 2.*I, 6. + 3.*I, 7. + 4.*I, 8. + 5.*I},
////		{1. + 1.*I, 2. + 2.*I, 3. + 3.*I, 4. + 4.*I}
////	};
////	auto A1cpy = blas::scal(2., A[1]); // zscal (2. is promoted to complex later)
////	BOOST_REQUIRE( A1cpy[2] == 14. + 8.*I );

//////	cuda::array<complex, 1> a = {1. + 10.*I, 2. + 20.*I, 3. + 30.*I};
//////	blas::scal(2., a);
//////	BOOST_REQUIRE(( a[1] == complex{4, 40} ));

//////	blas::scal(3., blas::imag(a));
//////	BOOST_REQUIRE(( a[1] == complex{4, 120} ));
////}

//#if 0
//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_scal_cuda_managed){
//	cuda::managed::array<complex, 2> A = {
//		{1. + 2.*I, 2. + 3.*I, 3. + 4.*I, 4. + 5.*I},
//		{5. + 2.*I, 6. + 3.*I, 7. + 4.*I, 8. + 5.*I},
//		{1. + 1.*I, 2. + 2.*I, 3. + 3.*I, 4. + 4.*I}
//	};
//	using blas::scal;
//	scal(2., A[1]);
//	BOOST_REQUIRE( A[1][2] == 14. + 8.*I );

//	scal(2., blas::imag(A[1]));
//	BOOST_REQUIRE( A[1][2] == 14. + 16.*I );
//}
//#endif

