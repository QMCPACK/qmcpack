// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS numeric"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "config.hpp"

#include "../../../array.hpp"
#include "../../blas/numeric.hpp"
#include "../../blas/operations.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag) {
	using complex = std::complex<double>; constexpr complex I{0, 1};

	namespace blas = multi::blas;
	multi::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
	BOOST_REQUIRE( blas::imag(a)[2] == 2. );
	BOOST_REQUIRE( blas::real(a)[2] == 9. );
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_conjugated) {
	using complex = std::complex<double>; complex const I{0, 1};

	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( B[0][0] == 1. - 3.*I );

	multi::array<complex, 2> const Bconst = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( Bconst[0][0] == 1. - 3.*I );

	namespace blas = multi::blas;
	auto BdataC = blas::make_conjugater(B.data_elements());

	decltype(blas::make_conjugater(Bconst.data_elements())) ppp;// = BdataC;
	ppp = BdataC;

	BOOST_REQUIRE( *ppp == 1. + 3.*I );

//	static_assert(    multi::blas::is_complex_array<multi::array<thrust::complex<double>, 2>>{}, "!");
	static_assert(    blas::is_complex_array<decltype(B)>{} );
	static_assert(not blas::is_conjugated<decltype(B)>{} );

	auto&& Bconj = blas::conj(B);
	static_assert( blas::is_conjugated<decltype(Bconj)>{} );

	BOOST_REQUIRE( Bconj[0][0] == 1. + 3.*I );
	BOOST_REQUIRE( imag(*base(Bconj)) == +3 );

//	BOOST_TEST_REQUIRE( base(Bconj)->imag() == +3 );
	BOOST_REQUIRE( rotated(Bconj)[1][0] == Bconj[0][1] );

//	BOOST_REQUIRE( base(Bconj) == -3.*I );
	static_assert( blas::is_complex_array<decltype(Bconj)>{} );

	BOOST_REQUIRE( blas::conj(Bconj) == B );

	BOOST_REQUIRE( blas::conj(B)[1][0] == std::conj(B[1][0]) );
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_decay) {
	using complex = std::complex<double>; complex const I{0, 1};

	multi::array<complex, 2> B = {
		{ 1. - 3.*I, 6. + 2.*I, 9. + 3.*I},
		{ 8. + 2.*I, 2. + 4.*I, 9. + 3.*I},
		{ 2. - 1.*I, 1. + 1.*I, 9. + 3.*I},
		{ 9. + 3.*I, 9. + 3.*I, 9. + 3.*I}
	};

	namespace blas = multi::blas;
	multi::array<complex, 2> conjB = blas::conj(B);

	BOOST_REQUIRE( conjB[2][1] == std::conj(B[2][1]) );
	BOOST_REQUIRE( blas::conj(B)[2][1] == std::conj(B[2][1]) );

	BOOST_REQUIRE( blas::transposed(B)[1][2] == B[2][1] );
	BOOST_REQUIRE( blas::transposed(B) == ~B );

	BOOST_REQUIRE( blas::hermitized(B)[2][1] == blas::conj(B)[1][2] );
	BOOST_REQUIRE( blas::hermitized(B)       == blas::conj(blas::transposed(B)) );

	BOOST_REQUIRE( blas::real(B)[2][1] == std::real(B[2][1]) );
	BOOST_REQUIRE( blas::imag(B)[2][1] == std::imag(B[2][1]) );

	multi::array<double, 2> B_real_doubled = {
		{ 1., -3., 6., 2., 9., 3.},
		{ 8.,  2., 2., 4., 9., 3.},
		{ 2., -1., 1., 1., 9., 3.},
		{ 9.,  3., 9., 3., 9., 3.}
	};
	BOOST_REQUIRE( sizes(blas::real_doubled(B)) == sizes(B_real_doubled) );
	BOOST_REQUIRE(       blas::real_doubled(B)  ==       B_real_doubled  );
}

#if defined(CUDA_FOUND) and CUDA_FOUND
#include<thrust/complex.h>

BOOST_AUTO_TEST_CASE(multi_blas_numeric_decay_thrust) {
	using complex = thrust::complex<double>; complex const I{0, 1};

	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	namespace blas = multi::blas;
	multi::array<complex, 2> conjB = blas::conj(B);
	BOOST_REQUIRE( conjB[1][2] == conj(B[1][2]) );
}
#endif

//#if defined(CUDA_FOUND) and CUDA_FOUND
//#include "../../blas/cuda.hpp"
//#include "../../../adaptors/cuda.hpp"
//namespace cuda = multi::cuda;

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag_cuda){
//	cuda::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
//	namespace blas = multi::blas;
//	BOOST_REQUIRE( blas::imag(a)[2] == 2. );
//}

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag_cuda_managed){
//	cuda::managed::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
//	using multi::blas::imag;
//	BOOST_REQUIRE( imag(a)[2] == 2. );
//}

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_hermitized_cuda){
//	cuda::array<complex, 2> const a = {
//		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
//		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
//		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
//	};
//	using multi::blas::hermitized;
//	hermitized(a);
//}
//#endif

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_imag_part) {
	using complex = std::complex<double>; complex const I{0., 1.};

	multi::array<double, 2> A = {
		{1., 3., 4.},
		{9., 7., 1.}
	};
	multi::array<complex, 2> Acplx = A;
	BOOST_REQUIRE( Acplx[1][1] == A[1][1] );

	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	multi::array<double, 2> Breal = {
		{1., 6.},
		{8., 2.},
		{2., 1.}
	};
	multi::array<double, 2> Bimag = {
		{-3., +2.},
		{+2., +4.},
		{-1., +1.}
	};

	using multi::blas::real;
	using multi::blas::imag;

	BOOST_REQUIRE( Breal == real(B) );
	BOOST_REQUIRE( real(B) == Breal );
	BOOST_REQUIRE( imag(B) == Bimag );

	BOOST_REQUIRE( B[1][0] == 8. + 2.*I );
	BOOST_REQUIRE( B[1][0].imag() == 2. );

	namespace blas = multi::blas;

	BOOST_REQUIRE( blas::hermitized(B)[1][2] == std::conj( B[2][1] ) );

	blas::hermitized(B)[1][2] = 20. + 30.*I;
	BOOST_REQUIRE( B[2][1] == 20. - 30.*I );
}

