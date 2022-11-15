// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS numeric"
#include<boost/test/unit_test.hpp>

#include "../../../array.hpp"
#include "../../blas/numeric.hpp"
#include "../../blas/operations.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag) {
	using complex = std::complex<double>; constexpr complex I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	namespace blas = multi::blas;
	multi::array<complex, 1> array = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
	BOOST_REQUIRE( blas::imag(array)[2] == 2. );
	BOOST_REQUIRE( blas::real(array)[2] == 9. );
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_conjugated) {
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> array = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( array[0][0] == 1. - 3.*I );

	multi::array<complex, 2> const carray = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	BOOST_REQUIRE( carray[0][0] == 1. - 3.*I );

	namespace blas = multi::blas;
	auto conjr = blas::make_conjugater(array.data_elements());

	decltype(blas::make_conjugater(carray.data_elements())) ppp;// = BdataC;
	ppp = conjr;

	BOOST_REQUIRE( *ppp == 1. + 3.*I );

//	static_assert(    multi::blas::is_complex_array<multi::array<thrust::complex<double>, 2>>{}, "!");
	static_assert(    blas::is_complex_array<decltype(array)>{} );
	static_assert(not blas::is_conjugated<decltype(array)>{} );

	auto&& conjd_array = blas::conj(array);
	static_assert( blas::is_conjugated<decltype(conjd_array)>{} );

	BOOST_REQUIRE( conjd_array[0][0] == 1. + 3.*I );
	BOOST_REQUIRE( imag(*base(conjd_array)) == +3 );

//	BOOST_TEST_REQUIRE( base(Bconj)->imag() == +3 );
	BOOST_REQUIRE( rotated(conjd_array)[1][0] == conjd_array[0][1] );

//	BOOST_REQUIRE( base(Bconj) == -3.*I );
	static_assert( blas::is_complex_array<decltype(conjd_array)>{} );

	BOOST_REQUIRE( blas::conj(conjd_array) == array );

	BOOST_REQUIRE( blas::conj(array)[1][0] == std::conj(array[1][0]) );
}

BOOST_AUTO_TEST_CASE(multi_blas_numeric_decay) {
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> arr = {
		{ 1. - 3.*I, 6. + 2.*I, 9. + 3.*I},
		{ 8. + 2.*I, 2. + 4.*I, 9. + 3.*I},
		{ 2. - 1.*I, 1. + 1.*I, 9. + 3.*I},
		{ 9. + 3.*I, 9. + 3.*I, 9. + 3.*I}
	};

	namespace blas = multi::blas;
	multi::array<complex, 2> conj_arr = blas::conj(arr);

	BOOST_REQUIRE( conj_arr[2][1] == std::conj(arr[2][1]) );
	BOOST_REQUIRE( blas::conj(arr)[2][1] == std::conj(arr[2][1]) );

	BOOST_REQUIRE( blas::transposed(arr)[1][2] == arr[2][1] );
	BOOST_REQUIRE( blas::transposed(arr) == ~arr );

	BOOST_REQUIRE( blas::hermitized(arr)[2][1] == blas::conj(arr)[1][2] );
	BOOST_REQUIRE( blas::hermitized(arr)       == blas::conj(blas::transposed(arr)) );

	BOOST_REQUIRE( blas::real(arr)[2][1] == std::real(arr[2][1]) );
	BOOST_REQUIRE( blas::imag(arr)[2][1] == std::imag(arr[2][1]) );

	multi::array<double, 2> B_real_doubled = {
		{ 1., -3., 6., 2., 9., 3.},
		{ 8.,  2., 2., 4., 9., 3.},
		{ 2., -1., 1., 1., 9., 3.},
		{ 9.,  3., 9., 3., 9., 3.}
	};
	BOOST_REQUIRE( sizes(blas::real_doubled(arr)) == sizes(B_real_doubled) );
	BOOST_REQUIRE(       blas::real_doubled(arr)  ==       B_real_doubled  );
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

BOOST_AUTO_TEST_CASE(multi_blas_numeric_real_imag_part) {
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<double, 2> arr = {
		{1., 3., 4.},
		{9., 7., 1.}
	};
	multi::array<complex, 2> complex_arr = arr;
	BOOST_REQUIRE( complex_arr[1][1] == arr[1][1] );

	multi::array<complex, 2> arr2 = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	multi::array<double, 2> arr2_real = {
		{1., 6.},
		{8., 2.},
		{2., 1.}
	};
	multi::array<double, 2> arr2_imag = {
		{-3., +2.},
		{+2., +4.},
		{-1., +1.}
	};

	using multi::blas::real;
	using multi::blas::imag;

	BOOST_REQUIRE( arr2_real == real(arr2) );
	BOOST_REQUIRE( real(arr2) == arr2_real );
	BOOST_REQUIRE( imag(arr2) == arr2_imag );

	BOOST_REQUIRE( arr2[1][0] == 8. + 2.*I );
	BOOST_REQUIRE( arr2[1][0].imag() == 2. );

	namespace blas = multi::blas;

	BOOST_REQUIRE( blas::hermitized(arr2)[1][2] == std::conj( arr2[2][1] ) );

	blas::hermitized(arr2)[1][2] = 20. + 30.*I;
	BOOST_REQUIRE( arr2[2][1] == 20. - 30.*I );
}
