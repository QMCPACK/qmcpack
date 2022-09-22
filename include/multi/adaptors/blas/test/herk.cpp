// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS herk"
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/blas/gemm.hpp"
#include "../../../adaptors/blas/herk.hpp"
#include "../../../adaptors/blas/nrm2.hpp"

#include "../../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_herk) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; constexpr complex I{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const a = {  // NOLINT(readability-identifier-length) conventional name in BLAS
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({2, 2}, 9999.);  // NOLINT(readability-identifier-length) conventional name in BLAS
		blas::herk(a, c);
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );

		multi::array<complex, 2> const c_copy = blas::herk(1., a);
		BOOST_REQUIRE( c == c_copy );

		BOOST_REQUIRE( +blas::gemm(1., a, blas::H(a)) == blas::herk(a) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case) {
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {  // NOLINT(readability-identifier-length) conventional name in BLAS
		{0.,  1.,  2.},
		{3.,  4.,  5.},
		{6.,  7.,  8.},
		{9., 10., 11.}
	};
	BOOST_REQUIRE( (+blas::gemm(1., a, blas::T(a)))[1][2] == 86. );
	{
		multi::array<double, 2> c({4, 4});  // NOLINT(readability-identifier-length) conventional name in BLAS
		blas::herk(1.0, a, c);
		BOOST_REQUIRE( c[1][2] == (+blas::gemm(1., a, blas::T(a)))[1][2] );
	//  BOOST_REQUIRE( c[2][1] == (+blas::gemm(1., a, blas::T(a)))[2][1] );
	}
	{
		multi::array<double, 2> c = blas::herk(1.0, a);  // NOLINT(readability-identifier-length) conventional name in BLAS
		BOOST_REQUIRE( c == +blas::gemm(1., a, blas::T(a)) );
		BOOST_REQUIRE( blas::herk(a) == +blas::gemm(1., a, blas::T(a)) );
		BOOST_REQUIRE( blas::herk(2.0, a) == +blas::gemm(2.0, a, blas::T(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({2, 2}, 9999);  // NOLINT(readability-identifier-length) BLAS naming
		blas::herk(1., a, c);
//		BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case) {
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {{1., 2., 3.}};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 2> b = blas::herk(a);  // NOLINT(readability-identifier-length) BLAS naming
	
	BOOST_REQUIRE( size(b) == 1 );
	BOOST_REQUIRE( b[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case_scale) {
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {{1., 2., 3.}};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 2> b = blas::herk(0.1, a);  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE( size(b) == 1 );
	BOOST_TEST( b[0][0] == (1.*1. + 2.*2. + 3.*3.)*0.1 );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const a = { {1., 2., 3.} };  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2> b = blas::herk(1.0, a);  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE( size(b) == 1 );
	BOOST_REQUIRE( b[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case_scale, *boost::unit_test::tolerance(0.00001)) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const a = {{1., 2., 3.}};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2> b = blas::herk(0.1, a);  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE( size(b) == 1 );
	BOOST_TEST( real( b[0][0]/0.1 ) == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const a = {{1. + 2.*I, 2.+3.*I, 3. + 4.*I}};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2> b = blas::herk(a);  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE( size(b) == 1 );
	BOOST_REQUIRE( b[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(a)[0][0])) == blas::nrm2(a[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_out_param) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const a = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2> b({1, 1});  // NOLINT(readability-identifier-length) BLAS naming
	BOOST_REQUIRE( size(b) == 1 );

	blas::herk(blas::filling::upper, 1.0, blas::H(a), 0.0, b);

	BOOST_REQUIRE( b[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(b[0][0])) == blas::nrm2(blas::T(a)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized) {
	using complex = std::complex<double>; auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) conventional name in BLAS

	multi::array<complex, 2> a = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};  // NOLINT(readability-identifier-length) BLAS naming

	namespace blas = multi::blas;
	multi::array<complex, 2> b = blas::herk(blas::H(a));  // NOLINT(readability-identifier-length) BLAS naming

	BOOST_REQUIRE( size(b) == 1 );
	BOOST_REQUIRE( b[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(a))[0][0])) == blas::nrm2(rotated(a)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_auto) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> arr = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	auto arr2 = blas::herk(1., blas::hermitized(arr));
	static_assert( std::is_same<decltype(arr2), multi::array<complex, 2>>{}, "!" );
	BOOST_REQUIRE( size(arr2) == 1 );
	BOOST_REQUIRE( arr2[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(arr))[0][0])) == blas::nrm2(rotated(arr)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_identity) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const arr = {  // NOLINT(readability-identifier-length) : conventional one-letter operation BLASs
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};

	{
		multi::array<complex, 2> arr2({2, 2}, 9999.);  // NOLINT(readability-identifier-length) conventional one-letter operation BLASs
		blas::herk(blas::filling::lower, 1., arr, 0., arr2); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( arr2[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( arr2[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);  // NOLINT(readability-identifier-length) conventional one-letter operation BLASs
		static_assert(blas::is_conjugated<decltype(blas::H(c))>{}, "!" );

		blas::herk(blas::filling::lower, 1., arr, 0., blas::H(c));  // c†=c=aa†=(aa†)†, `c` in upper triangular

		BOOST_REQUIRE( blas::H(c)[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( blas::H(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);  // NOLINT(readability-identifier-length) : conventional one-letter operation BLASs
		herk(blas::filling::lower, 1., blas::T(arr), 0., blas::T(c));  // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( transposed(c)[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);  // NOLINT(readability-identifier-length) : conventional one-letter operation BLASs
		blas::herk(blas::filling::lower, 1., blas::T(arr), 0., blas::H(blas::T(c)));  // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( blas::H(blas::T(c))[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( blas::H(blas::T(c))[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
		blas::herk(blas::U, 1., arr, 0., c);  // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
		BOOST_REQUIRE( c[1][0] == 9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
		blas::herk(1., arr, c); // c†=c=aa†=(aa†)†
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);  // NOLINT(readability-identifier-length) : conventional one-letter operation BLAS
		blas::herk(blas::L, 1., blas::H(arr), 0., c);  // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0] == complex(52., 90.) );
		BOOST_REQUIRE( c[0][1] == 9999. );
	}
}
