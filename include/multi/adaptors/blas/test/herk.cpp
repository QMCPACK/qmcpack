// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/blas/gemm.hpp"
#include "../../../adaptors/blas/herk.hpp"
#include "../../../adaptors/blas/nrm2.hpp"

#include "../../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_herk) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; constexpr complex I{0, 1};

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
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
	multi::array<double, 2> const a = {
		{0.,  1.,  2.},
		{3.,  4.,  5.},
		{6.,  7.,  8.},
		{9., 10., 11.}
	};
	BOOST_REQUIRE( (+blas::gemm(1., a, blas::T(a)))[1][2] == 86. );
	{
		multi::array<double, 2> c({4, 4});
		blas::herk(1.0, a, c);
		BOOST_REQUIRE( c[1][2] == (+blas::gemm(1., a, blas::T(a)))[1][2] );
	//  BOOST_REQUIRE( c[2][1] == (+blas::gemm(1., a, blas::T(a)))[2][1] );
	}
	{
		multi::array<double, 2> c = blas::herk(1.0, a);
		BOOST_REQUIRE( c == +blas::gemm(1., a, blas::T(a)) );
		BOOST_REQUIRE( blas::herk(a) == +blas::gemm(1., a, blas::T(a)) );
		BOOST_REQUIRE( blas::herk(2.0, a) == +blas::gemm(2.0, a, blas::T(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({2, 2}, 9999);
		blas::herk(1., a, c);
//		BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case_scale) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( B[0][0] == (1.*1. + 2.*2. + 3.*3.)*0.1 );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = { {1., 2., 3.} };
	multi::array<complex, 2> B = blas::herk(1.0, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case_scale, *boost::unit_test::tolerance(0.00001)) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = {{1., 2., 3.}};
	multi::array<complex, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( real( B[0][0]/0.1 ) == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const A = {{1. + 2.*I, 2.+3.*I, 3. + 4.*I}};
	multi::array<complex, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(A)[0][0])) == blas::nrm2(A[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_out_param) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	multi::array<complex, 2> B({1, 1});
	BOOST_REQUIRE( size(B) == 1 );

	blas::herk(blas::filling::upper, 1.0, blas::H(A), 0.0, B);

	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(B[0][0])) == blas::nrm2(blas::T(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized) {
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};

	namespace blas = multi::blas;
	multi::array<complex, 2> B = blas::herk(blas::H(A));

	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_auto) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	auto B = blas::herk(1., blas::hermitized(A));
	static_assert( std::is_same<decltype(B), multi::array<complex, 2>>{}, "!" );
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_identity) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};

	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::filling::lower, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		static_assert(blas::is_conjugated<decltype(blas::H(c))>{}, "!" );

		blas::herk(blas::filling::lower, 1., a, 0., blas::H(c)); // c†=c=aa†=(aa†)†, `c` in upper triangular

		BOOST_REQUIRE( blas::H(c)[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( blas::H(c)[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({2, 2}, 9999.);
	//	blas::herk(blas::filling::lower, 1., a, 0., blas::T(c)); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( transposed(c)[1][0]==complex(50., -49.) );
	//	BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aT(aT)† not supported
	//	print(c);
	//	BOOST_REQUIRE( c[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., hermitized(c)); // c†=c=aT(aT)† not supported
	//	BOOST_REQUIRE( hermitized(c)[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( hermitized(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(blas::filling::lower, 1., blas::T(a), 0., blas::T(c)); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( transposed(c)[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::filling::lower, 1., blas::T(a), 0., blas::H(blas::T(c))); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( blas::H(blas::T(c))[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( blas::H(blas::T(c))[0][1]==9999. );
	}
	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		using namespace multi::blas;
//		blas::herk(blas::filling::lower, 1., blas::T(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
//		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::U, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(1., a, c); // c†=c=aa†=(aa†)†
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::L, 1., blas::H(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({3, 3}, 9999.);
	//	using namespace multi::blas;
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( c[0][1]==9999. );
	//	BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
	}
}
