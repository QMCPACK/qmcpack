// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include<boost/test/unit_test.hpp>

#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/trsm.hpp>

#include <boost/multi/array.hpp>

#include <config.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_0x0) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A;  // NOLINT(readability-identifier-length) BLAS naming

	{
		multi::array<double, 2> B;  // NOLINT(readability-identifier-length) BLAS naming
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1.0, A, B);
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_1x1) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{10.0, },
	};
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3.0, },
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1.0, A, B);
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[0][0] , 3.0/10.0 , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1.0, A, B))[0][0] , B_cpy[0][0] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3.0, },
		};
		auto const B_cpy = B;
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 2.0, A, B);
		BOOST_REQUIRE_CLOSE( B[0][0] , 2.0*3.0/10.0 , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1.0, A, B))[0][0] , 2.*B_cpy[0][0] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3.0, 4.0, 5.0},
		};
		auto const B_cpy = B;
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1.0, A, B);
		BOOST_REQUIRE_CLOSE( B[0][0] , 3./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( B[0][1] , 4./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( B[0][2] , 5./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A, B))[0][1] , B_cpy[0][1] , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_square) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0, 3.0, 4.0},
		{NAN, 7.0, 1.0},
		{NAN, NAN, 8.0}
	};
	auto const A_cpy = triangular(blas::filling::upper, A);
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0},
			{2.0, 7.0, 1.0},
			{3.0, 4.0, 2.0}
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, B);  // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143 , 0.001 );
		BOOST_REQUIRE( (+blas::gemm(1., A_cpy, B))[1][2] == B_cpy[1][2] );
	}
	{
		auto const AT =+ ~A;
		auto const AT_cpy = triangular(blas::filling::lower, AT);
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0},
			{2.0, 7.0, 1.0},
			{3.0, 4.0, 2.0}
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., blas::T(AT), B);
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143 , 0.001 );
		BOOST_REQUIRE( (+blas::gemm(1., blas::T(AT_cpy), B))[1][2] == B_cpy[1][2] );
	}
	{
		auto const AT =+ ~A;
		auto const AT_cpy = triangular(blas::filling::lower, AT);
		multi::array<double, 2> const B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0},
			{2.0, 7.0, 1.0},
			{3.0, 4.0, 2.0}
		};
		auto BT =+ ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., blas::T(AT), blas::T(BT));
		BOOST_REQUIRE_CLOSE( blas::T(BT)[1][2] , 0.107143 , 0.001 );
		BOOST_REQUIRE( (+blas::gemm(1., blas::T(AT_cpy), blas::T(BT)))[1][2] == B[1][2] );
	}
	{
		multi::array<double, 2> const B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0},
			{2.0, 7.0, 1.0},
			{3.0, 4.0, 2.0}
		};
		auto BT =+ ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[1][2] , 0.107143 , 0.001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 2.0*I,  3.0 - 1.0*I,  4.0 + 9.0*I},
		{NAN        ,  7.0 + 4.0*I,  1.0 + 8.0*I},
		{NAN        , NAN         ,  8.0 + 2.0*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 - 9.0*I, 3.0 + 2.0*I, 4.0 + 3.0*I},
		{2.0 - 2.0*I, 7.0 - 2.0*I, 1.0 - 1.0*I},
		{3.0 + 1.0*I, 4.0 + 8.0*I, 2.0 + 7.0*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.+1.*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[1][2]) ,  2.33846   , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.0923077 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_rectangular) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1.0 + 2.0*I,  3.0 - 1.0*I,  4.0 + 9.0*I},
		{NAN         ,  7.0 + 4.0*I,  1.0 + 8.0*I},
		{NAN         , NAN         ,  8.0 + 2.0*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I, 3. + 2.*I},
		{2. - 2.*I, 7. - 2.*I},
		{3. + 1.*I, 4. + 8.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.0 + 1.0*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 2.0*I,  3.0 - 1.0*I,  4.0 + 9.0*I},
		{NAN        ,  7.0 + 4.0*I,  1.0 + 8.0*I},
		{NAN        , NAN         ,  8.0 + 2.0*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.0 + 1.0*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001);
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001);
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cpu) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imaginary unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1.0 + 2.0*I,  3.0 - 1.0*I,  4.0 + 9.0*I},
		{NAN         ,  7.0 + 4.0*I,  1.0 + 8.0*I},
		{NAN         , NAN         ,  8.0 + 2.0*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 - 9.0*I},
		{2.0 - 2.0*I},
		{3.0 + 1.0*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.0 + 1.0*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_hydrogen_inq_case_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{2.0, }, };  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 2> B = {{1.0, 2.0, 3.0}, };  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( B.size() == 1 );
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, 1.0, A, B);
		BOOST_REQUIRE( B[0][1] == B_cpy[0][1]/A[0][0] );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0},
			{2.0},
			{3.0},
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, 1.0, A, blas::T(B));
		BOOST_REQUIRE( blas::T(B)[0][1] == blas::T(B_cpy)[0][1]/A[0][0] );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_hydrogen_inq_case_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ {2.0, 0.0}, },
	};

	 {
		multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{ {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0} },
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, {1.0, 0.0}, A, B);
		BOOST_REQUIRE( B[0][1] == B_cpy[0][1]/A[0][0] );
	}
	multi::array<complex, 2> B1 = {
		{ {1.0, 0.0} },
		{ {2.0, 0.0} },
		{ {3.0, 0.0} },
	};
	multi::array<complex, 2> B2 = {
		{ {1.0, 0.0} },
		{ {2.0, 0.0} },
		{ {3.0, 0.0} },
	};

	blas::trsm(blas::side::left, blas::filling::lower, {1.0, 0.0}, A, blas::H(B1));

	{
		auto const B_cpy = B2;
		blas::trsm(blas::side::right, blas::filling::upper, {1.0, 0.0}, blas::H(A), B2);
	//  BOOST_REQUIRE( (+blas::gemm(1., A, blas::H(B)))[0][1] == blas::H(B_cpy)[0][1] );
		BOOST_REQUIRE( (+blas::gemm(1., B2, blas::H(A)))[1][0] == B_cpy[1][0] );
	}
	BOOST_REQUIRE( B1 == B2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0, 3.0, 40.0},
		{NAN, 7.0,  1.0},
		{NAN, NAN,  8.0}
	};
	auto const A_cpy = triangular(blas::filling::upper, A);
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0, 8.0},
			{2.0, 7.0, 1.0, 9.0},
			{3.0, 4.0, 2.0, 1.0},
		};
		auto const B_cpy =+ B;
		multi::array<double, 2> BT =+ ~B;
		BOOST_REQUIRE( BT == ~B );
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, B);  // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143                              , 0.001);
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, B))[1][2] , B_cpy[1][2] , 0.001);

		auto const BT_cpy = BT;
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( blas::T(BT)[1][2],  0.107143, 0.001 );

		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, blas::T(BT)))[1][2] , blas::T(BT_cpy)[1][2] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0, 3.0, 4.0, 8.0},
			{2.0, 7.0, 1.0, 9.0},
			{3.0, 4.0, 2.0, 1.0},
		};
		multi::array<double, 2> AT = ~A;
		multi::array<double, 2> BT = ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, blas::T(AT), B);  // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143 , 0.001 );

		blas::trsm(blas::side::left, blas::filling::upper, 1.0, blas::T(AT), blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[1][2] , 0.107143, 0.001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0},
			{2.0},
			{3.0},
		};
		auto const B_cpy =+ B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, B);  // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[2][0] , 0.375 , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, B))[1][0] , B_cpy[1][0] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0},
			{2.0},
			{3.0},
		};
		auto const B_cpy =+ B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.2, A, B);
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1.0, A_cpy, B))[1][0] , 1.2*B_cpy[1][0] , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1.0/1.2, A_cpy, B))[1][0] , B_cpy[1][0] , 0.00001 );
	}
	 {
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.0},
			{2.0},
			{3.0},
		};
		multi::array<double, 2> BT = rotated(B);
		blas::trsm(blas::side::left, blas::filling::upper, 1.0, A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[2][0] , 0.375 , 0.00001);
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
		{ 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
		{ 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
		{5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
	};

	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;
	blas::trsm(blas::side::left, blas::filling::upper, {1.0, 0.0}, A, blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.147059 , 0.001);
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
		{ 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
		{ 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
		{5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::left, {1.0, 0.0}, blas::U(A), blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.147059 , 0.001);
}
