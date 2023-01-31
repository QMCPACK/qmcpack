#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXX $0 -o $0x -lcudart -lcublas `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS trsm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../../memory/adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/blas/gemm.hpp"
#include "../../../adaptors/blas/trsm.hpp"
//#include "../../../adaptors/blas/cuda.hpp"

//#include "../../../adaptors/cuda.hpp"
#include "../../../array.hpp"

#include <config.hpp>

namespace multi = boost::multi;

template<class Matrix>
auto triangular(multi::blas::filling f, Matrix const& m) {  // NOLINT(readability-identifier-length) BLAS naming
	auto ret =+ m;
	switch(f) {
	case multi::blas::filling::upper:
		for(multi::size_type i = 0; i != size( ret); ++i) {
			for(multi::size_type j = 0; j != std::min(i, size(~ret)); ++j) {
				ret[i][j] = 0.;
			}
		}
		break;
	case multi::blas::filling::lower:
		for(multi::size_type j = 0; j != size(~ret); ++j) {
			for(multi::size_type i = 0; i != std::min(j, size( ret)); ++i) {
				ret[i][j] = 0.;
			}
		}
		break;
	}
	return ret;
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_0x0) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A;  // NOLINT(readability-identifier-length) BLAS naming

	{
		multi::array<double, 2> B;  // NOLINT(readability-identifier-length) BLAS naming
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1., A, B);
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_1x1) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{10., },
	};
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3., },
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1., A, B);
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[0][0] , 3./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A, B))[0][0] , B_cpy[0][0] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3., },
		};
		auto const B_cpy = B;
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 2., A, B);
		BOOST_REQUIRE_CLOSE( B[0][0] , 2.*3./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A, B))[0][0] , 2.*B_cpy[0][0] , 0.00001 );
	}
	 {
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{3., 4., 5.},
		};
		auto const B_cpy = B;
		// B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 1., A, B);
		BOOST_REQUIRE_CLOSE( B[0][0] , 3./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( B[0][1] , 4./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( B[0][2] , 5./10. , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A, B))[0][1] , B_cpy[0][1] , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_square) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{     1.,      3.,  4.},
		{    NAN,      7.,  1.},
		{    NAN,     NAN,  8.}
	};
	auto const A_cpy = triangular(blas::filling::upper, A);
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, B);  // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143 , 0.001 );
		BOOST_REQUIRE( (+blas::gemm(1., A_cpy, B))[1][2] == B_cpy[1][2] );
	}
	{
		auto const AT =+ ~A;
		auto const AT_cpy = triangular(blas::filling::lower, AT);
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
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
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		auto BT =+ ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., blas::T(AT), blas::T(BT));
		BOOST_REQUIRE_CLOSE( blas::T(BT)[1][2] , 0.107143 , 0.001 );
		BOOST_REQUIRE( (+blas::gemm(1., blas::T(AT_cpy), blas::T(BT)))[1][2] == B[1][2] );
	}
	{
		multi::array<double, 2> const B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		auto BT =+ ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[1][2] , 0.107143 , 0.001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.+1.*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[1][2]) ,  2.33846   , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.0923077 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_rectangular) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I, 3. + 2.*I},
		{2. - 2.*I, 7. - 2.*I},
		{3. + 1.*I, 4. + 8.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.+1.*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.+1.*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001);
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001);
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cpu) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imaginary unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	blas::trsm(blas::side::left, blas::filling::lower, 2.+1.*I, blas::H(A), B);  // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE_CLOSE( real(B[2][0]) , -4.16471 , 0.0001 );
	BOOST_REQUIRE_CLOSE( imag(B[2][0]) ,  8.25882 , 0.0001 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_hydrogen_inq_case_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{2., }, };  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 2> B = {{1., 2., 3.}, };  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( B.size() == 1 );
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, 1., A, B);
		BOOST_REQUIRE( B[0][1] == B_cpy[0][1]/A[0][0] );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.},
			{2.},
			{3.},
		};
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, 1., A, blas::T(B));
		BOOST_REQUIRE( blas::T(B)[0][1] == blas::T(B_cpy)[0][1]/A[0][0] );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_hydrogen_inq_case_complex) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = {{2., }, };  // NOLINT(readability-identifier-length) BLAS naming

	 {
		multi::array<complex, 2> B = {{1., 2., 3.}, };  // NOLINT(readability-identifier-length) BLAS naming
		auto const B_cpy = B;
		blas::trsm(blas::side::left, blas::filling::lower, 1., A, B);
		BOOST_REQUIRE( B[0][1] == B_cpy[0][1]/A[0][0] );
	}
	multi::array<complex, 2> B1 = {
		{1.},
		{2.},
		{3.},
	};
	multi::array<complex, 2> B2 = {
		{1.},
		{2.},
		{3.},
	};
	{
	//	auto const B_cpy = B1;
		blas::trsm(blas::side::left, blas::filling::lower, 1., A, blas::H(B1));
	//	BOOST_REQUIRE( (+blas::gemm(1., A, blas::H(B1)))[0][1] == blas::H(B_cpy)[0][1] );
	}

	 {
		auto const B_cpy = B2;
		blas::trsm(blas::side::right, blas::filling::upper, 1., blas::H(A), B2);
	//	BOOST_REQUIRE( (+blas::gemm(1., A, blas::H(B)))[0][1] == blas::H(B_cpy)[0][1] );
		BOOST_REQUIRE( (+blas::gemm(1., B2, blas::H(A)))[1][0] == B_cpy[1][0] );
	}
	BOOST_REQUIRE( B1 == B2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{  1.,   3.,  4.},
		{ NAN,   7.,  1.},
		{ NAN,  NAN,  8.}
	};
	auto const A_cpy = triangular(blas::filling::upper, A);
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1., 3., 4., 8.},
			{2., 7., 1., 9.},
			{3., 4., 2., 1.},
		};
		auto const B_cpy =+ B;
		multi::array<double, 2> BT =+ ~B;
		BOOST_REQUIRE( BT == ~B );
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143                              , 0.001);
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, B))[1][2] , B_cpy[1][2] , 0.001);

		auto const BT_cpy = BT;
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( blas::T(BT)[1][2],  0.107143, 0.001 );

		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, blas::T(BT)))[1][2] , blas::T(BT_cpy)[1][2] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1., 3., 4., 8.},
			{2., 7., 1., 9.},
			{3., 4., 2., 1.},
		};
		multi::array<double, 2> AT = ~A;
		multi::array<double, 2> BT = ~B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., blas::T(AT), B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[1][2] , 0.107143 , 0.001 );

		blas::trsm(blas::side::left, blas::filling::upper, 1., blas::T(AT), blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[1][2] , 0.107143, 0.001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.},
			{2.},
			{3.},
		};
		auto const B_cpy =+ B;
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_REQUIRE_CLOSE( B[2][0] , 0.375 , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, B))[1][0] , B_cpy[1][0] , 0.00001 );
	}
	{
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.},
			{2.},
			{3.},
		};
		auto const B_cpy =+ B;
		blas::trsm(blas::side::left, blas::filling::upper, 1.2, A, B);
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1., A_cpy, B))[1][0] , 1.2*B_cpy[1][0] , 0.00001 );
		BOOST_REQUIRE_CLOSE( (+blas::gemm(1./1.2, A_cpy, B))[1][0] , B_cpy[1][0] , 0.00001 );
	}
	 {
		multi::array<double, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
			{1.},
			{2.},
			{3.},
		};
		multi::array<double, 2> BT = rotated(B);
		blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::T(BT));
		BOOST_REQUIRE_CLOSE( (~BT)[2][0] , 0.375 , 0.00001);
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{ 1. + 4.*I,  3.      ,  4.- 10.*I},
		{ 0.       ,  7.- 3.*I,  1.       },
		{ 0.       ,  0.      ,  8.-  2.*I}
	};
	multi::array<complex, 2> B = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
		{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
	};

	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;
	blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.147059 , 0.001);
}

#if 0
BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check) {
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.,  4.- 10.*I},
		{ 0.,  7.- 3.*I,  1.},
		{ 0.,  0.,  8.- 2.*I}
	};
	namespace blas = multi::blas;

	 {
		{
			multi::array<complex, 2> B = {
				{1. + 1.*I, 5. + 3.*I},
				{2. + 1.*I, 9. + 3.*I},
				{3. + 1.*I, 1. - 1.*I},
			};
			auto S = blas::trsm(blas::side::left, blas::filling::lower, 1., blas::H(A), B); // S = A⁻¹†.B, S† = B†.A⁻¹
			BOOST_REQUIRE_CLOSE( real(S[2][1]) , 1.71608 , 0.001 );
		}
		{
			multi::array<complex, 2> B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_REQUIRE_CLOSE( imag(S[2][1]) , +0.147059 , 0.001);
			BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.147059 , 0.001);
		}
		 {
			multi::array<complex, 2> B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 2., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_REQUIRE_CLOSE( imag(S[2][1]) , +0.147059*2. , 0.001 );
			BOOST_REQUIRE_CLOSE( imag(B[1][2]) , -0.147059*2. , 0.001 );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_1x1_check) {
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {
		{ 4.},
	};
	 {
		{
			multi::array<double, 2> B = {
				{5.},
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::general, 3., A, B);
			BOOST_REQUIRE( S[0][0] == 3.*5./4. );
		}
		{
			multi::array<double, 2> B = {
				{5.},
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 1., A, B);
			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
		}
		 {
			multi::array<double, 2> B = {
				{5.},
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 1., A, B);
			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_1x1_check) {
	using complex = std::complex<double>; complex const I = complex{0, 1};
	multi::array<complex, 2> const A = {
		{ 4. + 2.*I},
	};
	namespace blas = multi::blas;
	 {
		multi::array<complex, 2> B = {
			{5. + 1.*I},
		};
		auto const B_cpy =+ B;

		blas::trsm(blas::side::left, blas::filling::upper, 3.+5.*I, A, B);
		BOOST_REQUIRE_CLOSE( real((+blas::gemm(1., A, B))[0][0]) , real((3.+5.*I)*B_cpy[0][0]) , 0.00001 );
		BOOST_REQUIRE_CLOSE( imag((+blas::gemm(1., A, B))[0][0]) , imag((3.+5.*I)*B_cpy[0][0]) , 0.00001 );

		BOOST_REQUIRE_CLOSE( real((+blas::gemm(1./(3.+5.*I), A, B))[0][0]) , real(B_cpy[0][0]) , 0.00001 );
		BOOST_REQUIRE_CLOSE( imag((+blas::gemm(1./(3.+5.*I), A, B))[0][0]) , imag(B_cpy[0][0]) , 0.00001 );
	}
}

#if defined(CUDA_FOUND) and CUDA_FOUND
#include<thrust/complex.h>

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_thrust_nonsquare_default_diagonal_hermitized_gemm_check) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.      ,  4.- 10.*I},
		{ 0.       ,  7.- 3.*I,  1.       },
		{ 0.       ,  0.      ,  8.-  2.*I}
	};
	 {
		{
			multi::array<complex, 2> B = {
				{1. + 1.*I, 5. + 3.*I},
				{2. + 1.*I, 9. + 3.*I},
				{3. + 1.*I, 1. - 1.*I},
			};
			auto S = blas::trsm(blas::side::left, blas::filling::lower, 1., blas::H(A), B); // S = A⁻¹†.B, S† = B†.A⁻¹
			BOOST_REQUIRE_CLOSE( S[2][1].real() , 1.71608 , 0.001 );
			BOOST_REQUIRE( S == B );
		}
		{
			multi::array<complex, 2> B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_REQUIRE_CLOSE( B[1][2].imag() , -0.147059 , 0.001 );
			BOOST_REQUIRE( S == blas::H(B) );
		}
		{
			multi::array<complex, 2> B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =+ blas::trsm(blas::side::left, blas::filling::upper, 2., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_REQUIRE_CLOSE( B[1][2].imag() , -0.147059*2. , 0.001  );
			BOOST_REQUIRE( S == blas::H(B) );
		}
	}
}
//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cuda, *utf::tolerance(0.00001)){
//	namespace cuda = multi::cuda;
//	cuda::array<complex, 2> A = {
//		{ 1.,  3.,  4.},
//		{NAN,  7.,  1.},
//		{NAN, NAN,  8.}
//	};
////	multi::cuda::array<complex, 2> const B = {
////		{1.},
////		{2.},
////		{3.}
////	};
//	namespace blas = multi::blas;
////	auto Bcpy = blas::trsm(blas::filling::upper, 1., A, B); // B ⬅ α Inv[A].B, B† ⬅ B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
////	multi::array<complex, 2> Bcpu = Bcpy;
////	BOOST_TEST_REQUIRE( std::real(Bcpu[2][0]) == 0.375 );
////	BOOST_TEST_REQUIRE( std::imag(Bcpu[2][0]) == 0.    );
//}
#endif

#endif
#if 0

//template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_column_cuda, *utf::tolerance(0.00001)) {
	multi::cuda::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	multi::cuda::array<double, 2> B = {
		{1.},
		{2.},
		{3.}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	trsm(filling::upper, 1., A, B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_REQUIRE( B[2][0] == 0.375 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cuda2, *utf::tolerance(0.00001)) {
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::array<complex, 2> B = {
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	multi::array<complex, 2> Bcpu = B;
	BOOST_REQUIRE( real(Bcpu[2][0]) == -4.16471 );
	BOOST_REQUIRE( imag(Bcpu[2][0]) ==  8.25882 );
}

BOOST_AUTO_TEST_CASE(multi_blas_cuda_trsm_complex, *utf::tolerance(0.00001)) {
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::array<complex, 2> const B = {
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};

	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
//	auto C = trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	auto C = trsm(filling::lower, 1., hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit
}

BOOST_AUTO_TEST_CASE(multi_blas_cuda_managed_trsm_complex, *utf::tolerance(0.00001)) {
	multi::cuda::managed::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::managed::array<complex, 2> const B = {
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};

	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	auto C = trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
}
#endif
