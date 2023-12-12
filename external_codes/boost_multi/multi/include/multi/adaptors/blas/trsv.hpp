#ifndef MULTI_ADAPTORS_BLAS_TRSV_HPP  // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
#define MULTI_ADAPTORS_BLAS_TRSV_HPP
// Copyright 2019-2021 Alfredo A. Correa

#include "../blas/core.hpp"

#include "../blas/operations.hpp" // uplo
#include "../blas/filling.hpp"
#include "../blas/side.hpp"

#include "../../config/NODISCARD.hpp"

namespace boost::multi::blas {

//enum DIAG : char{U='U', N='N'};

enum class diagonal : char {//typename std::underlying_type<char>::type{
	unit = 'U',
	non_unit = 'N', general = non_unit
};

using core::trsv;

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0> 
auto trsv_base(A&& a) {return base(a);}

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0> 
auto trsv_base(A&& a) {return underlying(base(a));}

template<class A2D, class X1D>
auto trsv(filling a_nonzero_side, diagonal a_diag, A2D const& a, X1D&& x)
->decltype(trsv(static_cast<char>(flip(a_nonzero_side)), 'N', static_cast<char>(a_diag), size(x), trsv_base(a), stride(rotated(a)), trsv_base(x), stride(x)), std::forward<X1D>(x))
{
//	if(is_conjugated(x)) trsv(a_nonzero_side, a_diag, conjugated(a), conjugated(std::forward<X1D>(x)));
	{
		auto base_a = trsv_base(a);
		auto base_x = trsv_base(x);
		if(not is_conjugated<A2D>{}) {
				 if(stride(        a )==1) {trsv(static_cast<char>(flip(a_nonzero_side)), 'N', static_cast<char>(a_diag), size(x), base_a, stride(rotated(a)), base_x, stride(x));}
			else if(stride(rotated(a))==1) {trsv(static_cast<char>(     a_nonzero_side ), 'T', static_cast<char>(a_diag), size(x), base_a, stride(        a ), base_x, stride(x));}
			else                           {assert(0);}
		}else{
				 if(stride(        a )==1) {assert(0);} //TODO fallback to trsm?
			else if(stride(rotated(a))==1) {trsv(static_cast<char>(     a_nonzero_side ), 'C', static_cast<char>(a_diag), size(x), base_a, stride(        a ), base_x, stride(x));}
			else                           {assert(0);}
		}
	}
	return std::forward<X1D>(x);
}

template<class A2D, class X1D>
auto trsv(filling a_nonzero_side, A2D const& a, X1D&& x)
->decltype(trsv(a_nonzero_side, diagonal::general, a, std::forward<X1D>(x))) {
	return trsv(a_nonzero_side, diagonal::general, a, std::forward<X1D>(x)); }

#if 0


#if 1
template<class A2D, class X1D, class Ret = typename X1D::decay_type>
Ret trsv(filling a_nonzero_side, diagonal a_diag, A2D const& a, X1D const& x, void* = 0){
	return trsv(a_nonzero_side, a_diag, a, Ret{x});}

template<class A2D, class X1D, class Ret = typename X1D::decay_type>
Ret trsv(filling a_nonzero_side, A2D const& a, X1D const& x, void* = 0){
	return trsv(a_nonzero_side, a, Ret{x});}
#endif
#endif

}  // end namespace boost::multi::blas

//#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_TRSV

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi.BLAS trsv"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../blas/gemm.hpp"

//#include "../../array.hpp"

//#include<iostream>

//namespace multi = boost::multi;

//template<class M> decltype(auto) print_1D(M const& C){
//	using boost::multi::size; using std::cout;
//	for(int i = 0; i != size(C); ++i)
//		cout<< C[i] <<' ';
//	cout<<std::endl;
//}

//template<class M> decltype(auto) print(M const& C){
//	using boost::multi::size; using std::cout;
//	for(int i = 0; i != size(C); ++i){
//		for(int j = 0; j != size(C[i]); ++j)
//			cout<< C[i][j] <<' ';
//		cout<<std::endl;
//	}
//	return cout<<std::endl;
//}

//namespace utf = boost::unit_test;
//namespace blas = multi::blas;

//BOOST_AUTO_TEST_CASE(multi_blas_trsv_real_square, *utf::tolerance(0.0001)){


//	{
//		multi::array<double, 2> const A = {
//			{ 1.,  3.,  4.},
//			{ NAN,  7.,  1.},
//			{ NAN,  NAN,  8.}
//		};
//		multi::array<double, 1> b = {1., 3., 4.};
//		blas::trsv(blas::filling::upper, blas::diagonal::general, A, b); // B<-Solve(A.X==B), B<-A⁻¹.B, B⊤<-(A⁻¹.B)⊤, B<-B⊤.A⁻¹⊤
//		BOOST_TEST( b[0] == -2.07143 );
//		BOOST_TEST( b[1] ==  0.357143 );
//		BOOST_TEST( b[2] ==  0.5 );
//	}
//	{
//		multi::array<double, 2> const A = {
//			{ 1.,  3.,  4.},
//			{ NAN,  7.,  1.},
//			{ NAN,  NAN,  8.}
//		};
//		multi::array<double, 1> b = {1., 3., 4.};
//		blas::trsv(blas::filling::lower, blas::diagonal::general, blas::T(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//		BOOST_TEST( b[0] ==  1. );
//		BOOST_TEST( b[1] ==  0. );
//		BOOST_TEST( b[2] ==  0. );
//	}
//#if 0
//	{
//		multi::array<double, 1> b = {3., 3., 1.};
//	//	trsv(filling::lower, diagonal::general, hermitized(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//	//	BOOST_TEST( b[0] ==  3. );
//	//	BOOST_TEST( b[1] == -0.857143 );
//	//	BOOST_TEST( b[2] == -1.26786 );
//	}
//#endif
//}

//#if 0
//using complex = std::complex<double>;

//BOOST_AUTO_TEST_CASE(multi_blas_trsv_complex_real_case_square, *utf::tolerance(0.00001)){
//	multi::array<complex, 2> const A = {
//		{ 1.,  3.,  4.},
//		{NAN,  7.,  1.},
//		{NAN, NAN,  8.}
//	};
//	using blas::filling;
//	using blas::diagonal;
//	using blas::transposed;
//	using blas::hermitized;
//	using blas::conjugated;
//	using blas::trsv;
//	{
//		multi::array<complex, 1> b = {1., 3., 4.};
//		blas::trsv(filling::upper, diagonal::general, A, b); // B<-Solve(A.X==B), B<-A⁻¹.B, B⊤<-(A⁻¹.B)⊤, B<-B⊤.A⁻¹⊤
//		BOOST_TEST( real(b[0]) == -2.07143 );
//		BOOST_TEST( real(b[1]) ==  0.357143 );
//		BOOST_TEST( real(b[2]) ==  0.5 );
//	}
//	{
//		multi::array<complex, 1> const b = {1., 3., 4.};
//		auto b_copy = blas::trsv(filling::upper, A, b); // B<-Solve(A.X==B), B<-A⁻¹.B, B⊤<-(A⁻¹.B)⊤, B<-B⊤.A⁻¹⊤
//		BOOST_TEST( real(b[0]) == 1. );
//		BOOST_TEST( real(b_copy[0]) == -2.07143 );
//		BOOST_TEST( real(b_copy[1]) ==  0.357143 );
//		BOOST_TEST( real(b_copy[2]) ==  0.5 );
//	}
//	{
//		multi::array<complex, 1> const b = {1., 3., 4.};
//		auto b_copy = blas::trsv(filling::upper, diagonal::general, A, b); // B<-Solve(A.X==B), B<-A⁻¹.B, B⊤<-(A⁻¹.B)⊤, B<-B⊤.A⁻¹⊤
//		BOOST_TEST( real(b[0]) == 1. );
//		BOOST_TEST( real(b_copy[0]) == -2.07143 );
//		BOOST_TEST( real(b_copy[1]) ==  0.357143 );
//		BOOST_TEST( real(b_copy[2]) ==  0.5 );
//	}
//	{
//		multi::array<complex, 1> b = {3., 3., 1.};
//		trsv(filling::lower, diagonal::general, transposed(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//		BOOST_TEST( real(b[0]) ==  3. );
//		BOOST_TEST( real(b[1]) == -0.857143 );
//		BOOST_TEST( real(b[2]) == -1.26786 );
//	}
//	{
//		multi::array<complex, 1> b = {3., 3., 1.};
//	//	trsv(filling::lower, diagonal::general, hermitized(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//	//	BOOST_TEST( real(b[0]) ==  3. );
//	//	BOOST_TEST( real(b[1]) == -0.857143 );
//	//	BOOST_TEST( real(b[2]) == -1.26786 );
//	}
//	{
//		multi::array<complex, 1> b = {3., 3., 1.};
////		trsv(filling::lower, diagonal::general, hermitized(A), conjugated(b)); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
////		BOOST_TEST( real(b[0]) ==  3. );
////		BOOST_TEST( real(b[1]) == -0.857143 );
////		BOOST_TEST( real(b[2]) == -1.26786 );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsv_complex_square, *utf::tolerance(0.00001)){
//	namespace blas = multi::blas;

//	multi::array<complex, 2> const A = {
//		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
//		{NAN       ,  7. - 10.*I,  1. + 2.*I},
//		{NAN       , NAN        ,  8. + 1.*I}
//	};
//	using blas::filling;
//	using blas::diagonal;
//	using blas::transposed;
//	using blas::hermitized;
//	using blas::conjugated;
//	using blas::trsv;
//	{
//		multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
//		blas::trsv(filling::upper, diagonal::general, A, b); // B<-Solve(A.X==B), B<-A⁻¹.B, B⊤<-(A⁻¹.B)⊤, B<-B⊤.A⁻¹⊤
//		BOOST_TEST( real(b[0]) == -1.37259 );
//		BOOST_TEST( real(b[1]) ==  0.2127 );
//		BOOST_TEST( real(b[2]) ==  0.569231 );
//	}
//	{
//		multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
//		trsv(filling::lower, diagonal::general, transposed(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//		BOOST_TEST( real(b[0]) ==  1.5      ); BOOST_TEST( imag(b[0]) == 0.5        );
//		BOOST_TEST( real(b[1]) == -0.285235 ); BOOST_TEST( imag(b[1]) == -0.0503356 );
//		BOOST_TEST( real(b[2]) == -0.129272 ); BOOST_TEST( imag(b[2]) == 0.28126    );
//	}
//	{
//		multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
//		trsv(filling::upper, diagonal::general, blas::H(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//		print_1D(b);
//		BOOST_TEST( real(b[0]) == -0.661693 ); BOOST_TEST( imag(b[0]) == -1.13934   );
//		BOOST_TEST( real(b[1]) ==  0.135261 ); BOOST_TEST( imag(b[1]) == -0.0283944 ); 
//		BOOST_TEST( real(b[2]) ==  0.415385 ); BOOST_TEST( imag(b[2]) ==  0.676923  ); 
//	}
//	{
//		multi::array<complex, 1> b = {1. - 2.*I, 3. - 1.*I, 4. - 5.*I};
//		trsv(filling::upper, diagonal::general, blas::H(A), blas::conj(b)); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
////		print_1D(b);
////		BOOST_TEST( real(conjugated(b)[0]) == -0.661693 ); BOOST_TEST( imag(conjugated(b)[0]) == -1.13934   );
////		BOOST_TEST( real(conjugated(b)[1]) ==  0.135261 ); BOOST_TEST( imag(conjugated(b)[1]) == -0.0283944 ); 
////		BOOST_TEST( real(conjugated(b)[2]) ==  0.415385 ); BOOST_TEST( imag(conjugated(b)[2]) ==  0.676923  ); 
//	}
//	{
//		multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
//	//	trsv(filling::lower, diagonal::general, hermitized(A), b); // B<-Solve(A.X==B), B<-A⊤⁻¹.B, B⊤<-(A⊤⁻¹.B)⊤, B<-B⊤.A⁻¹
//	//	BOOST_TEST( real(b[0]) == -0.5      ); BOOST_TEST( imag(b[0]) ==  1.5        );
//	//	BOOST_TEST( real(b[1]) ==  0.184564 ); BOOST_TEST( imag(b[1]) == -0.620805  ); 
//	//	BOOST_TEST( real(b[2]) ==  0.691791 ); BOOST_TEST( imag(b[2]) ==  0.0227155 ); 
//	}
//}


//#if 0
//BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_1x1, *utf::tolerance(0.00001)){
//	multi::array<double, 2> const A = {
//		{10.,},
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<double, 2> B = {
//			{3.,},
//		};
//		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[0][0] == 3./10. );
//	}
//	{
//		multi::array<double, 2> B = {
//			{3.,},
//		};
//		trsm(filling::upper, diagonal::general, 2., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[0][0] == 2.*3./10. );
//	}
//	{
//		multi::array<double, 2> B = {
//			{3., 4., 5.},
//		};
//		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[0][1] == 4./10. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_0x0, *utf::tolerance(0.00001)){
//	multi::array<double, 2> const A;
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<double, 2> B;
//		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare, *utf::tolerance(0.00001)){
//	multi::array<double, 2> const A = {
//		{ 1.,  3.,  4.},
//		{ 0.,  7.,  1.},
//		{ 0.,  0.,  8.}
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<double, 2> B = {
//			{1., 3., 4., 8.},
//			{2., 7., 1., 9.},
//			{3., 4., 2., 1.},
//		};
//		multi::array<double, 2> BT = rotated(B);
//		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[1][2] == 0.107143 );

//		trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
//		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
//	}
//	{
//		multi::array<double, 2> B = {
//			{1., 3., 4., 8.},
//			{2., 7., 1., 9.},
//			{3., 4., 2., 1.},
//		};
//		multi::array<double, 2> AT = rotated(A);
//		multi::array<double, 2> BT = rotated(B);
//		trsm(filling::upper, diagonal::general, 1., rotated(AT), B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[1][2] == 0.107143 );

//		trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
//		print(rotated(BT));
//		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
//	}
//	{
//		multi::array<double, 2> B = {
//			{1.},
//			{2.},
//			{3.},
//		};
//		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
//		BOOST_TEST( B[2][0] == 0.375 );
//	}
//	{
//		multi::array<double, 2> B = {
//			{1.},
//			{2.},
//			{3.},
//		};
//		multi::array<double, 2> BT = rotated(B);
//		trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
//		BOOST_TEST( rotated(BT)[2][0] == 0.375 );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare_default_diagonal_gemm_check, *utf::tolerance(0.00001)){
//	multi::array<double, 2> const A = {
//		{ 1.,  3.,  4.},
//		{ 0.,  7.,  1.},
//		{ 0.,  0.,  8.}
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<double, 2> const B = {
//			{1.},// 3., 4.},
//			{2.},// 7., 1.},
//			{3.},// 4., 2.},
//		};
//		using multi::blas::gemm;
//		{
//			auto S = trsm(filling::upper, diagonal::general, 1., A, B);
//			BOOST_REQUIRE( S[2][0] == 0.375 );
//			auto Bck=gemm(1., A, S);
//			BOOST_REQUIRE( Bck[2][0] == 3. );
//			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//		{
//			multi::array<double, 2> const BT = rotated(B);
//			auto Bck=gemm(1., A, trsm(filling::upper, diagonal::general, 1., A, rotated(BT)));
//			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//		{
//			auto const AT = rotated(A);
//			auto Bck=gemm(1., rotated(AT), trsm(filling::upper, diagonal::general, 1., rotated(AT), B));
//			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//		{
//			auto const AT =* rotated(A);
//			auto const BT =* rotated(B);
//			auto const Bck=gemm(1., A, trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT)));
//			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_REQUIRE_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//		{
//			auto const AT =* rotated(A);
//			auto const BT =* rotated(B);
//			using multi::blas::trsm;
//		//	auto const Bck=gemm(A, trsm(rotated(AT), rotated(BT)));
//		//	for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//		{
//			using multi::blas::trsm;
//		//	auto const Bck=gemm(A, trsm(A, B));
//		//	for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_1x1_check, *utf::tolerance(0.00001)){
//	multi::array<double, 2> const A = {
//		{ 4.},
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<double, 2> const B = {
//			{5.},
//		};
//		{
//			auto S = trsm(filling::upper, diagonal::general, 3., A, B);
//			BOOST_REQUIRE( S[0][0] == 3.*5./4. );
//		}
//		{
//			auto S = trsm(filling::upper, 1., A, B);
//			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
//		}
//		{
//			auto S = trsm(filling::upper, A, B);
//			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_1x1_check, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I = complex{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 4. + 2.*I},
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<complex, 2> const B = {
//			{5. + 1.*I},
//		};
//		using multi::blas::gemm;
//		{
//			auto S = trsm(filling::upper, diagonal::general, 3.+5.*I, A, B);
//			BOOST_TEST( real(S[0][0]) == real((3.+5.*I)*B[0][0]/A[0][0]) );
//			BOOST_TEST( imag(S[0][0]) == imag((3.+5.*I)*B[0][0]/A[0][0]) );
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_one_check, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 1. + 4.*I,  3.,  4.- 10.*I},
//		{ 0.,  7.- 3.*I,  1.},
//		{ 0.,  0.,  8.- 2.*I}
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<complex, 2> const B = {
//			{1. + 1.*I},
//			{2. + 1.*I},
//			{3. + 1.*I},
//		};
//		using multi::blas::gemm;
//		{
//			auto S = trsm(filling::upper, diagonal::general, 1., A, B);
//			BOOST_TEST( real(S[2][0]) == 0.323529 );
//		}
//		{
//			auto const BT = +rotated(B);
//			auto S = trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
//			BOOST_TEST( real(S[2][0]) == 0.323529 );
//		}
//		{
//			auto const AT = +rotated(A);
//			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), B);
//			BOOST_TEST( real(S[2][0]) == 0.323529 );
//		}
//		{
//			auto const AT = +rotated(A);
//			auto const BT = +rotated(B);
//			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
//			BOOST_TEST( real(S[2][0]) == 0.323529 );
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_gemm_check, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 1. + 4.*I,  3.,  4.- 10.*I},
//		{ 0.,  7.- 3.*I,  1.},
//		{ 0.,  0.,  8.- 2.*I}
//	};
//	using multi::blas::side;
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<complex, 2> const B = {
//			{1. + 1.*I, 5. + 3.*I},
//			{2. + 1.*I, 9. + 3.*I},
//			{3. + 1.*I, 1. - 1.*I},
//		};
//		using multi::blas::gemm;
//		{
//			auto S = trsm(filling::upper, diagonal::general, 1., A, B); // S = Ainv.B
//			BOOST_TEST( real(S[2][1]) == 0.147059  );
//		}
//		{
//			auto const BT = +rotated(B);
//			auto S = trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
//			BOOST_TEST( real(S[2][1]) == 0.147059  );
//		}
//		{
//			auto const AT = +rotated(A);
//			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), B);
//			BOOST_TEST( real(S[2][1]) == 0.147059  );
//		}
//		{
//			auto const AT = +rotated(A);
//			auto const BT = +rotated(B);
//			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
//			BOOST_TEST( real(S[2][1]) == 0.147059  );
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 1. + 4.*I,  3.,  4.- 10.*I},
//		{ 0.,  7.- 3.*I,  1.},
//		{ 0.,  0.,  8.- 2.*I}
//	};
//	using multi::blas::filling;
//	using multi::blas::diagonal;
//	{
//		multi::array<complex, 2> const B = {
//			{1. + 1.*I, 5. + 3.*I},
//			{2. + 1.*I, 9. + 3.*I},
//			{3. + 1.*I, 1. - 1.*I},
//		};
//		using multi::blas::hermitized;
//		{
//			auto S = trsm(filling::lower, diagonal::general, 1., hermitized(A), B); // S = A⁻¹†.B, S† = B†.A⁻¹
//			BOOST_TEST( real(S[2][1]) == 1.71608  );
//		}
//		{
//			multi::array<complex, 2> const B = {
//				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
//				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
//			};
//			auto S =* trsm(filling::upper, 1., A, hermitized(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
//			BOOST_TEST( imag(S[2][1]) == +0.147059 );
//			BOOST_TEST( imag(B[1][2]) == -0.147059 );
//		}
//		{
//			multi::array<complex, 2> const B = {
//				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
//				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
//			};
//			auto S =* trsm(filling::upper, 2., A, hermitized(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
//			BOOST_TEST( imag(S[2][1]) == +0.147059*2. );
//			BOOST_TEST( imag(B[1][2]) == -0.147059*2. );
//		}
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 1. + 4.*I,  3.,  4.- 10.*I},
//		{ 0.,  7.- 3.*I,  1.},
//		{ 0.,  0.,  8.- 2.*I}
//	};
//	multi::array<complex, 2> B = {
//		{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
//		{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
//	};
//	using multi::blas::trsm;
//	using multi::blas::filling;
//	using multi::blas::hermitized;
//	trsm(filling::upper, A, hermitized(B)); // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
//	BOOST_TEST( imag(B[1][2]) == -0.147059 );
//}
//#endif
//#endif


//#endif
#endif
