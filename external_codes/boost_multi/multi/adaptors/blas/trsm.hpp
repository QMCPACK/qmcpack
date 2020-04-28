#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&`#nvcc -x cu --expt-relaxed-constexpr`$CXX -D_TEST_MULTI_ADAPTORS_BLAS_TRSM $0.cpp -o $0x -lboost_unit_test_framework \
`pkg-config --cflags --libs blas` \
`#-Wl,-rpath,/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -L/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core` \
-lboost_timer &&$0x&& rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_TRSM_HPP
#define MULTI_ADAPTORS_BLAS_TRSM_HPP

#include "../blas/core.hpp"

#include "../blas/operations.hpp" // uplo
#include "../blas/filling.hpp"
#include "../blas/side.hpp"

#include "../../config/NODISCARD.hpp"

namespace boost{
namespace multi{namespace blas{

enum class DIAG : char{U='U', N='N'};

enum class diagonal : char{
	unit = static_cast<char>(DIAG::U),
	non_unit = static_cast<char>(DIAG::N), general = non_unit
};

template<class A> auto trsm_base_aux(A&& a, std::false_type){return base(a);}
template<class A> auto trsm_base_aux(A&& a, std::true_type){return underlying(base(a));}

using core::trsm;

template<class A2D, class B2D>
auto trsm_move(filling a_nonz, diagonal a_diag, typename A2D::element_type alpha, A2D const& a, B2D&& b, side s = side::left)
->decltype(
	trsm('R', static_cast<char>(+a_nonz), 'N', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, trsm_base_aux(a, is_hermitized<A2D>{}), stride(a)         , trsm_base_aux(b, is_hermitized<B2D>{}), stride(b))
	, std::forward<B2D>(b)
)
{
	if(s==side::left) assert( size(rotated(a)) == size(b) );
	else              assert( size(rotated(b)) == size(a) );

	if(size(b)==0) return std::forward<B2D>(b);

	auto base_a = trsm_base_aux(a, is_hermitized<A2D>{}); (void)base_a;
	auto base_b = trsm_base_aux(b, is_hermitized<B2D>{}); (void)base_b;

	using core::trsm;
	if(size(rotated(b))==1){
		if(stride(rotated(a))==1) trsm(static_cast<char>(s), static_cast<char>(+a_nonz), 'N', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(a)         , base_b, stride(b));
		else if(stride(a)==1){		
			if(!is_conjugated(a)) trsm(static_cast<char>(s), static_cast<char>(-a_nonz), 'T', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(rotated(a)), base_b, stride(b));
			else                  trsm(static_cast<char>(s), static_cast<char>(-a_nonz), 'C', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(rotated(a)), base_b, stride(b));
		}
	}else{
		      if(stride(rotated(a))==1 and stride(rotated(b))==1){
			assert(!is_conjugated(a));
			trsm(static_cast<char>(s), static_cast<char>(+a_nonz), 'N', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(a)         , base_b, stride(b));
		}else if(stride(a)==1 and stride(rotated(b))==1){
			if(!is_conjugated(a)) trsm(static_cast<char>(s), static_cast<char>(-a_nonz), 'T', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(rotated(a)), base_b, stride(b));
			else                  trsm(static_cast<char>(s), static_cast<char>(-a_nonz), 'C', static_cast<char>(a_diag), size(rotated(b)), size(b), alpha, base_a, stride(rotated(a)), base_b, stride(b));
		}
		else if(stride(a)==1 and stride(b)==1){
			assert(!is_conjugated(a));
			                      trsm(static_cast<char>(swap(s)), static_cast<char>(-a_nonz), 'N', static_cast<char>(a_diag), size(rotated(a)), size(rotated(b)), alpha, base_a, stride(rotated(a)), base_b, stride(rotated(b)));
		}else if(stride(rotated(a))==1 and stride(b)==1){
			if(!is_conjugated(a)) trsm(static_cast<char>(swap(s)), static_cast<char>(+a_nonz), 'T', static_cast<char>(a_diag), size(a), size(rotated(b)), alpha, base_a, stride(a), base_b, stride(rotated(b)));
			else                  trsm(static_cast<char>(swap(s)), static_cast<char>(+a_nonz), 'C', static_cast<char>(a_diag), size(a), size(rotated(b)), alpha, base_a, stride(a), base_b, stride(rotated(b)));
		}
		else assert(0); // case is not part of BLAS
	}
	return std::forward<B2D>(b);
}

template<typename AA, class A2D, class B2D>
auto trsm(filling a_nonz, diagonal a_diag, AA alpha, A2D const& a, B2D&& b)
->decltype(trsm_move(a_nonz, a_diag, alpha, a, std::forward<B2D>(b)), std::declval<B2D&&>())
{
	if(!is_conjugated(b)) trsm_move( a_nonz, a_diag, alpha,            a, std::forward<B2D>(b));
	else                  trsm_move(-a_nonz, a_diag, alpha, hermitized(a), rotated(b), side::right);
	return std::forward<B2D>(b);
}

template<typename AA, class A2D, class B2D>
NODISCARD("because last argument is const")
auto trsm(filling a_nonz, diagonal a_diag, AA alpha, A2D const& a, B2D const& b)
->std::decay_t<decltype(trsm_move(a_nonz, a_diag, alpha, a, decay(b)))>{
	return trsm_move(a_nonz, a_diag, alpha, a, decay(b));}

template<class AA, class A2D, class B2D>
auto trsm(filling a_nonz, AA alpha, A2D const& a, B2D&& b)
->decltype(trsm(a_nonz, diagonal::general, alpha, a, std::forward<B2D>(b))){
	return trsm(a_nonz, diagonal::general, alpha, a, std::forward<B2D>(b));}

template<class AA, class A2D, class B2D>
NODISCARD("because last argument is const")
auto trsm(filling a_nonz, AA alpha, A2D const& a, B2D const& b)
->decltype(trsm(a_nonz, diagonal::general, alpha, a, std::forward<B2D>(b))){
	return trsm(a_nonz, diagonal::general, alpha, a, std::forward<B2D>(b));}

template<class A2D, class B2D, class T = typename A2D::element_type>
auto trsm(filling a_nonz, A2D const& a, B2D&& b)
->decltype(trsm(a_nonz, T{1.}, a, std::forward<B2D>(b))){
	return trsm(a_nonz, T{1.}, a, std::forward<B2D>(b));}

template<class A2D, class B2D, class T = typename A2D::element_type>
NODISCARD("because last argument is const")
auto trsm(filling a_nonz, A2D const& a, B2D const& b)
->decltype(trsm(a_nonz, T{1.}, a, std::forward<B2D>(b))){
	return trsm(a_nonz, T{1.}, a, std::forward<B2D>(b));}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_TRSM

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi.BLAS trsm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include<boost/test/floating_point_comparison.hpp>

#include "../blas/gemm.hpp"

#include "../../array.hpp"

namespace multi = boost::multi;

#include<iostream>
#include<vector>

template<class M> decltype(auto) print(M const& C){
	using boost::multi::size; using std::cout;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j)
			cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<std::endl;
}

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_1x1, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A = {
		{10.,},
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> B = {
			{3.,},
		};
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[0][0] == 3./10. );
	}
	{
		multi::array<double, 2> B = {
			{3.,},
		};
		trsm(filling::upper, diagonal::general, 2., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[0][0] == 2.*3./10. );
	}
	{
		multi::array<double, 2> B = {
			{3., 4., 5.},
		};
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[0][1] == 4./10. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_0x0, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A;
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> B;
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_square, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> B = {
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[1][2] == 0.107143 );
	}
	{
		multi::array<double, 2> AT = rotated(A);
		multi::array<double, 2> B = {
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		trsm(filling::upper, diagonal::general, 1., rotated(AT), B);
		BOOST_TEST( B[1][2] == 0.107143 );
	}
	{
		multi::array<double, 2> AT = rotated(A);
		multi::array<double, 2> B = {
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		multi::array<double, 2> BT = rotated(B);
		trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
	}
	{
		multi::array<double, 2> AT = rotated(A);
		multi::array<double, 2> B = {
			{1., 3., 4.},
			{2., 7., 1.},
			{3., 4., 2.}
		};
		multi::array<double, 2> BT = rotated(B);
		trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{ 0.,  7.,  1.},
		{ 0.,  0.,  8.}
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> B = {
			{1., 3., 4., 8.},
			{2., 7., 1., 9.},
			{3., 4., 2., 1.},
		};
		multi::array<double, 2> BT = rotated(B);
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[1][2] == 0.107143 );

		trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
	}
	{
		multi::array<double, 2> B = {
			{1., 3., 4., 8.},
			{2., 7., 1., 9.},
			{3., 4., 2., 1.},
		};
		multi::array<double, 2> AT = rotated(A);
		multi::array<double, 2> BT = rotated(B);
		trsm(filling::upper, diagonal::general, 1., rotated(AT), B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[1][2] == 0.107143 );

		trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
		print(rotated(BT));
		BOOST_TEST( rotated(BT)[1][2] == 0.107143 );
	}
	{
		multi::array<double, 2> B = {
			{1.},
			{2.},
			{3.},
		};
		trsm(filling::upper, diagonal::general, 1., A, B); // B=Solve(A.X=alpha*B, X) B=A⁻¹B, B⊤=B⊤.(A⊤)⁻¹, A upper triangular (implicit zeros below)
		BOOST_TEST( B[2][0] == 0.375 );
	}
	{
		multi::array<double, 2> B = {
			{1.},
			{2.},
			{3.},
		};
		multi::array<double, 2> BT = rotated(B);
		trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
		BOOST_TEST( rotated(BT)[2][0] == 0.375 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_nonsquare_default_diagonal_gemm_check, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{ 0.,  7.,  1.},
		{ 0.,  0.,  8.}
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> const B = {
			{1.},
			{2.},
			{3.}
		};
		using multi::blas::gemm;
		{
			auto S = trsm(filling::upper, diagonal::general, 1., A, B);
			BOOST_REQUIRE( S[2][0] == 0.375 );
			auto Bck=gemm(1., A, S);
			BOOST_REQUIRE( Bck[2][0] == 3. );
			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
		{
			multi::array<double, 2> const BT = rotated(B);
			auto Bck=gemm(1., A, trsm(filling::upper, diagonal::general, 1., A, rotated(BT)));
			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
		{
			auto const AT = rotated(A);
			auto Bck=gemm(1., rotated(AT), trsm(filling::upper, diagonal::general, 1., rotated(AT), B));
			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
		{
			auto const AT =* rotated(A);
			auto const BT =* rotated(B);
			auto const Bck=gemm(1., A, trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT)));
			for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_REQUIRE_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
		{
			auto const AT =* rotated(A);
			auto const BT =* rotated(B);
			using multi::blas::trsm;
		//	auto const Bck=gemm(A, trsm(rotated(AT), rotated(BT)));
		//	for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
		{
			using multi::blas::trsm;
		//	auto const Bck=gemm(A, trsm(A, B));
		//	for(int i{};i<3;++i)for(int j{};j<size(rotated(B));++j) BOOST_CHECK_SMALL(Bck[i][j]-B[i][j], 0.00001);
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_real_1x1_check, *utf::tolerance(0.00001)){
	multi::array<double, 2> const A = {
		{ 4.},
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<double, 2> const B = {
			{5.},
		};
		{
			auto S = trsm(filling::upper, diagonal::general, 3., A, B);
			BOOST_REQUIRE( S[0][0] == 3.*5./4. );
		}
		{
			auto S = trsm(filling::upper, 1., A, B);
			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
		}
		{
			auto S = trsm(filling::upper, A, B);
			BOOST_REQUIRE( S[0][0] == 1.*5./4. );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_1x1_check, *utf::tolerance(0.00001)){
	using complex = std::complex<double>; complex const I = complex{0, 1};
	multi::array<complex, 2> const A = {
		{ 4. + 2.*I},
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<complex, 2> const B = {
			{5. + 1.*I},
		};
		using multi::blas::gemm;
		{
			auto S = trsm(filling::upper, diagonal::general, 3.+5.*I, A, B);
			BOOST_TEST( real(S[0][0]) == real((3.+5.*I)*B[0][0]/A[0][0]) );
			BOOST_TEST( imag(S[0][0]) == imag((3.+5.*I)*B[0][0]/A[0][0]) );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_one_check, *utf::tolerance(0.00001)){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.,  4.- 10.*I},
		{ 0.,  7.- 3.*I,  1.},
		{ 0.,  0.,  8.- 2.*I}
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<complex, 2> const B = {
			{1. + 1.*I},
			{2. + 1.*I},
			{3. + 1.*I},
		};
		using multi::blas::gemm;
		{
			auto S = trsm(filling::upper, diagonal::general, 1., A, B);
			BOOST_TEST( real(S[2][0]) == 0.323529 );
		}
		{
			auto const BT = +rotated(B);
			auto S = trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
			BOOST_TEST( real(S[2][0]) == 0.323529 );
		}
		{
			auto const AT = +rotated(A);
			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), B);
			BOOST_TEST( real(S[2][0]) == 0.323529 );
		}
		{
			auto const AT = +rotated(A);
			auto const BT = +rotated(B);
			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
			BOOST_TEST( real(S[2][0]) == 0.323529 );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_gemm_check, *utf::tolerance(0.00001)){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.,  4.- 10.*I},
		{ 0.,  7.- 3.*I,  1.},
		{ 0.,  0.,  8.- 2.*I}
	};
	using multi::blas::side;
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<complex, 2> const B = {
			{1. + 1.*I, 5. + 3.*I},
			{2. + 1.*I, 9. + 3.*I},
			{3. + 1.*I, 1. - 1.*I},
		};
		using multi::blas::gemm;
		{
			auto S = trsm(filling::upper, diagonal::general, 1., A, B); // S = Ainv.B
			BOOST_TEST( real(S[2][1]) == 0.147059  );
		}
		{
			auto const BT = +rotated(B);
			auto S = trsm(filling::upper, diagonal::general, 1., A, rotated(BT));
			BOOST_TEST( real(S[2][1]) == 0.147059  );
		}
		{
			auto const AT = +rotated(A);
			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), B);
			BOOST_TEST( real(S[2][1]) == 0.147059  );
		}
		{
			auto const AT = +rotated(A);
			auto const BT = +rotated(B);
			auto S = trsm(filling::upper, diagonal::general, 1., rotated(AT), rotated(BT));
			BOOST_TEST( real(S[2][1]) == 0.147059  );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check, *utf::tolerance(0.00001)){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.,  4.- 10.*I},
		{ 0.,  7.- 3.*I,  1.},
		{ 0.,  0.,  8.- 2.*I}
	};
	using multi::blas::filling;
	using multi::blas::diagonal;
	{
		multi::array<complex, 2> const B = {
			{1. + 1.*I, 5. + 3.*I},
			{2. + 1.*I, 9. + 3.*I},
			{3. + 1.*I, 1. - 1.*I},
		};
		using multi::blas::hermitized;
		{
			auto S = trsm(filling::lower, diagonal::general, 1., hermitized(A), B); // S = A⁻¹†.B, S† = B†.A⁻¹
			BOOST_TEST( real(S[2][1]) == 1.71608  );
		}
		{
			multi::array<complex, 2> const B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =* trsm(filling::upper, 1., A, hermitized(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_TEST( imag(S[2][1]) == +0.147059 );
			BOOST_TEST( imag(B[1][2]) == -0.147059 );
		}
		{
			multi::array<complex, 2> const B = {
				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
			};
			auto S =* trsm(filling::upper, 2., A, hermitized(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
			BOOST_TEST( imag(S[2][1]) == +0.147059*2. );
			BOOST_TEST( imag(B[1][2]) == -0.147059*2. );
		}
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const, *utf::tolerance(0.00001)){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{ 1. + 4.*I,  3.,  4.- 10.*I},
		{ 0.,  7.- 3.*I,  1.},
		{ 0.,  0.,  8.- 2.*I}
	};
	multi::array<complex, 2> B = {
		{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
		{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
	};
	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;
	trsm(filling::upper, A, hermitized(B)); // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_TEST( imag(B[1][2]) == -0.147059 );
}


#endif
#endif

