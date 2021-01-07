#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework \
`#-Wl,-rpath,/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -L/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core` \
&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_GEMV_HPP
#define MULTI_ADAPTORS_BLAS_GEMV_HPP

#include "../blas/core.hpp"

#include "./../../detail/../utility.hpp"
#include "./../../detail/../array_ref.hpp"

//#include<cblas64.h>

#include "../blas/operations.hpp"

namespace boost{
namespace multi{namespace blas{

using multi::blas::core::gemv;

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0>
auto gemv_base_aux(A&& a){return base(a);}

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto gemv_base_aux(A&& a){return underlying(base(a));}

template<class A, class X, class Y>
auto gemv(typename A::element alpha, A const& a, X const& x, typename A::element beta, Y&& y)
->decltype(gemv('N', x.size(), y.size(), alpha, gemv_base_aux(a), a.rotated().stride(), x.base(), x.stride(), beta, y.base(), y.stride()), std::forward<Y>(y)){
	assert(a.rotated().size() == x.size());
	assert(a.size() == y.size());
	
	auto base_A = gemv_base_aux(a);
	
	assert(stride(a)==1 or stride(~a)==1);
	
	assert(y.size() == a.size());
	     if(stride( a)==1 and not is_conjugated<A>{})        gemv('N'                          , size(y), size(x),  alpha, base_A, stride(~a), base(x), stride(x),  beta, base(y), stride(y));
	else if(stride(~a)==1 and not is_conjugated<A>{})        gemv('T'                          , size(x), size(y),  alpha, base_A, stride( a), base(x), stride(x),  beta, base(y), stride(y));
	else if(stride(~a)==1 and     is_conjugated<A>{})        gemv('C'                          , size(x), size(y),  alpha, base_A, stride( a), base(x), stride(x),  beta, base(y), stride(y));
	else if(stride( a)==1 and     is_conjugated<A>{}){
		assert(0&&"case not supported by blas");
	//	cblas_zgemv(CblasRowMajor, CblasConjTrans, size(x), size(y), &alpha, base_A, stride(~a), base(x), stride(x), &beta, base(y), stride(y));
	}
	
	return std::forward<Y>(y);
}

template<class A2D, class X1D, class Y1D>
auto gemv(A2D const& A, X1D const& x, Y1D&& y)
->decltype(gemv(1., A, x, 0., std::forward<Y1D>(y))){
	return gemv(1., A, x, 0., std::forward<Y1D>(y));}

template<class Alloc, class A2D, class X1D, typename T = typename A2D::element>
NODISCARD("")
auto gemv(A2D const& A, X1D const& x, Alloc const& alloc = {})->std::decay_t<
decltype(gemv(1., A, x, 0., multi::array<T, 1, Alloc>(A.size(), alloc)))>{
	return gemv(1., A, x, 0., multi::array<T, 1, Alloc>(A.size(), alloc));}

template<class A2D, class X1D>
NODISCARD("")
auto gemv(A2D const& A, X1D const& x){
	return gemv(A, x, get_allocator(x));
}

namespace operators{
	template<class A2D, class X1D> auto operator%(A2D const& A, X1D const& x) DECLRETURN(gemv(A, x))
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_GEMV

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi blas gemv"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include "../blas/dot.hpp"
#include "../blas/axpy.hpp"
#include "../blas/nrm2.hpp"
#include "../blas/gemm.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>
#include<random>

using std::cout;
namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_blas_gemv_real){
	namespace blas = multi::blas;

	using std::abs;
	multi::array<double, 2> const M = {
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<double, 1> const X = {1.1,2.1,3.1, 4.1};
	{
		multi::array<double, 1> Y = {4.,5.,6.};
		double const a = 1.1, b = 1.2;
		blas::gemv(a, M, X, b, Y); // y = a*M*x + b*y

		multi::array<double, 1> const Y3 = {214.02, 106.43, 188.37};
		BOOST_REQUIRE( abs(Y[1] - Y3[1]) < 2e-14 );
	}
	{
		auto Y = blas::gemv(M, X);
		BOOST_REQUIRE((
			Y == decltype(Y){
				blas::dot(M[0], X),
				blas::dot(M[1], X),
				blas::dot(M[2], X)
			}
		));
		BOOST_REQUIRE(std::equal(begin(Y), end(Y), begin(M), [&X](auto&& y, auto&& m){return y==blas::dot(m, X);}));
	}
	{
		multi::array<double, 1> const a = {1., 2., 3.};
		multi::array<double, 1> const b = {4., 5., 6.};
		BOOST_REQUIRE(
			blas::gemv(multi::array<double, 2>({a}), b)[0] == blas::dot(a, b)
		);
	}
	{
		multi::array<double, 2> const MT = ~M;
		using boost::test_tools::tolerance;
	//	using blas::gemv; BOOST_TEST_REQUIRE( nrm2(blas::axpy(-1., gemv(~+~M, X), gemv(M, X)))() == 0.,  tt::tolerance(1e-13) );
		using namespace blas;
		using namespace blas::operators;
	//	blas::nrm2(X);
		BOOST_TEST_REQUIRE( (((~+~M)%X - M%X)^2) == 0., tolerance(1e-13) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_real_complex){
	namespace blas = multi::blas;
	using complex = std::complex<double>; //#define I *std::complex<double>(0, 1)
	using std::abs;
	multi::array<complex, 2> const M = {
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<complex, 1> const X = {1.1, 2.1, 3.1, 4.1};
	{
		multi::array<complex, 1> Y = {4., 5., 6.};
		double const a = 1.1, b = 1.2;
		blas::gemv(a, M, X, b, Y); // y = a*M*x + b*y
		
		multi::array<complex, 1> const Y3 = {214.02, 106.43, 188.37};
		
		using namespace blas::operators;
		BOOST_TEST_REQUIRE( ((Y - Y3)^2)  == 0. , boost::test_tools::tolerance(1e-13) );
	
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_complex){
	
	namespace blas = multi::blas;
	using complex = std::complex<double>; std::complex<double> const I{0, 1};
	
	using std::abs;
	multi::array<complex, 2> const M = {{2. + 3.*I, 2. + 1.*I, 1. + 2.*I}, {4. + 2.*I, 2. + 4.*I, 3. + 1.*I}, 
 {7. + 1.*I, 1. + 5.*I, 0. + 3.*I}};
	multi::array<complex, 1> const X = {1. + 2.*I, 2. + 1.*I, 9. + 2.*I};
	using namespace blas::operators;
	BOOST_REQUIRE(( blas::gemv(   M, X) == multi::array<complex, 1>{4. + 31.*I, 25. + 35.*I, -4. + 53.*I} ));
	
	auto MT = +~M;
	BOOST_REQUIRE(( blas::gemv(~MT, X) == multi::array<complex, 1>{4. + 31.*I, 25. + 35.*I, -4. + 53.*I} ));
	
	auto MH = +*~M;
	BOOST_REQUIRE( blas::gemv(~M, X) == (multi::array<complex, 1>{63. + 38.*I, -1. + 62.*I, -4. + 36.*I}) );
	BOOST_REQUIRE( blas::gemv(~M, X) == blas::gemv(MT, X) );// == multi::array<complex, 1>{4. + 31.*I, 25. + 35.*I, -4. + 53.*I} ));
	
	BOOST_REQUIRE( blas::gemv(*M, X) == (multi::array<complex, 1>{26. - 15.*I, 45. - 3.*I, 22. - 23.*I}) );
	
//	BOOST_REQUIRE( blas::gemv(~*M, X) == (multi::array<complex, 1>{83. + 6.*I, 31. - 46.*I, 18. - 26.*I}) ); // not supported by blas

}

BOOST_AUTO_TEST_CASE(multi_blas_gemv_temporary){

	using complex = std::complex<double>;
	
	multi::array<complex, 2> const A = {
		{1., 0., 0.}, 
		{0., 1., 0.},
		{0., 0., 1.}
	};
	
	auto const B = []{
		multi::array<complex, 2> B({3, 3});
		auto rand = [d=std::normal_distribution<>{}, g=std::mt19937{}]()mutable{return complex{d(g), d(g)};};
		std::generate(B.elements().begin(), B.elements().end(), rand);
		return B;
	}();
	
	namespace blas = multi::blas;
	
	using namespace blas::operators;
	BOOST_TEST( ( (A%(~B)[0] - ( ~(A*B)  )[0])^2) == 0. );
	BOOST_TEST( ( (A%(~B)[0] - ((~B)*(~A))[0])^2) == 0. );

}

#if 0
	{
		auto Y = blas::gemv(M, X);
		BOOST_REQUIRE((
			Y == decltype(Y){
				blas::dot(M[0], X),
				blas::dot(M[1], X),
				blas::dot(M[2], X)
			}
		));
		BOOST_REQUIRE(std::equal(begin(Y), end(Y), begin(M), [&X](auto&& y, auto&& m){return y==blas::dot(m, X);}));
	}
	{
		multi::array<double, 1> const a = {1., 2., 3.};
		multi::array<double, 1> const b = {4., 5., 6.};
		BOOST_REQUIRE(
			blas::gemv(multi::array<double, 2>({a}), b)[0] == blas::dot(a, b)
		);
	}
	{
		multi::array<double, 2> const MT = ~M;
		using boost::test_tools::tolerance;
	//	using blas::gemv; BOOST_TEST_REQUIRE( nrm2(blas::axpy(-1., gemv(~+~M, X), gemv(M, X)))() == 0.,  tt::tolerance(1e-13) );
		using namespace blas;
		using namespace blas::operators;
	//	blas::nrm2(X);
		BOOST_TEST_REQUIRE( (((~+~M|X)-(M|X))^2) == 0., tolerance(1e-13) );
	}
#endif


#if 0
	{
		double const M[3][4] = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		double const X[4] = {1.1,2.1,3.1, 4.1};
		double Y[3] = {4.,5.,6.};
		double const a = 1.1;
		double const b = 1.2;
		gemv('T', a, M, X, b, Y); // y = a*M*x + b*y
		double Y3[3] = {214.02, 106.43, 188.37};
		assert( abs(Y[1] - Y3[1]) < 2e-14 );
	}

	{
		multi::array<double, 2> const M = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		multi::array<double, 1> const X = {1.1,2.1,3.1};
		multi::array<double, 1> Y = {4.,5.,6., 7.};
		double a = 1.8, b = 1.6;
		gemv('N', a, M, X, b, Y); // y = a*(M^T)*x + b*y, y^T = a*(x^T)*M + b*y^T
		multi::array<double, 1> const Y3 = {117.46, 182.6, 315.24, 61.06}; // =1.8 Transpose[{{9., 24., 30., 9.}, {4., 10., 12., 7.}, {14., 16., 36., 1.}}].{1.1, 2.1, 3.1} + 1.6 {4., 5., 6., 7.}
		assert( abs(Y[2] - Y3[2]) < 1e-13 );
	}
	{
		multi::array<double, 2> const M = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		multi::array<double, 1> const X = {1.1,2.1,3.1, 4.1};
		multi::array<double, 1> Y = {4.,5.,6.};
		double a = 1.1, b = 1.2;
		gemv(a, M, X, b, Y); // y = a*M*x + b*y
		multi::array<double, 1> const Y3 = {214.02, 106.43, 188.37}; // = 1.1 {{9., 24., 30., 9.}, {4., 10., 12., 7.}, {14., 16., 36., 1.}}.{1.1, 2.1, 3.1, 4.1} + 1.2 {4., 5., 6.}
		assert( std::abs(Y[1] - Y3[1]) < 2e-14 );
	}
	{
		double const M[3][4] = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		double const X[4] = {1.1,2.1,3.1, 4.1};
		double Y[3] = {4.,5.,6.};
		double a = 1.1, b = 1.2;
		gemv(a, M, X, b, Y); // y = a*M*x + b*y
		double const Y3[3] = {214.02, 106.43, 188.37};
		assert( std::abs(Y[1] - Y3[1]) < 2e-14 );
	}
	{
		multi::array<double, 2> const M = {
			{ 9., 4., 14.},
			{24., 10., 16.},
			{30., 12., 36.},
			{9., 7., 1.}
		}; assert( M[0][2] == 14. );
		multi::array<double, 1> const X = {1.1,2.1,3.1, 4.1};
		multi::array<double, 1> Y = {4.,5.,6.};
		double a = 1.1, b = 1.2;
		gemv(a, rotated(M), X, b, Y); // y = a*M*x + b*y

		multi::array<double, 1> const Y3 = {214.02, 106.43, 188.37}; // = 1.1 {{9., 24., 30., 9.}, {4., 10., 12., 7.}, {14., 16., 36., 1.}}.{1.1, 2.1, 3.1, 4.1} + 1.2 {4., 5., 6.}
		assert( abs(Y[1] - Y3[1]) < 2e-14 );
	}
#endif

#endif
#endif

