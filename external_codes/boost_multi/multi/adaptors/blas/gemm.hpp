#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS -O3 $0 -o $0x -lboost_unit_test_framework -lboost_timer \
`pkg-config --libs blas` \
`#-Wl,-rpath,/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -L/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core` \
&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_GEMM_HPP
#define MULTI_ADAPTORS_BLAS_GEMM_HPP

#include "../blas/core.hpp"

#include "../blas/numeric.hpp"
#include "../blas/operations.hpp"

#include "../blas/gemv.hpp"

#include "../../config/NODISCARD.hpp"
#include "../../config/MARK.hpp"

namespace boost{
namespace multi{
namespace blas{

using multi::blas::core::gemm;
using multi::blas::core::gemv;

template<class M> bool is_c_ordering(M const& m){
	return stride(rotated(m))==1 and size(m)!=1;
}

double const& conj(double const& d){return d;}
float const& conj(float const& f){return f;}

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0>
auto gemm_base_aux(A&& a){return base(a);}

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto gemm_base_aux(A&& a){return underlying(base(a));}

template<class Context, class A, class B, class C>
C&& gemm(Context& ctx, typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c)
//->decltype(ctx.gemm('N', 'T', size(~c), size(a), size(b), &alpha, gemm_base_aux(b), stride( b), gemm_base_aux(a), stride(~a), &beta, gemm_base_aux(c), size(b)) , std::forward<C>(c))
{

	MULTI_MARK_SCOPE("multi::blas::gemm with context");

	if(c.is_empty()){
		assert(a.is_empty() and b.is_empty());
		return std::forward<C>(c);
	}

	assert( size(~a) == size( b) );
	assert( size( a) == size( c) );
	assert( size(~b) == size(~c) );
	
	auto base_a = gemm_base_aux(a);
	auto base_b = gemm_base_aux(b);
	auto base_c = gemm_base_aux(c);

	assert( stride(a)==1 or stride(~a)==1 );
	assert( stride(b)==1 or stride(~b)==1 );
	assert( stride(c)==1 or stride(~c)==1 );
	
	     if(stride(c)==1 and stride(~c)!=1) blas::gemm(ctx, alpha, ~b, ~a, beta, ~c);
	else if(is_conjugated<C>{}) blas::gemm(ctx, conj(alpha), conj(a), conj(b), conj(beta), conj(c));
	else{
		;;;;; if(stride(~a)==1 and stride(~b)==1 and not is_conjugated<A>{} and not is_conjugated<B>{}){
			if(size(a)==1) ctx.gemm('N', 'N', size(~c), size(a), size(b), &alpha, base_b, stride( b), base_a, size(b)   , &beta, base_c, size(b)  );
			else           ctx.gemm('N', 'N', size(~c), size(a), size(b), &alpha, base_b, stride( b), base_a, stride( a), &beta, base_c, stride(c));
		}else if(stride( a)==1 and stride(~b)==1 and     is_conjugated<A>{} and not is_conjugated<B>{}) ctx.gemm('N', 'C', size(~c), size(a), size(b), &alpha, base_b, stride( b), base_a, stride(~a), &beta, base_c, stride(c));
		else if(stride( a)==1 and stride(~b)==1 and not is_conjugated<A>{} and not is_conjugated<B>{}){
			if(size(a)==1) ctx.gemm('N', 'T', size(~c), size(a), size(b), &alpha, base_b, stride( b), base_a, stride(~a), &beta, base_c, size(b));
			else           ctx.gemm('N', 'T', size(~c), size(a), size(b), &alpha, base_b, stride( b), base_a, stride(~a), &beta, base_c, stride(c));
		}
		else if(stride(~a)==1 and stride( b)==1 and not is_conjugated<A>{} and     is_conjugated<B>{}) ctx.gemm('C', 'N', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride( a), &beta, base_c, stride(c));
		else if(stride( a)==1 and stride( b)==1 and     is_conjugated<A>{} and     is_conjugated<B>{}) ctx.gemm('C', 'C', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride(~a), &beta, base_c, stride(c));
		else if(stride( a)==1 and stride( b)==1 and not is_conjugated<A>{} and     is_conjugated<B>{}) ctx.gemm('C', 'T', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride(~a), &beta, base_c, stride(c));
		else if(stride(~a)==1 and stride( b)==1 and not is_conjugated<A>{} and not is_conjugated<B>{}){
			if(size(a)==1) ctx.gemm('T', 'N', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, size(b)   , &beta, base_c, stride(c));
			else           ctx.gemm('T', 'N', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride( a), &beta, base_c, stride(c));
		}
		else if(stride( a)==1 and stride( b)==1 and     is_conjugated<A>{} and not is_conjugated<B>{}) ctx.gemm('T', 'C', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride(~a), &beta, base_c, stride(c));
		else if(stride( a)==1 and stride( b)==1 and not is_conjugated<A>{} and not is_conjugated<B>{}) ctx.gemm('T', 'T', size(~c), size(a), size(b), &alpha, base_b, stride(~b), base_a, stride(~a), &beta, base_c, stride(c));
		else                                                                                           assert(0&&" case not implemented in blas");
	}
	return std::forward<C>(c);
}

template<class A, class B, class C>
C&& gemm(typename A::element alpha, A const& a, B const& b, typename A::element beta, C&& c){
//->decltype(gemm('N', 'T', size(~c), size(a), size(b), &alpha, gemm_base_aux(b), stride( b), gemm_base_aux(a), stride(~a), &beta, gemm_base_aux(c), size(b)) , std::forward<C>(c)){
	using multi::blas::default_allocator_of;
	auto ctx = default_context_of(gemm_base_aux(a)); // ADL
	return gemm(ctx, alpha, a, b, beta, std::forward<C>(c));
}

template<class A2D, class B2D, class C2D = typename A2D::decay_type>
NODISCARD("because input arguments are const")
auto gemm(typename A2D::element a, A2D const& A, B2D const& B){
	assert(get_allocator(A) == get_allocator(B));
	return gemm(a, A, B, 0., C2D({size(A), size(rotated(B))}, get_allocator(A)));
}

template<class Context, class A2D, class B2D, class C2D = typename A2D::decay_type>
NODISCARD("because input arguments are const")
auto gemm(Context&& ctx, typename A2D::element a, A2D const& A, B2D const& B)
->std::decay_t<decltype(gemm(std::forward<Context>(ctx), a, A, B, 0., C2D({size(A), size(rotated(B))}, get_allocator(A))))>{
	assert(get_allocator(A) == get_allocator(B));
	return gemm(std::forward<Context>(ctx), a, A, B, 0., C2D({size(A), size(rotated(B))}, get_allocator(A)));
}

template<class A2D, class B2D> 
auto gemm(A2D const& A, B2D const& B)
->decltype(gemm(1., A, B)){
	return gemm(1., A, B);}

template<class Context, class A2D, class B2D> 
auto gemm(Context&& ctx, A2D const& A, B2D const& B)
->decltype(gemm(std::forward<Context>(ctx), 1., A, B)){
	return gemm(std::forward<Context>(ctx), 1., A, B);}

namespace operators{
	template<class A2D, class B2D> 
	auto operator*(A2D const& A, B2D const& B)
	->decltype(gemm(1., A, B)){
		return gemm(1., A, B);}
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_GEMM

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>
#include<random>

#include <boost/timer/timer.hpp>

#include "../blas/axpy.hpp"
#include "../blas/dot.hpp"
#include "../blas/nrm2.hpp"

namespace multi = boost::multi;

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j)
			cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<std::endl;
}

using complex = std::complex<double>; complex const I{0,1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1, 3},
		{ 9, 7},
	};
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
	};
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 148 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a,  b, 0.,  c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 169 and c[1][0] == 82 ));
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1.,  a, ~b, 0.,  c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 183 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a, ~b, 0.,  c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 117 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a, ~b, 0., ~c); // c⸆=a⸆b⸆, c=ba
		BOOST_REQUIRE( c[0][1] == 117 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1, 3, 1},
		{ 9, 7, 1},
	};
	multi::array<double, 2> const b = {
		{ 11, 12, 1},
		{  7, 19, 1},
		{  1,  1, 1}
	};
	{
		multi::array<double, 2> c({2, 3});
		using multi::blas::gemm;
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{2.+3.*I, 1.-2.*I}
	};
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::hermitized(a), 0., c); // c=aa†, c†=aa†
		BOOST_TEST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_TEST_REQUIRE( c[0][1] == 7.+10.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated){
	using complex = std::complex<double>; complex const I{0,1};
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I}
	};
	{
		multi::array<complex, 2> c({1, 1});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(a), 0., c); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[0][0] == 87. + 0.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bisbis){
	using complex = std::complex<double>; complex const I{0,1};
	multi::array<complex, 2> const a = {
		{1. + 2.*I}, 
		{9. - 1.*I}, 
		{1. + 1.*I}
	};
	multi::array<complex, 2> const b = {	
		{ 11. - 2.*I, 7. - 3.*I, 8. - 1.*I}
	};
	{
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({1, 1});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		using multi::blas::is_c_ordering;

		BOOST_REQUIRE( size(hermitized(a)) == 1 );
		BOOST_REQUIRE( size(hermitized(b)[0]) == 1 );

		gemm(1., hermitized(a), hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.+7.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_empty){
	multi::array<double, 2> const a({0, 5});
	BOOST_REQUIRE( size( a) == 0 );
	BOOST_REQUIRE( size(~a) == 5 );
	BOOST_REQUIRE( a.is_empty() );
	multi::array<double, 2> const b({5, 0});
	BOOST_REQUIRE( size( b) == 0 );
	BOOST_REQUIRE( size(~b) == 0 );
	BOOST_REQUIRE( b.is_empty() );
	{
		multi::array<double, 2> c;
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare2){
	multi::array<double, 2> const a = {
		{ 1, 3},
		{ 9, 7},
		{ 1, 1}
	};
	multi::array<double, 2> const b = {	
		{ 11, 12},
		{  7, 19}
	};
	{
		multi::array<double, 2> c({3, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );
	}
	{
		multi::array<double, 2> c({2, 3});
		using multi::blas::gemm;
		gemm(1., a, b, 0., ~c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );		
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({3, 2});
		using multi::blas::gemm;
		gemm(1., ~ar, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );		
	}
	{
		auto ar = rotated(a).decay();
		multi::array<double, 2> c({2, 3});
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., rotated(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );				
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x2){
	multi::array<double, 2> const a = {
		{ 1, 3},
		{ 9, 4},
		{ 1, 5}
	};
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., ~a, b, 0.,  c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 101 );
		gemm(1., ~a, b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 101 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x2){
	multi::array<double, 2> const a = {
		{1, 9, 1}
	};
	BOOST_TEST_REQUIRE( stride(~a) == 1 );
	BOOST_TEST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({size(a), size(~b)});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		auto ar = +~a;
		using multi::blas::gemm;
		multi::array<double, 2> c({size(~b), size(~ar)});
		gemm(1., ~ar, b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_3x2){
	multi::array<complex, 2> const a = {
		{1, 9, 1}
	};
	BOOST_TEST_REQUIRE( stride(~a) == 1 );
	BOOST_TEST_REQUIRE( stride( a) == 3 );
	multi::array<complex, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<complex, 2> c({size(a), size(~b)});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		auto ar = +~a;
		using multi::blas::gemm;
		multi::array<complex, 2> c({size(~b), size(~ar)});
		gemm(1., ~ar, b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_part_3x2){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 3, 3}
	};
	BOOST_TEST_REQUIRE( stride(~a) == 1 );
	BOOST_TEST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({size(a({0, 1})), size(~b)});
		using multi::blas::gemm;
		gemm(1., a({0, 1}), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		auto ar = +~a;
		using multi::blas::gemm;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_part_3x2){
	multi::array<complex, 2> const a = {
		{1., 9., 1.},
		{3., 3., 3.}
	};
	BOOST_TEST_REQUIRE( stride(~a) == 1 );
	BOOST_TEST_REQUIRE( stride( a) == 3 );
	multi::array<complex, 2> const b = {
		{ 11., 12.},
		{  7., 19.},
		{  8.,  1.}
	};
	{
		multi::array<complex, 2> c({size(a({0, 1})), size(~b)});
		using multi::blas::gemm;
		gemm(1., a({0, 1}), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		auto ar = +~a;
		using multi::blas::gemm;
		multi::array<complex, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184. );
	}
}


BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 3, 3}
	};
	BOOST_TEST_REQUIRE( stride(~a) == 1 );
	BOOST_TEST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11},
		{  7},
		{  8}
	};
	{
		multi::array<double, 2> c({size(a), size(~b)});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][0] == 82 );
		BOOST_REQUIRE( c[1][0] == 78 );
	}
	{
		auto ar = +~a;
		using multi::blas::gemm;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][0] == 82 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1_bis){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 4, 5}
	};
	multi::array<double, 2> const b = {
		{ 11},
		{  7},
		{  8}
	};

	{
		multi::array<double, 2> c({1, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., ~c); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][0] ==  82 );
		BOOST_REQUIRE( (~c)[1][0] == 101 );
	}
	{
		multi::array<double, 2> c({2, 1});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][1] == 101 );
		BOOST_REQUIRE(   c [1][0] == 101 );
	}
	{
		multi::array<double, 2> c({1, 2});
		auto ar = +~a;
		using multi::blas::gemm;
		gemm(1., ~ar, b, 0., ~c); // c⸆=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 101 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x1){
	multi::array<double, 2> const a = {
		{1, 9, 1}
	};
	multi::array<double, 2> const b = {	
		{ 11},
		{  7},
		{  8}
	};
	using multi::blas::gemm;
	{
		multi::array<double, 2> c({1, 1});
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto ar = +~a;
		gemm(1., ~ar, b, 0., c); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = +~b;
		gemm(1., a, ~br, 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = +~b;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(br), 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square){
	multi::array<complex, 2> const a = {
		{ 1.+3.*I, 3.+2.*I},
		{ 9.+1.*I, 7.+1.*I},
	};
	multi::array<complex, 2> const b = {
		{11.+2.*I, 12.+4.*I},
		{ 7.+1.*I, 19.-9.*I},
	};
	using multi::blas::gemm;
	{
		multi::array<complex, 2> c({2, 2});
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		gemm(1., ~a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 170.-8.*I and c[1][0] == 77.+42.*I ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		gemm(1., a, ~b, 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 177.+69.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using namespace multi::blas;
		gemm(1., T(a), T(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 109. + 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), transposed(b), 0., transposed(c)); // c⸆=a⸆b⸆, c=ba
		BOOST_REQUIRE( c[0][1] == 109.+68.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x1){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{1. + 2.*I, 9. - 1.*I, 1. + 1.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I},
		{  7. - 3.*I},
		{  8. - 1.*I}
	};
	using multi::blas::gemm;
	{
		multi::array<complex, 2> c({1, 1});
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto ar = +~a;
		gemm(1., ~ar, b, 0., c); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		using multi::blas::gemm;
		gemm(1., a, ~br, 0., c);
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(br), 0., ~c);
		BOOST_TEST_REQUIRE( c[0][0] == 80. + 53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_hermitized_square){
	multi::array<complex, 2> const a = {
		{ 1.+3.*I, 3.+2.*I},
		{ 9.+1.*I, 7.+1.*I},
	};
	multi::array<complex, 2> const b = {
		{11.+2.*I, 12.+4.*I},
		{ 7.+1.*I, 19.-9.*I},
	};
	namespace blas = multi::blas;
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c†=b†a†
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., blas::H(c)); // c†=a†b†, c=ba
		BOOST_REQUIRE( c[1][0] == 184. - 40.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 87. - 16.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I);
	}
	{
	//	multi::array<complex, 2> c({2, 2});
	//n	blas::gemm(1., blas::hermitized(a), blas::hermitized(b), 0., ~c); // c⸆=a†b†, c=b*a*
	//	BOOST_REQUIRE( c[0][1] == 109. - 68.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I},
		{9. - 1.*I},
		{1. + 1.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I},
		{  7. - 3.*I},
		{  8. - 1.*I}
	};
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto ha = +blas::hermitized(a);
		blas::gemm(1., ha, b, 0., c);
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
		blas::gemm(1., blas::H(b), a, 0., c);
		BOOST_REQUIRE( c[0][0] == 80.+53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x2){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I, 9. - 1.*I, 1. + 1.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I, 5. + 2.*I},
		{  7. - 3.*I, 2. + 1.*I},
		{  8. - 1.*I, 1. + 1.*I}
	};
	{
		multi::array<complex, 2> c({1, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 20.+21.*I );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({1, 2});
		blas::gemm(1., blas::H(ar), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x2){
	multi::array<complex, 2> const a = {
		{1. + 2.*I}, 
		{9. - 1.*I}, 
		{1. + 1.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I, 5. + 2.*I},
		{  7. - 3.*I, 2. + 1.*I},
		{  8. - 1.*I, 1. + 1.*I}
	};
	namespace blas = multi::blas;
	{
		auto ar = +~a;
		multi::array<complex, 2> c({1, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x2){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I, 5. + 2.*I}, 
		{9. - 1.*I, 9. + 1.*I}, 
		{1. + 1.*I, 2. + 2.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I, 5. + 2.*I},
		{  7. - 3.*I, 2. + 1.*I},
		{  8. - 1.*I, 1. + 1.*I}
	};
	{
		auto ar = +~a;
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x1){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I, 5. + 2.*I}, 
		{9. - 1.*I, 9. + 1.*I}, 
		{1. + 1.*I, 2. + 2.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I},
		{  7. - 3.*I},
		{  8. - 1.*I}
	};
	{
		auto ar = +~a;
		multi::array<complex, 2> c({2, 1});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bis){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I}, 
		{9. - 1.*I}, 
		{1. + 1.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I},
		{  7. - 3.*I},
		{  8. - 1.*I}
	};
	{
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80. - 53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square_automatic){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3.},
		{ 9., 7.},
	};
	multi::array<double, 2> const b = {
		{ 11., 12.},
		{  7., 19.},
	};
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 148 and c[1][1] == 241 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., a, blas::T(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][1] == 196. );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., blas::T(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == 169. and c[1][0] == 82. ));
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][1] == 154. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::array<double, 2> const b = {	
		{ 11., 12., 1.},
		{  7., 19., 1.},
		{  1.,  1., 1.}
	};
	{
		multi::array<double, 2> c({2, 3});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic){
	multi::array<complex, 2> const a = {
		{ 1. + 2.*I, 3. - 3.*I},
		{ 9. + 1.*I, 7. + 4.*I},
	};
	multi::array<complex, 2> const b = {
		{ 11. + 1.*I, 12. + 1.*I},
		{  7. + 8.*I, 19. - 2.*I},
	};
	namespace blas = multi::blas;
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::T(b), 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == complex(178, 75) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == complex(180, 29) and c[1][0] == complex(53, 54) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == complex(186, 65) and c[1][0] == complex(116, 25) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == complex(111, 64) and c[1][1] == complex(158, -51) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == complex(188, 43) and c[1][1] == complex(196, 25) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == complex(116, -25) and c[1][1] == complex(186, -65) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::H(b), 0., c); // c=a⸆b†, c†=ba⸆†
		BOOST_REQUIRE( c[1][0] == complex(118, 5) and c[1][1] == complex(122, 45) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == complex(116, 25) and c[1][1] == complex(186, 65) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	multi::array<complex, 2> const b = {
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
}

BOOST_AUTO_TEST_CASE(submatrix_result_issue_97){
	multi::array<complex, 2> M = {
		{2. + 3.*I, 2. + 1.*I, 1. + 2.*I},
		{4. + 2.*I, 2. + 4.*I, 3. + 1.*I},
		{7. + 1.*I, 1. + 5.*I, 0. + 3.*I}
	};
	
	multi::array<complex, 2> V = {
		{1. + 2.*I},
		{2. + 1.*I},
		{9. + 2.*I}
	};
	
	using multi::blas::gemm;
	using multi::blas::hermitized;
	using std::get;
	
	auto M2 = +M({0, 3}, {0, 1});
	BOOST_REQUIRE( M2 == M({0, 3}, {0, 1}) );
	
	BOOST_TEST( gemm(hermitized(M2               ), V)[0][0] == 83. + 6.*I );
	BOOST_TEST( gemm(hermitized(M({0, 3}, {0, 1})), V)[0][0] == 83. + 6.*I );

	using namespace multi::blas::operators;
	BOOST_REQUIRE( (hermitized(M)*V)[0][0] == 83. + 6.*I );
}

BOOST_AUTO_TEST_CASE(blas_context_gemm){

	auto rand = [d=std::normal_distribution<>{}, g=std::mt19937{}]()mutable{return complex{d(g), d(g)};};

	multi::array<complex, 2> A({30, 40});
	multi::array<complex, 2> B({40, 50});
	
	std::generate(A.elements().begin(), A.elements().end(), rand);
	std::generate(B.elements().begin(), B.elements().end(), rand);

	namespace blas = multi::blas;
	auto C = blas::gemm(A, B);

	using namespace multi::blas::operators;

	{
		auto sum = 0.;
		for(auto i : extension(~C))
			sum += blas::nrm2((~C)[i] - blas::gemv(A, (~B)[i]))();

		BOOST_TEST_REQUIRE(sum == 0, boost::test_tools::tolerance(1e-12));
	}
	
	BOOST_TEST_REQUIRE( std::inner_product(
		begin(~C), end(~C), begin(~B), 0., std::plus<>{}, [&A](auto const& Ccol, auto const& Bcol){
			return multi::blas::nrm2( Ccol - multi::blas::gemv(A, Bcol) );
	}) == 0. , boost::test_tools::tolerance(1e-12) );

	BOOST_TEST_REQUIRE( std::equal(
		begin(~C), end(~C), begin(~B), [&A](auto const& Ccol, auto const& Bcol){
			return multi::blas::nrm2( Ccol - multi::blas::gemv(A, Bcol) ) < 1e-12;
		}
	) );

	blas::context ctxt;
	auto C2 = gemm(ctxt, A, B);
	
	BOOST_TEST_REQUIRE( std::equal(
		begin(C), end(C), begin(C2), [](auto const& crow, auto const& c2row){return ((crow - c2row)^2) < 1e-13;}
	) );

}

#if 1
BOOST_AUTO_TEST_CASE(blas_gemm_timing){

	multi::array<complex, 2> A({1000, 2000});
	multi::array<complex, 2> B({2000, 3000});
	multi::array<complex, 2> C({size(A), size(~B)});
	multi::array<complex, 2> C2(extensions(C), complex{NAN, NAN});
	multi::array<complex, 2> C3(extensions(C), complex{NAN, NAN});
	multi::array<complex, 2> C4(extensions(C), complex{NAN, NAN});
	A[99][99] = B[11][22] = C[33][44] = 1.0;
	std::cerr<< "memory " << (A.num_elements()+ B.num_elements() + C.num_elements())*sizeof(complex)/1e6 <<" MB"<<std::endl;
	
	{
		boost::timer::auto_cpu_timer t;
		auto rand = [d=std::uniform_real_distribution<>{0., 10.}, g=std::mt19937{}]() mutable{return complex{d(g), d(g)};};
		std::generate(A.elements().begin(), A.elements().end(), rand);
		std::generate(B.elements().begin(), B.elements().end(), rand);
	}
	namespace blas = multi::blas;
	{
		boost::timer::auto_cpu_timer t; // 0.237581s
		C = blas::gemm(A, B);
	}
	{
		boost::timer::auto_cpu_timer t; // 4.516157s
		for(auto i : extension(~B)) (~C2)[i] = blas::gemv(A, (~B)[i]);
	}
	{
		boost::timer::auto_cpu_timer t; // 4.516157s
		for(auto i : extension(A)) C2[i] = blas::gemv(~B, A[i]);
	}
	{
		boost::timer::auto_cpu_timer t; // 32.705804s
		for(auto i:extension(A)) for(auto j:extension(~B)) C3[i][j] = blas::dot(A[i], (~B)[j]);
	}
	using namespace blas::operators;

	BOOST_TEST_REQUIRE( std::equal(
		begin(C), end(C), begin(C2), [](auto const& crow, auto const& c2row){
			return ((crow - c2row)^2) < 1e-13;
		}
	) );

	BOOST_TEST_REQUIRE( std::equal(
		begin(C), end(C), begin(C3), [](auto const& crow, auto const& c3row){
			return ((crow - c3row)^2) < 1e-13;
		}
	) );

}
#endif

#endif
#endif

