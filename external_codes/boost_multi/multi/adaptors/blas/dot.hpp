#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS -I/usr/local/cuda/include $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_DOT_HPP
#define MULTI_ADAPTORS_BLAS_DOT_HPP

#include "../blas/core.hpp"
#include "../blas/operations.hpp"

#include "../../array.hpp"
#include "../../config/NODISCARD.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0> 
auto dot_base_aux(A&& a){return base(a);}

template<class A, std::enable_if_t<  is_conjugated<A>{}, int> =0> 
auto dot_base_aux(A&& a){return underlying(base(a));}

using core::dot;

template<class X1D, class Y1D, class R, std::enable_if_t<not is_complex_array<X1D>{}, int> =0>
auto dot(X1D const& x, Y1D const& y, R&& r){
	return dot(size(x), base(x), stride(x), base(y), stride(y), &r), std::forward<R>(r);}

template<class X1D, class Y1D, class R, std::enable_if_t<    is_complex_array<X1D>{}, int> =0>
auto dot(X1D const& x, Y1D const& y, R&& r){

	auto base_x = dot_base_aux(x);
	auto base_y = dot_base_aux(y);

	using core::dotu;
	using core::dotc;

	     if(not is_conjugated<X1D>{} and not is_conjugated<Y1D>{}) dotu(size(x), base_x, stride(x), base_y, stride(y), &r);
	else if(not is_conjugated<X1D>{} and     is_conjugated<Y1D>{}) dotc(size(x), base_y, stride(y), base_x, stride(x), &r);
	else if(    is_conjugated<X1D>{} and not is_conjugated<Y1D>{}) dotc(size(x), base_x, stride(x), base_y, stride(y), &r);
//	else if(    is_conjugated<X1D>{} and     is_conjugated<Y1D>{}) dotu(size(x), base_y, stride(y), base_x, stride(x), &r);
	else                                                           assert(0); // case not implemented in blas
	return std::forward<R>(r);
}

template<class X1D, class Y1D, class Alloc>
NODISCARD("when last argument is an allocator")
auto alloc_dot(X1D const& x, Y1D const& y, Alloc const& alloc){
	return dot(x, y, multi::array<typename X1D::value_type, 0, Alloc>(0, alloc) );
}

template<class X1D, class Y1D>
NODISCARD("pure function")
auto dot(X1D const& x, Y1D const& y){
	return alloc_dot(x, y, common(get_allocator(x), get_allocator(y)));
}

namespace operators{

	template<class X1D, class Y1D>
	NODISCARD("no side effect")
	auto operator,(X1D const& x, Y1D const& y)
	->decltype(dot(x, y)){
		return dot(x, y);}

}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_DOT

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include "../blas/nrm2.hpp"

#include<cassert>
#include<numeric> // inner_product

#include<thrust/complex.h>

namespace multi = boost::multi;
namespace blas = multi::blas;

template<class M, typename = decltype(std::declval<M const&>()[0]), typename = decltype(std::declval<M const&>()[0][0])> 
decltype(auto) print_2D(M const& C){
	using std::cout;
	using multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j)
			cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<std::endl;
}

template<class M, typename = decltype(std::declval<M const&>()[0])>//, typename = decltype(std::declval<M const&>()[0])>
decltype(auto) print_1D(M const& C){
	using std::cout; using multi::size;
	for(int i = 0; i != size(C); ++i) cout<< C[i] <<' ';
	return cout<<std::endl;
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real){

	multi::array<float, 1> V = {1., 2., 3.};
	multi::array<float, 1> W = {1., 2., 3.};
	
	using blas::dot;
	BOOST_REQUIRE( 14. == dot(V, W) );
	BOOST_REQUIRE( dot(V, W) == 14. );

}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real){
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	{
		double d = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		blas::dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		auto d2 = blas::dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==d2 );
	}
	{
		multi::array<double, 0> d;
		auto d2 = blas::dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
		BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
	}
	{	
		double s;
		blas::dot(cA[1], cA[1], s);
		BOOST_REQUIRE( std::sqrt(s)==blas::nrm2(cA[1]) );
	}
	{
		auto d1 = blas::dot(cA[1], cA[1]);
	//	auto d2 = blas::dot(blas::conj(cA[1]), cA[1]);
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex){
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
		complex c; blas::dot(A[1], A[2], c);
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], blas::C(A[2]));
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
	}
	{
		complex c = blas::dot(blas::conj(A[1]), A[2]);
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
	}
	{
	//	complex c = blas::dot(blas::C(A[1]), blas::C(A[2]));
	//	BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto a, auto b){return conj(a)*conj(b);}) );
	}

	{
	//	complex c = blas::dot(blas::C(A[1]), A[2]);
	//	BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
	}
	{
//		complex c = blas::dot(blas::C(A[1]), blas::C(A[2]));
//		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return conj(a)*conj(b);}) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case){
	multi::array<double, 1> v1(10, +1.0);
	multi::array<double, 1> v2(10, -1.0);

	using blas::dot;
	using blas::hermitized;
	using blas::conj;
	
	auto a = dot(v1, v2);
	auto b = dot(hermitized(v1), v2);
	
	BOOST_REQUIRE(a == b);
	
	auto c = dot(blas::conj(v1), v2); // conjugation doesn't do anything for real array
	BOOST_REQUIRE(c == a);
	
	auto d_arr = dot(blas::C(v1), v2);
	BOOST_REQUIRE(d_arr == a);
	
	static_assert( not std::is_same<decltype(d_arr), double>{}, "!" );

	using blas::C;
	double d_doub = dot(C(v1), v2);
	
	BOOST_REQUIRE( d_doub == d_arr );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex_thrust){
	namespace blas = multi::blas;

	using complex = thrust::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
		complex c; blas::dot(A[1], A[2], c);
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], blas::C(A[2]));
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
	}
}

#endif
#endif

