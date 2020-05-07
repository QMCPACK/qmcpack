#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
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

template<class X1D, class Y1D, class R, std::enable_if_t<not is_complex<X1D>{}, int> =0>
auto dot(X1D const& x, Y1D const& y, R&& r){
	return core::dot(size(x), base(x), stride(x), base(y), stride(y), &r), std::forward<R>(r);}

template<class X1D, class Y1D, class R, std::enable_if_t<    is_complex<X1D>{}, int> =0>
auto dot(X1D const& x, Y1D const& y, R&& r){

	auto base_x = dot_base_aux(x);
	auto base_y = dot_base_aux(y);

	using core::dotu;
	using core::dotc;

	     if(not is_conjugated<X1D>{} and not is_conjugated<Y1D>{}) dotu(size(x), base_x, stride(x), base_y, stride(y), &r);
	else if(not is_conjugated<X1D>{} and     is_conjugated<Y1D>{}) dotc(size(x), base_y, stride(y), base_x, stride(x), &r);
	else if(    is_conjugated<X1D>{} and not is_conjugated<Y1D>{}) dotc(size(x), base_x, stride(x), base_y, stride(y), &r);
	else                                                           assert(0);
	return std::forward<R>(r);
}

template<class X1D, class Y1D, class Alloc>
NODISCARD("when last argument is an allocator")
auto alloc_dot(X1D const& x, Y1D const& y, Alloc const& alloc){
	return dot(x, y, multi::array<typename X1D::value_type, 0, Alloc>(0, alloc) );
}

template<class X1D, class Y1D>
NODISCARD("when last argument is read-only")
auto dot(X1D const& x, Y1D const& y){
	return alloc_dot(x, y, common(get_allocator(x), get_allocator(y)));
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_DOT

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include "../blas/nrm2.hpp"

#include<cassert>
#include<numeric> // inner_product

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

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real){
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	using blas::dot;
	{
		double d = dot(cA[1], cA[2]);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		auto d2 = dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==d2 );
	}
	{
		multi::array<double, 0> d;
		auto d2 = dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = dot(cA[1], cA[2]);
		BOOST_REQUIRE( d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
		BOOST_REQUIRE( dot(cA[1], cA[2]) == dot(cA[2], cA[1]) );
	}
	{	
		using blas::nrm2;
		using std::sqrt;
		{
			double s;
			dot(cA[1], cA[1], s);
			assert( sqrt(s)==nrm2(cA[1]) );
		}
	}
}

#if 1
BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex){
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
//	print_2D(A);
//	print_1D(A[1]);
	{
		complex c; blas::dot(A[1], A[1], c);
		BOOST_TEST( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[1]);
		BOOST_TEST( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}) );
	}
	{
//		conjugated(A[1]);
//		complex c; dot(A[1], conjugated(A[1]), c);
//		BOOST_TEST( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
	}
#if 0


	{




		{

		{
			multi::array<complex, 1> cc = {1., 2., 3.};
			dot(A[1], conjugated(A[1]), cc[0]);
			BOOST_TEST( cc[0] == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
		}
		{
			auto const c = dot(A[1], conjugated(A[1]));
			std::cout<< c() <<std::endl;
			BOOST_REQUIRE( c() == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
			BOOST_REQUIRE( dot(A[1], conjugated(A[1])) == dot(conjugated(A[1]), A[1]) );
		}
		{
			auto const c = dot(conjugated(A[1]), A[1]);
			std::cout<< c() <<std::endl;
			BOOST_REQUIRE( c() == std::inner_product(begin(A[1]), end(A[1]), begin(A[1]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
		}
		{
			multi::array<complex, 1> a = {1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I};
			multi::array<complex, 1> b = {5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I};
			BOOST_REQUIRE( dot(a            , b            )()== 19. - 27.*I );
			BOOST_REQUIRE( dot(a            , conjugated(b))()==121. - 43.*I );
			BOOST_REQUIRE( dot(conjugated(a), b            )()==121. + 43.*I );
		//	BOOST_REQUIRE( dot(conjugated(a), conjugated(b))() == 19. + 27.*I );
		}
	}
#endif
}
#endif


#endif
#endif

