// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/blas/gemm.hpp"

#include "config.hpp"

#include "../../../array.hpp"

#include<random>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_gemm_square_real){
	multi::array<double, 2> const a = {
		{1, 3, 4},
		{9, 7, 1},
		{1, 2, 3}
	};
	multi::array<double, 2> const b = {
		{11, 12, 4},
		{ 7, 19, 1},
		{11, 12, 4}
	};
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., a, b, 0., c);
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		BOOST_REQUIRE( size( a) == size( c) );
		BOOST_REQUIRE( size(~b) == size(~c) );
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 86 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., a, blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm_n(1., a.begin(), a.size(), blas::T(b).begin(), 0., c.begin());
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., blas::T(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., blas::T(a), blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 50 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 50 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., a, blas::T(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm_n(1., begin(a), size(a), begin(blas::T(b)), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 48 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm(1., blas::T(a), b, 0., c);
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)}, 9999);
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 103 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm(2., blas::H(a), blas::H(b), 0., c);
		BOOST_REQUIRE( c[2][1] == 100 );
	}
	{
		multi::array<double, 2> c = blas::gemm(2., blas::H(a), blas::H(b));
		BOOST_REQUIRE( c[2][1] == 100 );
	}
	{
		multi::array<double, 2> const c = blas::gemm(2., blas::H(a), blas::H(b));
		BOOST_REQUIRE( c[2][1] == 100 );
	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		c = blas::gemm(2., blas::H(a), blas::H(b));
		BOOST_REQUIRE( c[2][1] == 100 );
	}
	{
		multi::array<double, 2> c;
		c = blas::gemm(2., blas::H(a), blas::H(b));
		BOOST_REQUIRE( c[2][1] == 100 );
	}
//	{
//		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
//		blas::gemm(2., blas::H(a), blas::H(b), 0., c);
//		BOOST_REQUIRE( c[2][1] == 100 );

//		multi::array<double, 2> const c_copy = blas::gemm(2., blas::H(a), blas::H(b));
//		BOOST_REQUIRE( c == c_copy );
//		multi::array<double, 2> const c_copy2 = blas::gemm(1., blas::H(a), blas::H(b));
//		BOOST_REQUIRE( c_copy2[2][1] == 50 );
//	}
	{
		multi::array<double, 2> c({size(a), size(rotated(b))}, 9999);
		blas::gemm_n(2., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c));
		BOOST_REQUIRE( c[2][1] == 100 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square){
	multi::array<double, 2> const a = {
		{  1,  3},
		{  9,  7},
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
		blas::context ctxt;
		blas::gemm_n(ctxt, 1., begin(a), size(a), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[1][0] == 148 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a,  b, 0.,  c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 169 and c[1][0] == 82 ));
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::context ctxt;
		blas::gemm_n(ctxt, 1., begin(~a), size(~a), begin(b), 0., begin( c));
		BOOST_REQUIRE(( c[1][1] == 169 and c[1][0] == 82 ));
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::context ctxt;
		blas::gemm_n(ctxt, 1., begin(~a), size(~a), begin(b), 0., begin(~c));
		BOOST_REQUIRE( (~c)[1][1] == 169 );
		BOOST_REQUIRE( (~c)[1][0] ==  82 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1.,  a, ~b, 0.,  c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 183 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::context ctxt;
		blas::gemm_n(ctxt, 1., begin(a), size(a), begin(~b), 0., begin(c)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 183 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1.,  a, ~b, 0.,  ~c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( (~c)[1][0] == 183 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(~b), 0., begin(~c)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( (~c)[1][0] == 183 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a, ~b, 0.,  c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 117 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm_n(1., begin(~a), size(~a), begin(~b), 0., begin(c)); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 117 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a, ~b, 0., ~c); // c⸆=a⸆b⸆, c=ba
		BOOST_REQUIRE( c[0][1] == 117 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm_n(1., begin(~a), size(~a), begin(~b), 0., begin(~c)); // c⸆=a⸆b⸆, c=ba
		BOOST_REQUIRE( c[0][1] == 117 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare){
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
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17 );
	}
	{
		multi::array<double, 2> c({2, 3});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic){//, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::array<double, 2> const b = {	
		{ 11., 12., 4., 8.},
		{  7., 19., 2., 7.},
		{  5.,  3., 3., 1.}
	};
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm_n(0.1, begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
	{
		auto c =+ blas::gemm(0.1, a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
	{
		multi::array c = blas::gemm(0.1, a, b);
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh){
	using complex = std::complex<double>; complex const I{0,1};
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{2.+3.*I, 1.-2.*I}
	};
	{
		auto c =+ blas::gemm(1., a, blas::H(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array c = blas::gemm(1., a, blas::H(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c = blas::gemm(1., a, blas::H(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		c = blas::gemm(1., a, blas::H(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		c() = blas::gemm(1., a, blas::H(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::gemm(1., a, blas::H(a), 0., c); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(a)), 0., begin(c)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
}

#if CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh_thrust){
	using complex = thrust::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{2.+3.*I, 1.-2.*I}
	};
	{
		auto c =+ blas::gemm(1., a, blas::hermitized(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array c = blas::gemm(1., a, blas::hermitized(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c = blas::gemm(1., a, blas::hermitized(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		c = blas::gemm(1., a, blas::hermitized(a)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::hermitized(a), 0., c); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(a)), 0., begin(c)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
	}
}
#endif

BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I}
	};
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, blas::H(a), 0., c); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[0][0] == 87. + 0.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(a)), 0., begin(c)); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[0][0] == 87. + 0.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bisbis){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{1. + 2.*I},
		{9. - 1.*I},
		{1. + 1.*I}
	};
	multi::array<complex, 2> const b = {	
		{ 11. - 2.*I, 7. - 3.*I, 8. - 1.*I}
	};
	{
		multi::array<complex, 2> c({1, 1});

		BOOST_REQUIRE( size(blas::H(a)) == 1 );
		BOOST_REQUIRE( size(blas::H(b)[0]) == 1 );

		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.+7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
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
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
	}
	{
		multi::array<double, 2> c;
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
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
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );
	}
	{
		multi::array<double, 2> c({size(~b), size(a)});
		blas::gemm(1., a, b, 0., ~c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );
	}
	{
		multi::array<double, 2> c({size(~b), size(a)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(~c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({3, 2});
		blas::gemm(1., ~ar, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({3, 2});
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[2][1] == 31 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({2, 3});
		blas::gemm(1., ~ar, b, 0., ~c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({2, 3});
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(~c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x2_2x2){
	multi::array<double, 2> const a = {
		{ 1, 3},
		{ 9, 4},
	};
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
	};
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm(1., ~a, b, 0.,  c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 61 );

		blas::gemm(1., ~a, b, 0., ~c); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[0][1] == 61 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin( c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 61 );

		blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin(~c)); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[0][1] == 61 );
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
		blas::gemm(1., ~a, b, 0.,  c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 101 );

		blas::gemm(1., ~a, b, 0., ~c); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[0][1] == 101 );
	}
	{
		multi::array<double, 2> c({2, 2});
		blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin( c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 101 );

		blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin(~c)); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[0][1] == 101 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x2){
	multi::array<double, 2> const a = {
		{1, 9, 1}
	};
	BOOST_REQUIRE( stride(~a) == 1 );
	BOOST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({size(~b), size(~ar)});
		blas::gemm(1., ~ar, b, 0., ~c); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
	{
		auto ar = +~a;
		BOOST_REQUIRE( size(~ar) == 1 );
		BOOST_REQUIRE( begin(~ar).stride() == 1 );
		BOOST_REQUIRE( begin(~ar)->stride() == 1 );
		BOOST_REQUIRE( begin( ar)->stride() == 1 );

		multi::array<double, 2> c({size(~b), size(~ar)});
		BOOST_REQUIRE( begin( c).stride() == 1 );
		BOOST_REQUIRE( begin(~c).stride() == 1 );
		BOOST_REQUIRE( begin(c)->stride() == 1 );
		
		BOOST_REQUIRE( begin(b) );
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(~c)); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_3x2){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a = {
		{1, 9, 1}
	};
	BOOST_REQUIRE( stride(~a) == 1 );
	BOOST_REQUIRE( stride( a) == 3 );
	multi::array<complex, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<complex, 2> c({size(a), size(~b)});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		multi::array<complex, 2> c({size(a), size(~b)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({size(~b), size(~ar)});
		blas::gemm(1., ~ar, b, 0., ~c); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[1][0] == 184. );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({size(~b), size(~ar)});
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(~c)); // c⸆=a⸆b, c=b⸆a
		BOOST_REQUIRE( c[1][0] == 184. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_part_3x2){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 3, 3}
	};
	BOOST_REQUIRE( stride(~a) == 1 );
	BOOST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({size(a({0, 1})), size(~b)});
		blas::gemm(1., a({0, 1}), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		multi::array<double, 2> c({size(a({0, 1})), size(~b)});
		blas::gemm_n(1., begin(a({0, 1})), size(a({0, 1})), begin(b), 0., begin(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm_n(1., begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0., begin(~c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complexreal_1x3_part_3x2){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a = {
		{1., 9., 1.},
		{3., 3., 3.}
	};
	BOOST_REQUIRE( stride(~a) == 1 );
	BOOST_REQUIRE( stride( a) == 3 );
	multi::array<complex, 2> const b = {
		{ 11., 12.},
		{  7., 19.},
		{  8.,  1.}
	};
	{
		multi::array<complex, 2> c({size(a({0, 1})), size(~b)});
		blas::gemm(1., a({0, 1}), b, 0., c);
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		multi::array<complex, 2> c({size(a({0, 1})), size(~b)});
		blas::gemm_n(1., begin(a({0, 1})), size(a({0, 1})), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[0][1] == 184. );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c);
		BOOST_REQUIRE( c[1][0] == 184. );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm_n(1., begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0., begin(~c));
		BOOST_REQUIRE( c[1][0] == 184. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 3, 3}
	};
	BOOST_REQUIRE( stride(~a) == 1 );
	BOOST_REQUIRE( stride( a) == 3 );
	multi::array<double, 2> const b = {
		{ 11},
		{  7},
		{  8}
	};
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 82 );
		BOOST_REQUIRE( c[1][0] == 78 );
	}
	{
		multi::array<double, 2> c({size(a), size(~b)});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][0] == 82 );
		BOOST_REQUIRE( c[1][0] == 78 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm(1., ~(ar(extension(ar), {0, 1})), b, 0., ~c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		auto ar = +~a;
		multi::array<double, 2> c({size(~b), size(~ar(extension(ar), {0, 1}))});
		blas::gemm_n(1., begin(~(ar(extension(ar), {0, 1}))), size(~(ar(extension(ar), {0, 1}))), begin(b), 0., begin(~c)); // c=a⸆b, c⸆=b⸆a
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
		blas::gemm(1., a, b, 0., ~c); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][0] ==  82 );
		BOOST_REQUIRE( (~c)[1][0] == 101 );
	}
	{
		multi::array<double, 2> c({1, 2});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(~c)); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][0] ==  82 );
		BOOST_REQUIRE( (~c)[1][0] == 101 );
	}
	{
		multi::array<double, 2> c({2, 1});
		blas::gemm(1., a, b, 0., c); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][1] == 101 );
		BOOST_REQUIRE(   c [1][0] == 101 );
	}
	{
		multi::array<double, 2> c({2, 1});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( (~c)[0][1] == 101 );
		BOOST_REQUIRE(   c [1][0] == 101 );
	}
	{
		multi::array<double, 2> c({1, 2});
		auto ar = +~a;
		blas::gemm(1., ~ar, b, 0., ~c); // c⸆=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 101 );
	}
	{
		multi::array<double, 2> c({1, 2});
		auto ar = +~a;
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(~c)); // c⸆=ab, c⸆=b⸆a⸆
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
	{
		multi::array<double, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto ar = +~a;
		blas::gemm(1., ~ar, b, 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto ar = +~a;
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = +~b;
		blas::gemm(1., a, ~br, 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		BOOST_REQUIRE( begin(c). stride() == 1 );
		BOOST_REQUIRE( begin(c)->stride() == 1 );

		auto br = +~b;
	//	BOOST_REQUIRE( begin(br). stride() == 1 );
		BOOST_REQUIRE( begin( br)->stride() == 1 );

		BOOST_REQUIRE(begin(a)->stride() == 1);
		BOOST_REQUIRE( begin(~br). stride() == 1 );
	//	BOOST_REQUIRE( begin(~br)->stride() == 1 );
		BOOST_REQUIRE(begin(c)->stride() == 1);
		BOOST_REQUIRE(begin(c).stride() == 1);
		BOOST_REQUIRE(size(a) == 1);

		blas::gemm_n(1., begin(a), size(a), begin(~br), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = +~b;
		blas::gemm(1., a, blas::H(br), 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = +~b;
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(br)), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 82 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square){
	using complex = std::complex<double>; constexpr complex I{0, 1};
	multi::array<complex, 2> const a = {
		{ 1.+3.*I, 3.+2.*I},
		{ 9.+1.*I, 7.+1.*I},
	};
	multi::array<complex, 2> const b = {
		{11.+2.*I, 12.+4.*I},
		{ 7.+1.*I, 19.-9.*I},
	};
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., ~a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 170.-8.*I and c[1][0] == 77.+42.*I ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(~a), size(~a), begin(b), 0., begin(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 170.-8.*I and c[1][0] == 77.+42.*I ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, ~b, 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 177.+69.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(~b), 0., begin(c)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 177.+69.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 109. + 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c)); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 109. + 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., blas::T(c)); // c⸆=a⸆b⸆, c=ba
		BOOST_REQUIRE( c[0][1] == 109.+68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(blas::T(c))); // c⸆=a⸆b⸆, c=ba
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
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto ar = +~a;
		blas::gemm(1., ~ar, b, 0., c); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto ar = +~a;
		blas::gemm_n(1., begin(~ar), size(~ar), begin(b), 0., begin(c)); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		blas::gemm(1., a, ~br, 0., c);
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		blas::context ctxt;
		blas::gemm_n(ctxt, 1., begin(a), size(a), begin(~br), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		blas::gemm(1., a, blas::H(br), 0., ~c);
		BOOST_REQUIRE( c[0][0] == 80. + 53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = +~b;
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(br)), 0., begin(~c));
		BOOST_REQUIRE( c[0][0] == 80. + 53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_hermitized_square){
	using complex = std::complex<double>; constexpr complex I{0, 1};
	multi::array<complex, 2> const a = {
		{ 1.+3.*I, 3.+2.*I},
		{ 9.+1.*I, 7.+1.*I},
	};
	multi::array<complex, 2> const b = {
		{11.+2.*I, 12.+4.*I},
		{ 7.+1.*I, 19.-9.*I},
	};
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c†=b†a†
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c†=b†a†
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c)); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., blas::H(c)); // c†=a†b†, c=ba
		BOOST_REQUIRE( c[1][0] == 184. - 40.*I );
	}
//	{
//		multi::array<complex, 2> c({2, 2});
//		blas::context ctxt;
//		blas::gemm_n(ctxt, 1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(blas::H(c))); // c†=a†b†, c=ba
//		BOOST_REQUIRE( c[1][0] == 184. - 40.*I );
//	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 87. - 16.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 87. - 16.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		c = blas::gemm(1., a, blas::H(b)); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array<complex, 2> c = blas::gemm(1., a, blas::H(b)); // c=ab†, c†=ba†
		BOOST_REQUIRE( size(c) == 2 );
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array c = blas::gemm(1., a, blas::H(b)); // CTAD
		BOOST_REQUIRE( size(c) == 2 );
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		auto c = multi::array<complex, 2>(blas::gemm(1., a, blas::H(b))); // c=ab†, c†=ba†
		BOOST_REQUIRE( size(c) == 2 );
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c)); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I);
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c)); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I);
	}
//	{
//		multi::array<complex, 2> c({2, 2});
//		blas::gemm(1., blas::H(a), blas::H(b), 0., ~c); // case no implemented in blas
//		BOOST_REQUIRE( c[0][1] == 109. - 68.*I );
//	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1){
	using complex = std::complex<double>; complex const I{0, 1};
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
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., blas::H(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=a†b, c†=b†a
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
	{
		multi::array<complex, 2> c({1, 1});
		auto ha = +blas::hermitized(a);
		blas::gemm_n(1., begin(ha), size(ha), begin(b), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );

		blas::gemm_n(1., begin(blas::H(b)), size(blas::H(b)), begin(a), 0., begin(c));
		BOOST_REQUIRE( c[0][0] == 80.+53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x2){
	using complex = std::complex<double>; constexpr complex I{0, 1};
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
		multi::array<complex, 2> c({1, 2});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 20.+21.*I );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({1, 2});
		blas::gemm(1., blas::H(ar), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
	{
		auto ar = +~a;
		multi::array<complex, 2> c({1, 2});
		blas::gemm_n(1., begin(blas::H(ar)), size(blas::H(ar)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x2){
	using complex = std::complex<double>; complex const I{0, 1};
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
	{
		multi::array<complex, 2> c({1, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
	{
		multi::array<complex, 2> c({1, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x2){
	using complex = std::complex<double>; complex const I{0, 1};
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
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x1){
	using complex = std::complex<double>; complex const I{0, 1};
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
		multi::array<complex, 2> c({2, 1});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
	{
		multi::array<complex, 2> c({2, 1});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bis){
	using complex = std::complex<double>; complex const I{0, 1};
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
		BOOST_REQUIRE( c[0][0] == 80. - 53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80. - 53.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square_automatic){
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
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
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

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic){
	using complex = std::complex<double>; complex const I{0, 1};
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
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, blas::T(b), 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == complex(178, 75) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(blas::T(b)), 0., begin(c)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == complex(178, 75) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == complex(180, 29) and c[1][0] == complex(53, 54) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(b), 0., begin(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == complex(180, 29) and c[1][0] == complex(53, 54) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == complex(186, 65) and c[1][0] == complex(116, 25) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == complex(186, 65) and c[1][0] == complex(116, 25) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == complex(111, 64) and c[1][1] == complex(158, -51) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(b), 0., begin(c)); // c=a†b, c†=b†a
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
		blas::gemm_n(1., begin(blas::H(a)), size(blas::H(a)), begin(blas::H(b)), 0., begin(c)); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == complex(116, -25) and c[1][1] == complex(186, -65) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::H(b), 0., c); // c=a⸆b†, c†=ba⸆†
		BOOST_REQUIRE( c[1][0] == complex(118, 5) and c[1][1] == complex(122, 45) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::H(b)), 0., begin(c)); // c=a⸆b†, c†=ba⸆†
		BOOST_REQUIRE( c[1][0] == complex(118, 5) and c[1][1] == complex(122, 45) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::T(a), blas::T(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == complex(116, 25) and c[1][1] == complex(186, 65) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm_n(1., begin(blas::T(a)), size(blas::T(a)), begin(blas::T(b)), 0., begin(c)); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == complex(116, 25) and c[1][1] == complex(186, 65) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
	using complex = std::complex<double>; complex const I{0, 1};
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
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_realcomplex_complex_nonsquare_automatic){
	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const a = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};
	multi::array<complex, 2> const b = {
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	{
		multi::array<complex, 2> c = blas::gemm(1., a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(53, 24) );
	}
	{
		multi::array<complex, 2> c({2, 4});
		c = blas::gemm(1., a, b); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(53, 24) );
	}
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(53, 24) );
	}
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm_n(1., begin(a), size(a), begin(b), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(53, 24) );
	}
	{
		multi::array<double, 2> const a_real = {
			{ 1., 3., 1.},
			{ 9., 7., 1.},
		};
		multi::array<complex, 2> c({2, 4});
		blas::real_doubled(c) = blas::gemm(1., a_real, blas::real_doubled(b));

		BOOST_REQUIRE( c[1][2] == complex(53, 24) );
	}
}

BOOST_AUTO_TEST_CASE(submatrix_result_issue_97){
	using complex = std::complex<double>; constexpr complex I{0, 1};
	multi::array<complex, 2> M = {
		{2. + 3.*I, 2. + 1.*I, 1. + 2.*I},
		{4. + 2.*I, 2. + 4.*I, 3. + 1.*I},
		{7. + 1.*I, 1. + 5.*I, 0. + 3.*I}
	};
	
	
	auto M2 = +M({0, 3}, {0, 1});
	BOOST_REQUIRE( M2 == M({0, 3}, {0, 1}) );

//	multi::array<complex, 2> V = {
//		{1. + 2.*I},
//		{2. + 1.*I},
//		{9. + 2.*I}
//	};

//	BOOST_REQUIRE( (+blas::gemm(1., blas::H(M2               ), V))[0][0] == 83. + 6.*I );
//	BOOST_REQUIRE( (+blas::gemm(1., blas::H(M({0, 3}, {0, 1})), V))[0][0] == 83. + 6.*I );

//	using namespace multi::blas::operators;
//	BOOST_REQUIRE( (+(blas::H(M)*V))[0][0] == 83. + 6.*I );
}


BOOST_AUTO_TEST_CASE(blas_context_gemm){
	using complex = std::complex<double>; static constexpr complex I{0, 1};
	auto rand = [d=std::normal_distribution<>{}, g=std::mt19937{}]()mutable{return d(g) + d(g)*I;}; // NOLINT(cert-msc32-c, cert-msc51-cpp): test purposes

	multi::array<complex, 2> A({30, 40});
	multi::array<complex, 2> B({40, 50});
	
	std::generate(A.elements().begin(), A.elements().end(), rand);
	std::generate(B.elements().begin(), B.elements().end(), rand);

//	auto C = +blas::gemm(1., A, B);

//	using namespace multi::blas::operators;

//	{
//		auto sum = 0.;
//		for(auto i : extension(~C))
//			sum += blas::nrm2((~C)[i] - blas::gemv(A, (~B)[i]))();

//		BOOST_REQUIRE(sum == 0, boost::test_tools::tolerance(1e-12));
//	}
	
//	BOOST_REQUIRE( std::inner_product(
//		begin(~C), end(~C), begin(~B), 0., std::plus<>{}, [&A](auto const& Ccol, auto const& Bcol){
//			return multi::blas::nrm2( Ccol - multi::blas::gemv(A, Bcol) );
//	}) == 0. , boost::test_tools::tolerance(1e-12) );

//	BOOST_REQUIRE( std::equal(
//		begin(~C), end(~C), begin(~B), [&A](auto const& Ccol, auto const& Bcol){
//			return multi::blas::nrm2( Ccol - multi::blas::gemv(A, Bcol) ) < 1e-12;
//		}
//	) );

//	blas::context ctxt;
//	auto C2 = +blas::gemm(&ctxt, 1., A, B);
	
//	BOOST_REQUIRE( std::equal(
//		begin(C), end(C), begin(C2), [](auto const& crow, auto const& c2row){return ((crow - c2row)^2) < 1e-13;}
//	) );

}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_gemm_range){//, *utf::tolerance(0.00001)){
	multi::array<double, 2> const a({2, 3}, 0.);
	multi::array<double, 2> const b({4, 3}, 0.);
	{
		multi::array<double, 2> c({2, 4});
		c() = blas::gemm(0.1, a, blas::H(b));
		BOOST_REQUIRE_CLOSE( c[1][2],  0., 0.00001 );
	}
	{
		multi::array<double, 2> c = blas::gemm(0.1, a, blas::H(b)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][2] == 0. );
	}
	{
		multi::array<double, 2> const a = {
			{1, 3, 1},
			{9, 7, 1},
		};
		(void)a;
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_hermitized_second_gemm_range){//, *utf::tolerance(0.00001)){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a({2, 3}, 0.);
	multi::array<complex, 2> const b({4, 3}, 0.);
	{
		multi::array<complex, 2> c({2, 4}, 999.);
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
		BOOST_REQUIRE( c[1][2] != 999. );
	}
	{
		multi::array<complex, 2> c = blas::gemm(1., a, blas::H(b)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][2] == 0. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second){//, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	multi::array<double, 2> const b = {
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53. );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53. );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
	{
		multi::array<double, 2> c({2, 4});
		blas::gemm_n(0.1, begin(a), size(a), begin(blas::H(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
	{
		multi::array<double, 2> c({2, 4});
		c() = blas::gemm(0.1, a, blas::H(b));
	}
	{
		multi::array<double, 2> c = blas::gemm(0.1, a, blas::H(b)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE_CLOSE( c[1][2] , 5.3 , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_real_nonsquare_hermitized_second){//, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const a = {
		{1., 3., 1.},
		{9., 7., 1.},
	};
	multi::array<complex, 2> const b = {
		{11.,  7., 5.},
		{12., 19., 3.},
		{ 4.,  2., 3.},
		{ 8.,  7., 1.}
	};
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53. );
	}
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53. );
	}
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( real(c[1][2]) , 5.3 , 0.00001 );
	}
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm_n(0.1, begin(a), size(a), begin(blas::H(b)), 0., begin(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE_CLOSE( real(c[1][2]) , 5.3 , 0.00001 );
	}
	{
		multi::array<complex, 2> c({2, 4});
		c() = blas::gemm(0.1, a, blas::H(b));
	}
	{
		multi::array<complex, 2> c = blas::gemm(0.1, a, blas::H(b)); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE_CLOSE( real(c[1][2]) , 5.3 , 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(blas_gemm_1xn_complex){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a({1, 100}, 1.);
	multi::array<complex, 2> const b({1, 100}, 1.);

	multi::array<complex, 2> c({1, 1}, 999.);
	blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
	BOOST_REQUIRE( c[0][0] == 100. );
}

BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_complex_inq_hydrogen_case){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a({100, 1}, 2.);
	multi::array<complex, 2> const b({1,   1}, 3.);

	multi::array<complex, 2> c({100, 1}, 999.);
	blas::gemm_n(1., begin(a), size(a), begin(blas::H(b)), 0., begin(c));
	BOOST_REQUIRE( c[0][0] == 6. );
	BOOST_REQUIRE( c[1][0] == 6. );
}

BOOST_AUTO_TEST_CASE(blas_gemm_nx1_times_1x1_1x1_complex_inq_hydrogen_case){
	using complex = std::complex<double>;
	multi::array<complex, 2> const a({1, 1}, 2.);
	multi::array<complex, 2> const b({1, 1}, 3.);

	multi::array<complex, 2> c({1, 1}, 999.);
	c = blas::gemm(1., a, b);
	BOOST_REQUIRE( c[0][0] == 6. );
}

BOOST_AUTO_TEST_CASE(blas_gemm_inq_case){ // https://gitlab.com/correaa/boost-multi/-/issues/97

	using complex = std::complex<double>;
	multi::array<complex, 2> mat({10, 2},  1.0);
	multi::array<complex, 2> vec({10, 1}, -2.0);

	mat({0, 10}, {1, 2}) = vec;
	
	namespace blas = multi::blas;

	{
		auto olap1 =+ blas::gemm(1., blas::H(mat)                 , vec);
		auto olap2 =+ blas::gemm(1., blas::H(mat({0, 10}, {0, 1})), vec);

		multi::array<complex, 2> mat2 = mat({0, 10}, {0, 1});
		auto olap3 =+ blas::gemm(1., blas::H(mat2), vec);

		BOOST_REQUIRE(olap1[0][0] == olap2[0][0]);
		BOOST_REQUIRE(olap3[0][0] == olap2[0][0]);
	}
	{
		multi::array<complex, 2> mat2 = mat({0, 3}, {0, 1});
		auto olap3 =+ blas::gemm(1., blas::H(mat({0, 3}, {0, 1})), vec);
		BOOST_REQUIRE( (+blas::gemm(1., blas::H(mat2), vec))[0][0] == (+blas::gemm(1., blas::H(mat({0, 3}, {0, 1})), vec))[0][0] );
	}

}

#if 0
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_gpu){//, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	cuda::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	using multi::blas::gemm;using multi::blas::hermitized;
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	auto c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	cuda::array<double, 2> c;
	//	c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto c = blas::gemm(0.1, a, blas::H(b)); // c=ab, c⸆=b⸆a⸆
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_hermitized_second_managed){//, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::managed::array<double, 2> const a = {
		{1, 3, 1},
		{9, 7, 1},
	};
	cuda::managed::array<double, 2> const b = {	
		{11,  7, 5},
		{12, 19, 3},
		{ 4,  2, 3},
		{ 8,  7, 1}
	};
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(1., a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 53 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
		multi::cuda::managed::array<double, 2> c({2, 4});
		auto c_copy = blas::gemm(0.1, a, blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c_copy[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	auto c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
	//	auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
	//	cuda::managed::array<double, 2> c;
	//	c = f(a, b);
	//	BOOST_REQUIRE( c[1][2] == 5.3 );
	}
	{
		auto f = [](auto&& a, auto&& b){return blas::gemm(0.1, a, blas::H(b));};
		multi::cuda::managed::array<double, 2> c = a; BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 3);
		c = f(a, b);
		BOOST_REQUIRE( c[1][2] == 5.3 ); BOOST_REQUIRE(size(c) == 2 and size(rotated(c)) == 4);
	}
	{
		auto c = blas::gemm(0.1, a, blas::H(b)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 5.3 );
	}
}

using complex = std::complex<double>; complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
	namespace cuda = multi::cuda;
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
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		blas::gemm(complex(1.), acu, bcu, complex(0.), ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., amcu, bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}

struct multiplies_bind1st{
	multiplies_bind1st(multi::cuda::managed::array<complex, 2>&& m) : m_(std::move(m)){}
	template<class A>
	auto operator()(A const& a) const{
		using multi::blas::gemm;
		return gemm(m_, a);
	}
private:
	multi::cuda::managed::array<complex, 2> m_;
};

BOOST_AUTO_TEST_CASE(multi_constructors_inqnvcc_bug){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	cuda::managed::array<complex, 2> m = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	cuda::managed::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	auto c = blas::gemm(m, b);
	BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	BOOST_REQUIRE( b[1][2] == 2.+1.*I );

	auto m_as_operator2 = [&](auto const& B){return blas::gemm(m, B);};
	auto c2 = m_as_operator2(b);
	BOOST_REQUIRE( c == c2 );

	auto m_as_operator3 = [=](auto const& B){return blas::gemm(m, B);};
	auto c3 = m_as_operator3(b);
	BOOST_REQUIRE( c == c3 );

	multiplies_bind1st m_as_operator4(std::move(m));
	auto c4 = m_as_operator4(b);
	BOOST_REQUIRE( c == c4 );
	BOOST_REQUIRE( is_empty(m) );
}

BOOST_AUTO_TEST_CASE(multi_blas_gemm_elongated){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I}
	};
	BOOST_REQUIRE( size(a) == 1 and size(a[0]) == 2 );
	multi::array<complex, 2> const b = {
		{2. + 3.*I}, 
		{19.+11.*I}
	};
	BOOST_REQUIRE( size(b) == 2 and size(b[0]) == 1 );
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*b[0][0] + a[0][1]*b[1][0] );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., b, a, 0., c); // c=ba, c⸆=a⸆b⸆
		BOOST_REQUIRE( c[0][0] == b[0][0]*a[0][0] );
		BOOST_REQUIRE( c[1][1] == b[1][0]*a[0][1] );
		auto const c_copy = blas::gemm(1., b, a);
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, blas::H(a), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) + a[0][1]*conj(a[0][1]) );
		auto const c_copy = blas::gemm(1., a, blas::H(a));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> c({2, 2});
		blas::gemm(1., blas::H(a), a, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == a[0][0]*conj(a[0][0]) );
		auto const c_copy = blas::gemm(blas::H(a), a);
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> const a = {
			{1.-2.*I, 9.-1.*I}
		};
		multi::array<complex, 2> const b = {
			{2.+3.*I}, 
			{19.+11.*I}
		};
		multi::array<complex, 2> c({1, 1});
		blas::gemm(1., a, b, 0., c);
	}
	{
		multi::array<double, 2> const a = {{2., 3.}};
		multi::array<double, 2> const b = {{4., 5.}};
		multi::array<double, 2> c({1, 1});
		blas::gemm(1., a, blas::T(b), 0., c); // blas error
	}
	{
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4.},
			{5.},
			{6.}
		};
		multi::array<double, 2> c1({1, 1}), c2({1, 1});
		auto ra = rotated(a).decay();
		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

	//	gemm(1., rotated(a), b, 0., c2); // was blas error
	//	BOOST_REQUIRE(c1 == c2); // not reached
	}
	{
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4., 2.},
			{5., 1.},
			{6., 2.}
		};
		multi::array<double, 2> c1({1, 2}), c2({1, 2});
		auto ra = rotated(a).decay();

		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

		blas::gemm(1., blas::T(a), b, 0., c2);
		BOOST_REQUIRE(c1 == c2);
	}
	if(0){
		multi::array<double, 2> const a = {
			{2.},
			{3.},
			{5.}
		};
		multi::array<double, 2> const b = {
			{4.},
			{5.},
			{6.}
		};
		multi::array<double, 2> c1({1, 1}), c2({1, 1});
		auto ra = rotated(a).decay();
		blas::gemm(1., ra, b, 0., c1); // ok
		BOOST_REQUIRE( c1[0][0] == a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0] );

		blas::gemm(1., blas::T(a), b, 0., c2);
		BOOST_REQUIRE(c1 == c2);
	}
	if(0){
		multi::array<complex, 2> const a = {
			{2. + 1.*I},
			{3. + 2.*I}
		};
		multi::array<complex, 2> const b = {
			{4. + 3.*I}, 
			{5. + 4.*I}
		};
		multi::array<complex, 2> c1({1, 1}), c2({1, 1});
		auto ha = blas::hermitized(a).decay(); 
		blas::gemm(1., ha, b, 0., c1); // ok
		blas::gemm(1., blas::H(a), b, 0., c2); // was blas error
		print(c1);
		print(c2);
	//	BOOST_REQUIRE(c1 == c2);
	}
	{
		multi::array<complex, 2> const a = {
			{1.-2.*I}, 
			{9.-1.*I}
		};
		multi::array<complex, 2> const b = {
			{2.+3.*I, 2. + 999.*I}, 
			{19.+11.*I, 1. + 999.*I}
		};
		multi::array<complex, 2> c1({1, 2}), c2({2, 1}, 999.);
		auto ha = blas::hermitized(a).decay();
		blas::gemm(1., ha, b, 0., c1);
		blas::gemm(1., blas::H(b), a, 0., c2);
//		print(c1);
//		print(c2);
	//	std::cout << std::endl;
		BOOST_REQUIRE( c1 == blas::hermitized(c2) );
	}
	{
	//	multi::array<complex, 2> c({1, 1});
	//	using multi::blas::gemm;
	//	using multi::blas::hermitized;
//		gemm(1., hermitized(b), b, 0., c);
	//	BOOST_REQUIRE( c[0][0] == b[0][0]*conj(b[0][0]) + b[1][0]*conj(b[1][0]) );
	}
}

template<class A> void what(A&&) = delete;
template<class... A> void what() = delete;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic2){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{3.+3.*I, 7.-4.*I},
		{1.+9.*I, 1.+8.*I}
	};
	multi::array<complex, 2> const b = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );

		multi::array<complex, 2> const c_copy  = blas::gemm(1., blas::H(a), b);
		multi::array<complex, 2> const c_copy2 = blas::gemm(blas::H(a), b);
		BOOST_REQUIRE(( c == c_copy and c == c_copy2 ));
	}
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4}, acu.get_allocator());
		blas::gemm(1., blas::H(acu), bcu, 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );

	//	what(base(acu));
	//	what(blas::H(acu).get_allocator());
	//	what<decltype(blas::H(acu))::decay_type, decltype(acu)::decay_type>();
	//	cuda::array<complex, 2> const ccu_copy  = blas::gemm(1., blas::H(acu), bcu);
	//	cuda::array<complex, 2> const ccu_copy2 = blas::gemm(blas::H(acu), bcu);
	//	BOOST_REQUIRE(( ccu_copy == ccu and ccu_copy2 == ccu ));
	}
#if 0
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., blas::H(amcu), bmcu, 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );

	//	[](void*){}();

		cuda::managed::array<complex, 2> const cmcu_copy  = blas::gemm(1., blas::H(amcu), bmcu);
		cuda::managed::array<complex, 2> const cmcu_copy2 = blas::gemm(blas::H(amcu), bmcu);
		BOOST_REQUIRE(( cmcu_copy == cmcu and cmcu_copy2 == cmcu ));
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic3){
	namespace cuda = multi::cuda;
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{3.+3.*I, 7.-4.*I},
		{1.+9.*I, 1.+8.*I}
	};
	multi::array<complex, 2> const bH = {	
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
	multi::array<complex, 2> const b = multi::blas::hermitized(bH);
	{
		multi::array<complex, 2> c({2, 4});
		blas::gemm(1., blas::H(a), blas::H(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
	{
		cuda::array<complex, 2> const acu = a;
		cuda::array<complex, 2> const bcu = b;
		cuda::array<complex, 2> ccu({2, 4});
		blas::gemm(1., blas::H(acu), blas::H(bcu), 0., ccu);
		BOOST_REQUIRE( ccu[1][2] == complex(112, 12) );
	}
	{
		cuda::managed::array<complex, 2> const amcu = a;
		cuda::managed::array<complex, 2> const bmcu = b;
		cuda::managed::array<complex, 2> cmcu({2, 4});
		blas::gemm(1., blas::H(amcu), blas::H(bmcu), 0., cmcu);
		BOOST_REQUIRE( cmcu[1][2] == complex(112, 12) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic4){
	namespace blas = multi::blas;

	multi::array<complex, 2> c({12, 12});
	{
		multi::array<complex, 2> const a({12, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		using multi::blas::hermitized;
		using multi::blas::gemm;
		blas::gemm(1., a, blas::H(b), 0., c);
		BOOST_REQUIRE( real(c[0][0]) > 0);

		auto c_copy = blas::gemm(1., a, blas::H(b));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		multi::array<complex, 2> const a_block({24, 100}, 1.+2.*I);
		multi::array<complex, 2> const b({12, 100}, 1.+2.*I);
		multi::array<complex, 2> c2({12, 12});

		blas::gemm(1., a_block.strided(2), blas::H(b), 0., c2);

		BOOST_REQUIRE( real(c[0][0]) > 0);
		BOOST_REQUIRE( c == c2 );

		auto c2_copy = blas::gemm(1., a_block.strided(2), blas::H(b));
		BOOST_REQUIRE( c2_copy == c2 );
	}
}

//BOOST_AUTO_TEST_CASE(multi_blas_gemm_complex_issue68){
//	namespace cuda = multi::cuda;
//	namespace blas = multi::blas;
//	cuda::managed::array<complex, 2> const a = {
//		{1.-2.*I},
//		{3.+3.*I},
//		{1.+9.*I}
//	};
//	{
//		cuda::managed::array<complex, 2> c({1, 1});
//		blas::gemm(1., blas::H(a), a, 0., c);
//		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
//	}
//#if 0
//	{
//		auto c = blas::gemm(2., blas::H(a), a);
//		BOOST_REQUIRE( c[0][0] == 210. + 0.*I );
//	}
//	{
//		auto c = blas::gemm(blas::H(a), a);

//		BOOST_REQUIRE( c[0][0] == 105. + 0.*I );
//	}
//#endif
//}
#endif

#if 0
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
//	{
//		boost::timer::auto_cpu_timer t; // 4.516157s
//		for(auto i : extension(~B)) (~C2)[i] = blas::gemv(A, (~B)[i]);
//	}
//	{
//		boost::timer::auto_cpu_timer t; // 4.516157s
//		for(auto i : extension(A)) C2[i] = blas::gemv(~B, A[i]);
//	}
	{
		boost::timer::auto_cpu_timer t; // 32.705804s
		for(auto i:extension(A)) for(auto j:extension(~B)) C3[i][j] = blas::dot(A[i], (~B)[j]);
	}
	using namespace blas::operators;

	BOOST_REQUIRE( std::equal(
		begin(C), end(C), begin(C2), [](auto const& crow, auto const& c2row){
			return ((crow - c2row)^2) < 1e-13;
		}
	) );

	BOOST_REQUIRE( std::equal(
		begin(C), end(C), begin(C3), [](auto const& crow, auto const& c3row){
			return ((crow - c3row)^2) < 1e-13;
		}
	) );

}
#endif

