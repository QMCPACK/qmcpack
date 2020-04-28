#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -D_TEST_MULTI_ADAPTORS_BLAS_GEMM $0.cpp -o $0x -lboost_unit_test_framework \
`pkg-config --libs blas` \
`#-Wl,-rpath,/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -L/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core` \
&&$0x&& rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2019-2020
#ifndef MULTI_ADAPTORS_BLAS_GEMM_HPP
#define MULTI_ADAPTORS_BLAS_GEMM_HPP

#include "../blas/core.hpp"
#include "../blas/operations.hpp"
#include "../blas/numeric.hpp"

#include "../../config/NODISCARD.hpp"

namespace boost{
namespace multi{
namespace blas{

using multi::blas::core::gemm;

template<class M> bool is_c_ordering(M const& m){
	return stride(rotated(m))==1 and size(m)!=1;
}

double const& conj(double const& d){return d;}
float const& conj(float const& f){return f;}

template<class A> auto gemm_base_aux(A&& a, std::false_type){return base(a);}
template<class A> auto gemm_base_aux(A&& a, std::true_type){return underlying(base(a));}

template<class A2D, class B2D, class C2D>
C2D&& gemm(typename std::decay_t<C2D>::element_type alpha, A2D const& a, B2D const& b, typename std::decay_t<C2D>::element_type beta, C2D&& c){
	assert( size(rotated(a)) == size(b) );
	assert( size(c) == size(a) );
	assert( size(rotated(b)) == size(rotated(c)) );
	auto base_a = gemm_base_aux(a, is_hermitized<A2D>{}); (void)base_a;
	auto base_b = gemm_base_aux(b, is_hermitized<B2D>{}); (void)base_b;
	auto base_c = gemm_base_aux(c, is_hermitized<std::decay_t<C2D>>{}); (void)base_c;

	if(is_conjugated(c)) blas::gemm(conj(alpha), conjugated(a), conjugated(b), conj(beta), conjugated(c));
	else{
		if(is_c_ordering(c)){//gemm(alpha, transposed(b), transposed(a), beta, transposed(c));
			     if( is_c_ordering(a) and  is_c_ordering(b)){
				assert(!is_conjugated(a) and !is_conjugated(b));
				gemm('N', 'N', size(rotated(c)), size(a   ), size(b   ), &alpha, base_b, stride(b   ), base_a, stride(a   ), &beta, base_c, stride(c   ));}
			else if(!is_c_ordering(a) and  is_c_ordering(b)){
				assert(!is_conjugated(b));
				gemm('N', is_conjugated(a)?'C':'T', size(rotated(c)), size(a   ), size(b   ), &alpha, base_b, stride(b   ), base_a, stride(rotated(a)), &beta, base_c, stride(c   ));
			}
			else if( is_c_ordering(a) and !is_c_ordering(b)){
				assert(!is_conjugated(a));
				gemm(is_conjugated(b)?'C':'T', 'N', size(rotated(c)), size(a   ), size(b   ), &alpha, base_b, stride(rotated(b)), base_a, stride(a   ), &beta, base_c, stride(c   ));
			}else if(!is_c_ordering(a) and !is_c_ordering(b)){gemm(is_conjugated(b)?'C':'T', is_conjugated(a)?'C':'T', size(rotated(c)), size(a   ), size(b   ), &alpha, base_b, stride(rotated(b)), base_a, stride(rotated(a)), &beta, base_c, stride(c   ));}
		}else{
			using core::gemm;
				 if( is_c_ordering(a) and  is_c_ordering(b)){
				gemm(is_conjugated(a)?'C':'T', is_conjugated(b)?'C':'T', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(a   ), base_b, stride(b   ), &beta, base_c, stride(rotated(c)));
			}else if(!is_c_ordering(a) and  is_c_ordering(b)){
				if(is_conjugated(a) and size(a)==1){
					gemm('C', is_conjugated(b)?'C':'T', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, size(rotated(a)), base_b, stride(b   ), &beta, base_c, stride(rotated(c)));
				}else{
					assert(not is_conjugated(a));
					gemm('N', is_conjugated(b)?'C':'T', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(rotated(a)), base_b, stride(b   ), &beta, base_c, stride(rotated(c)));
				}
			}
			else if( is_c_ordering(a) and !is_c_ordering(b)){
				if(is_conjugated(b) and size(rotated(b))==1){
					gemm(is_conjugated(a)?'C':'T', 'C', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(a   ), base_b, size(rotated(b)), &beta, base_c, stride(rotated(c)));
				}else{
					assert(not is_conjugated(b));
					gemm(is_conjugated(a)?'C':'T', 'N', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(a   ), base_b, stride(rotated(b)), &beta, base_c, stride(rotated(c)));
				}
			}else if(!is_c_ordering(a) and !is_c_ordering(b)){
				      if(not is_conjugated(a) and is_conjugated(b) and size(rotated(b))==1){
					gemm('N', 'C', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(rotated(a)), base_b, stride(b), &beta, base_c, stride(rotated(c)));
				}else if(is_conjugated(a) and size(a)==1 and is_conjugated(b) and size(rotated(b))==1){
					gemm('C', 'C', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, size(rotated(a)), base_b, stride(b), &beta, base_c, stride(rotated(c)));
				}else{
					assert(not is_conjugated(a)); assert(not is_conjugated(b));
					gemm('N', 'N', size(c   ), size(rotated(b)), size(rotated(a)), &alpha, base_a, stride(rotated(a)), base_b, stride(rotated(b)), &beta, base_c, stride(rotated(c)));
				}
			}
		}
	}
	return std::forward<C2D>(c);
}

template<class AA, class A2D, class B2D, class C2D = typename A2D::decay_type>
NODISCARD("second argument is const")
auto gemm(AA a, A2D const& A, B2D const& B){
	assert(get_allocator(A) == get_allocator(B));
	return gemm(a, A, B, 0., C2D({size(A), size(rotated(B))}, get_allocator(A)));
}

template<class A2D, class B2D> auto gemm(A2D const& A, B2D const& B){return gemm(1., A, B);}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_GEMM

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include <boost/timer/timer.hpp>

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

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

BOOST_AUTO_TEST_CASE(multi_blas_gemm_nh){
	using complex = std::complex<double>; complex const I{0,1};
	multi::array<complex, 2> const a = {
		{1.-2.*I, 9.-1.*I},
		{2.+3.*I, 1.-2.*I}
	};
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(a), 0., c); // c=aa†, c†=aa†
		BOOST_REQUIRE( c[1][0] == 7.-10.*I );
		BOOST_REQUIRE( c[0][1] == 7.+10.*I );
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
		print(c);
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

		BOOST_REQUIRE( is_conjugated(hermitized(a)) );
		BOOST_REQUIRE( is_conjugated(hermitized(b)) );
		BOOST_REQUIRE( !is_c_ordering(hermitized(a)) );
		BOOST_REQUIRE( !is_c_ordering(hermitized(b)) );
		BOOST_REQUIRE( size(hermitized(a)) == 1 );
		BOOST_REQUIRE( size(hermitized(b)[0]) == 1 );

		gemm(1., hermitized(a), hermitized(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.+7.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_square){
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 148 );
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 169 and c[1][0] == 82 ));
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., a, transposed(b), 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 183 );		
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), transposed(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == 117 );		
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), transposed(b), 0., transposed(c)); // c⸆=a⸆b⸆, c=ba
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_empty){
	multi::array<double, 2> const a({0, 5});
	BOOST_REQUIRE( size(a) == 0 );
	BOOST_REQUIRE( size(rotated(a)) == 5 );
	BOOST_REQUIRE( a.empty() );
//	assert( rotated(a).empty() );
	multi::array<double, 2> const b;
	{
		multi::array<double, 2> c;//({3, 2});
		using multi::blas::gemm;
//		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
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
		gemm(1., a, b, 0., rotated(c)); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 31 );		
	}
	{
		auto ar = rotated(a).decay();
		multi::array<double, 2> c({3, 2});
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., c); // c=ab, c⸆=b⸆a⸆
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
		gemm(1., rotated(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 101 );
		gemm(1., rotated(a), b, 0., rotated(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 101 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_1x3_3x2){
	multi::array<double, 2> const a = {
		{1, 9, 1}
	};
	multi::array<double, 2> const b = {	
		{ 11, 12},
		{  7, 19},
		{  8, 1 }
	};
	{
		multi::array<double, 2> c({1, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[0][1] == 184 );
	}
	{
		multi::array<double, 2> c({2, 1});
		auto ar = rotated(a).decay();
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., rotated(c)); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE( c[1][0] == 184 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_2x3_3x1){
	multi::array<double, 2> const a = {
		{1, 9, 1},
		{3, 4, 5}
	};
	multi::array<double, 2> const b = {	
		{ 11},
		{  7},
		{  8}
	};
	using multi::blas::is_conjugated;
	BOOST_REQUIRE( not is_conjugated(a) );
	{
		multi::array<double, 2> c({1, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., rotated(c)); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( rotated(c)[1][0] == 101 );
		BOOST_REQUIRE( c[0][1] == 101 );
	}
	{
		multi::array<double, 2> c({2, 1});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c⸆=ab, c=b⸆a⸆
		BOOST_REQUIRE( rotated(c)[0][1] == 101 );
		BOOST_REQUIRE( c[1][0] == 101 );
	}
	{
		multi::array<double, 2> c({1, 2});
		auto ar = rotated(a).decay();
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., rotated(c)); // c⸆=ab, c⸆=b⸆a⸆
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto ar = rotated(a).decay();
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., c); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = rotated(b).decay();
		using multi::blas::gemm;
		gemm(1., a, rotated(br), 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
	{
		multi::array<double, 2> c({1, 1});
		auto br = rotated(b).decay();
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(br), 0., c);
		BOOST_REQUIRE( c[0][0] == 82 );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == 170.-8.*I and c[1][0] == 77.+42.*I ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., a, transposed(b), 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == 177.+69.*I );		
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::gemm;
		gemm(1., transposed(a), transposed(b), 0., c); // c=a⸆b⸆, c⸆=ba
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
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto ar = rotated(a).decay();
		using multi::blas::gemm;
		gemm(1., rotated(ar), b, 0., c); // c=ab, c⸆=ba
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = rotated(b).decay();
		using multi::blas::gemm;
		gemm(1., a, rotated(br), 0., c);
		BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		auto br = rotated(b).decay();
		using multi::blas::gemm;
		using multi::blas::hermitized;
	//	gemm(1., a, hermitized(br), 0., rotated(c));
	//	BOOST_REQUIRE( c[0][0] == 84.-7.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_hermitized_square){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::gemm;
		using multi::blas::conjugated;
		gemm(1., a, b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::conjugated;
		gemm(1., conjugated(a), conjugated(b), 0., conjugated(c)); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 145. + 43.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == 87. - 16.*I );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., a, hermitized(b), 0., c); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == 189. - 23.*I );		
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), hermitized(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == 109. - 68.*I);		
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
	//	gemm(1., hermitized(a), hermitized(b), 0., rotated(c)); // c⸆=a†b†, c=b*a*
	//	print(rotated(c));
	//	BOOST_REQUIRE( c[0][1] == 109. - 68.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
	}
	{
		multi::array<complex, 2> c({1, 1});
		using multi::blas::hermitized;
		auto ha = hermitized(a).decay();
		using multi::blas::gemm;
		gemm(1., ha, b, 0., c);
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );
		gemm(1., hermitized(b), a, 0., c);
		BOOST_REQUIRE( c[0][0] == 80.+53.*I );
	}
}
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_1x3_3x2){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 20.+21.*I );
	}
	{
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({1, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(ar), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );		
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x2){
	using complex = std::complex<double>; complex const I{0,1};
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
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({1, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][1] == 28.+3.*I );		
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x2){
	using complex = std::complex<double>; complex const I{0,1};
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
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({2, 2});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );		
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x2_3x1){
	using complex = std::complex<double>; complex const I{0,1};
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
		auto ar = rotated(a).decay();
		multi::array<complex, 2> c({2, 1});
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 125.-84.*I );		
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_3x1_3x1_bis){
	using complex = std::complex<double>; complex const I{0,1};
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
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., hermitized(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[0][0] == 80.-53.*I );		
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == 148 and c[1][1] == 241 );
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., a, rotated(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][1] == 196. );
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == 169. and c[1][0] == 82. ));
	}
	{
		multi::array<double, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(a), rotated(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][1] == 154. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_real_nonsquare_automatic){
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == 17. );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_square_automatic){
	using complex = std::complex<double>; constexpr complex I{0,1};
	multi::array<complex, 2> const a = {
		{ 1. + 2.*I, 3. - 3.*I},
		{ 9. + 1.*I, 7. + 4.*I},
	};
	multi::array<complex, 2> const b = {	
		{ 11. + 1.*I, 12. + 1.*I},
		{  7. + 8.*I, 19. - 2.*I},
	};
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., a, rotated(b), 0., c); // c=ab⸆, c⸆=ba⸆
		BOOST_REQUIRE( c[1][0] == complex(178, 75) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(a), b, 0., c); // c=a⸆b, c⸆=b⸆a
		BOOST_REQUIRE(( c[1][1] == complex(180, 29) and c[1][0] == complex(53, 54) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(a), rotated(b), 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE(( c[1][1] == complex(186, 65) and c[1][0] == complex(116, 25) ));
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][0] == complex(115, 104) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), b, 0., c); // c=a†b, c†=b†a
		BOOST_REQUIRE( c[1][0] == complex(111, 64) and c[1][1] == complex(158, -51) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., a, hermitized(b), 0., c); // c=ab†, c†=ba†
		BOOST_REQUIRE( c[1][0] == complex(188, 43) and c[1][1] == complex(196, 25) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::hermitized;
		gemm(1., hermitized(a), hermitized(b), 0., c); // c=a†b†, c†=ba
		BOOST_REQUIRE( c[1][0] == complex(116, -25) and c[1][1] == complex(186, -65) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::transposed;
		using multi::blas::hermitized;
		using multi::blas::gemm;
		gemm(1., transposed(a), hermitized(b), 0., c); // c=a⸆b†, c†=ba⸆†
		BOOST_REQUIRE( c[1][0] == complex(118, 5) and c[1][1] == complex(122, 45) );
	}
	{
		multi::array<complex, 2> c({2, 2});
		using multi::blas::gemm;
		using multi::blas::transposed;
		gemm(1., transposed(a), transposed(b), 0., c); // c=a⸆b⸆, c⸆=ba
		BOOST_REQUIRE( c[1][0] == complex(116, 25) and c[1][1] == complex(186, 65) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_gemm_complex_nonsquare_automatic){
	using complex = std::complex<double>; constexpr complex I{0,1};
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
		using multi::blas::gemm;
		gemm(1., a, b, 0., c); // c=ab, c⸆=b⸆a⸆
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
	}
}

#endif
#endif

