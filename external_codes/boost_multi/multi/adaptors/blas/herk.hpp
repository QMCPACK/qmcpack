#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework `pkg-config --libs blas` \
`#-Wl,-rpath,/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -L/usr/local/Wolfram/Mathematica/12.0/SystemFiles/Libraries/Linux-x86-64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5` \
-lboost_timer &&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_HERK_HPP
#define MULTI_ADAPTORS_BLAS_HERK_HPP

#include "../blas/core.hpp"
#include "../blas/copy.hpp" 
//#include "../blas/scal.hpp" 
#include "../blas/syrk.hpp" // fallback to real case

#include "../blas/side.hpp"
#include "../blas/filling.hpp"

#include "../blas/operations.hpp"

#include "../../config/NODISCARD.hpp"

//#include<iostream> //debug
//#include<type_traits> // void_t

namespace boost{
namespace multi{namespace blas{

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0> 
auto base_aux(A&& a)
->decltype(base(a)){
	return base(a);}

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0>
auto base_aux(A&& a)
->decltype(underlying(base(a))){
	return underlying(base(a));}

using core::herk;

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr, std::enable_if_t<is_complex_array<C2D>{}, int> =0>
C2D&& herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c)
//->decltype(herk('\0', '\0', c.size(), a.size(), &alpha, base_aux(a), stride(a.rotated()), &beta, base_aux(c), stride(c)), std::forward<C2D>(c))
{
	assert( a.size() == c.size() );
	assert( c.size() == rotated(c).size() );
	if(c.size()==0) return std::forward<C2D>(c);
	if constexpr(is_conjugated<C2D>{}){herk(flip(c_side), alpha, a, beta, hermitized(c)); return std::forward<C2D>(c);}
	{
		auto base_a = base_aux(a);
		auto base_c = base_aux(c); //  static_assert( not is_conjugated<C2D>{}, "!" );
		if constexpr(is_conjugated<A2D>{}){
		//	auto& ctxt = *blas::default_context_of(underlying(a.base()));
			// if you get an error here might be due to lack of inclusion of a header file with the backend appropriate for your type of iterator
				 if(stride(a)==1 and stride(c)!=1) herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(c));
			else if(stride(a)==1 and stride(c)==1){
				if(size(a)==1)                     herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(c));
				else assert(0);
			}
			else if(stride(a)!=1 and stride(c)==1) herk(c_side==filling::upper?'U':'L', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(rotated(c)));
			else if(stride(a)!=1 and stride(c)!=1) herk(c_side==filling::upper?'L':'U', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(        c ));
			else assert(0);
		}else{
		//	auto& ctxt = *blas::default_context_of(           a.base() );
			;;;; if(stride(a)!=1 and stride(c)!=1) herk(c_side==filling::upper?'L':'U', 'C', size(c), size(rotated(a)), &alpha, base_a, stride(        a ), &beta, base_c, stride(c));
			else if(stride(a)!=1 and stride(c)==1){
				if(size(a)==1)                     herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(rotated(c)));
				else assert(0);
			}
			else if(stride(a)==1 and stride(c)!=1) assert(0);//case not implemented, herk(c_side==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), alpha, base_a, stride(rotated(a)), beta, base(c), stride(c)); 
			else if(stride(a)==1 and stride(c)==1) herk(c_side==filling::upper?'U':'L', 'N', size(c), size(rotated(a)), &alpha, base_a, stride(rotated(a)), &beta, base_c, stride(rotated(c)));
			else assert(0);
		}
	}
	return std::forward<C2D>(c);
}

template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr, std::enable_if_t<not is_complex_array<C2D>{}, int> =0>
auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c)
->decltype(syrk(c_side, alpha, a, beta, std::forward<C2D>(c))){
	return syrk(c_side, alpha, a, beta, std::forward<C2D>(c));}

//template<class AA, class BB, class A2D, class C2D, class = typename A2D::element_ptr>
//auto herk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c)
//->decltype(herk_aux(c_side, alpha, a, beta, std::forward<C2D>(c), is_complex<C2D>{})){
//	return herk_aux(c_side, alpha, a, beta, std::forward<C2D>(c), is_complex<C2D>{});}

template<class AA, class A2D, class C2D, class = typename A2D::element_ptr>
auto herk(filling c_side, AA alpha, A2D const& a, C2D&& c)
->decltype(herk(c_side, alpha, a, 0., std::forward<C2D>(c))){
	return herk(c_side, alpha, a, 0., std::forward<C2D>(c));}

template<typename AA, class A2D, class C2D>
auto herk(AA alpha, A2D const& a, C2D&& c)
->decltype(herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(c)))){
	return herk(filling::lower, alpha, a, herk(filling::upper, alpha, a, std::forward<C2D>(c)));}

template<class A2D, class C2D>
auto herk(A2D const& a, C2D&& c)
->decltype(herk(1., a, std::forward<C2D>(c))){
	return herk(1., a, std::forward<C2D>(c));}

/*
template<class A2D, class C2D>
NODISCARD("when last argument is const")
auto herk(A2D const& a, C2D const& c)
->decltype(herk(1., a, decay(c))){
	return herk(1., a, decay(c));}
*/

template<class AA, class A2D, class Ret = typename A2D::decay_type>
NODISCARD("when argument is read-only")
auto herk(AA alpha, A2D const& a)//->std::decay_t<decltype(herk(alpha, a, Ret({size(a), size(a)}, get_allocator(a))))>{
{
	return herk(alpha, a, Ret({size(a), size(a)}));//Ret({size(a), size(a)}));//, get_allocator(a)));
}

template<class T> struct numeric_limits : std::numeric_limits<T>{};
template<class T> struct numeric_limits<std::complex<T>> : std::numeric_limits<std::complex<T>>{
	static std::complex<T> quiet_NaN(){auto n=numeric_limits<T>::quiet_NaN(); return {n, n};}
};

template<class AA, class A2D, class Ret = typename A2D::decay_type>
NODISCARD("because argument is read-only")
auto herk(filling cs, AA alpha, A2D const& a)
->std::decay_t<
decltype(herk(cs, alpha, a, Ret({size(a), size(a)}, 0., get_allocator(a))))>{
	return herk(cs, alpha, a, Ret({size(a), size(a)},
#ifdef NDEBUG
		numeric_limits<typename Ret::element_type>::quiet_NaN(),
#endif	
		get_allocator(a)
	));
}

template<class A2D> auto herk(filling s, A2D const& a)
->decltype(herk(s, 1., a)){
	return herk(s, 1., a);}

template<class A2D> auto herk(A2D const& a)
//->decltype(herk(1., a)){
{	return herk(1., a);}

}}

}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_HERK

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../adaptors/blas/gemm.hpp"
#include "../../adaptors/blas/nrm2.hpp"

#include<iostream>
#include<numeric>

namespace utf = boost::unit_test;
namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout << C[i][j] << ' ';
		cout << std::endl;
	}
	return cout << std::endl;
}

BOOST_AUTO_TEST_CASE(inq_case){
	using namespace multi::blas;
	multi::array<double, 2> const a = {
		{0, 1, 2},	
		{3, 4, 5},
		{6, 7, 8},
		{9, 10, 11}
	};
	BOOST_REQUIRE( gemm(a, T(a))[1][2] == 86. );
	{
		multi::array<double, 2> c({4, 4});
		herk(1.0, a, c);
		BOOST_REQUIRE( c == gemm(a, T(a)) );
	}
	{
		multi::array<double, 2> c = herk(1.0, a);
		BOOST_REQUIRE( c == gemm(a, T(a)) );
	}
	{
		BOOST_REQUIRE( herk(a) == gemm(a, T(a)) );
	}
	{
		BOOST_REQUIRE( herk(2.0, a) == gemm(2.0, a, T(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({2, 2}, 9999);
		blas::herk(1., a, c);
		BOOST_REQUIRE( c[1][0] == 34 );
		BOOST_REQUIRE( c[0][1] == 34 );

		multi::array<double, 2> const c_copy = blas::herk(1., a);
		BOOST_REQUIRE( c == c_copy );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case){
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case_scale){
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( B[0][0] == (1.*1. + 2.*2. + 3.*3.)*0.1 );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case){
	namespace blas = multi::blas;
	multi::array<complex, 2> const A = {{1., 2., 3.}};
	multi::array<complex, 2> B = blas::herk(1.0, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case_scale, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::array<complex, 2> const A = {{1., 2., 3.}};
	multi::array<complex, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( real( B[0][0]/0.1 ) == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case){
	namespace blas = multi::blas;
	multi::array<complex, 2> const A = {{1. + 2.*I, 2.+3.*I, 3. + 4.*I}};
	multi::array<complex, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(A)[0][0])) == blas::nrm2(A[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_out_param){
	namespace blas = multi::blas;
	multi::array<complex, 2> const A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	multi::array<complex, 2> B({1, 1});
	BOOST_REQUIRE( size(B) == 1 );

	blas::herk(blas::filling::upper, 1.0, blas::H(A), 0.0, B);

	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(B[0][0])) == blas::nrm2(blas::T(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized){
	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	namespace blas = multi::blas;
	multi::array<complex, 2> B = blas::herk(blas::H(A));
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_auto){
	namespace blas = multi::blas;

	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	auto B = blas::herk(1., blas::hermitized(A));
	static_assert( std::is_same<decltype(B), multi::array<complex, 2>>{}, "!" );
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

#if 1
#if 1

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_identity){
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};

	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::filling::lower, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		static_assert(blas::is_conjugated<decltype(blas::H(c))>{}, "!" );

		blas::herk(blas::filling::lower, 1., a, 0., blas::H(c)); // c†=c=aa†=(aa†)†, `c` in upper triangular

		BOOST_REQUIRE( blas::H(c)[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( blas::H(c)[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({2, 2}, 9999.);
	//	blas::herk(blas::filling::lower, 1., a, 0., blas::T(c)); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( transposed(c)[1][0]==complex(50., -49.) );
	//	BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aT(aT)† not supported
	//	print(c);
	//	BOOST_REQUIRE( c[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., hermitized(c)); // c†=c=aT(aT)† not supported
	//	BOOST_REQUIRE( hermitized(c)[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( hermitized(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(blas::filling::lower, 1., blas::T(a), 0., blas::T(c)); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( transposed(c)[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::filling::lower, 1., blas::T(a), 0., blas::H(blas::T(c))); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( blas::H(blas::T(c))[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( blas::H(blas::T(c))[0][1]==9999. );
	}
	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		using namespace multi::blas;
//		blas::herk(blas::filling::lower, 1., blas::T(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
//		BOOST_REQUIRE( c[0][1]==9999. );
	}
#if 1
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::U, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(1., a, c); // c†=c=aa†=(aa†)†
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::L, 1., blas::H(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({3, 3}, 9999.);
	//	using namespace multi::blas;
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( c[0][1]==9999. );
	//	BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
	}
#endif
}

#if 0

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_real_case){
	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::transposed;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);

		herk(filling::lower, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(19.,0.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][2]==complex(19.,0.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::upper, 1., hermitized(a), 0., transposed(c));//c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(transposed(c));
	//	BOOST_REQUIRE( c[1][2]==complex(19.,0.) );
	//	BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using blas::transposed;
	//	herk(filling::upper, 1., transposed(a), 0., c);//c_†=c_=a_†a_=(a_†a_)†, `c_` in lower triangular
	//	BOOST_REQUIRE( c[2][1] == 9999. );
	//	BOOST_REQUIRE( c[1][2] == 19. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_transparent_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, information in `c` lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using multi::blas::herk;
		herk(filling::upper, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†, `a` and `c` are c-ordering, information in c lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(filling::upper, 1., a, 0., c); //c†=c=aa†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_enum_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	using blas::transposed;
	{
	//	multi::array<complex, 2> c({2, 2}, 8888.);
	//	std::cerr << "here" << std::endl;
	//	herk(filling::lower, 1., hermitized(transposed(a)), 0., c); //c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(c) << std::endl;
	//	std::cerr << "there" << std::endl;
	//	BOOST_REQUIRE( c[0][1]==complex(41.,2.) );
	//	BOOST_REQUIRE( c[1][0]==8888. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); //c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using namespace multi::blas;
		herk(filling::upper, 1., hermitized(a), 0., c); //c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using namespace multi::blas;
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using namespace multi::blas;
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_explicit_enum_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	BOOST_REQUIRE( herk(hermitized(a)) == gemm(hermitized(a), a) );
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., hermitized(a), 0., transposed(c)); // c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(transposed(c));
	//	BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
	//	BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., hermitized(transposed(a)), 0., transposed(c)); // c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( transposed(c)[1][0]==complex(50.,+49.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
//	BOOST_REQUIRE( herk(hermitized(transposed(a))) == gemm(hermitized(transposed(a)), transposed(a)) );
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
		BOOST_REQUIRE( herk(hermitized(a)) == gemm(hermitized(a), a) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
		BOOST_REQUIRE( herk(a) == gemm(a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(a) == gemm(a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 2., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(100., 98.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(2., a) == gemm(2., a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(1., a) == gemm(1., a, hermitized(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_operator_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), 0., c); // c=c†=a†a, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi:: blas::filling;
		herk(filling::lower, 1., a, 0., c); // c=c†=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(1., a, c); // c=c†=aa†
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), 0., c); // c=c†=a†a, `c` in lower triangular
		herk(filling::upper, 1., hermitized(a), 0., c);
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_operator_interface_implicit_no_sum){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), c); // c=c†=a†a, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::filling;
		herk(filling::lower, 1., a, c); // c=c†=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_ordering_and_symmetrization){

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::herk;
	using blas::hermitized;
	using blas::filling;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==9999. );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(1., hermitized(a), c); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==complex(41., +2.) );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa† // c implicit hermitic in upper
		BOOST_REQUIRE( c[1][0] == 9999. );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(1., a, c); // c†=c=aa†
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(1., a); // c†=c=aa†
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, 1., hermitized(a)); // c†=c=a†a

		BOOST_REQUIRE( size(hermitized(a))==3 );		
	//	BOOST_REQUIRE( c[2][1] == complex(41., +2.) );
		BOOST_REQUIRE( c[1][2] == complex(41., -2.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, a); // c†=c=a†a
//		what(multi::pointer_traits<decltype(base(a))>::default_allocator_of(base(a)));
	//	BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, hermitized(a)); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1] == complex(41., +2.) );
		BOOST_REQUIRE( c[1][2] == complex(41., -2.) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size1_real_case){
	multi::array<complex, 2> const a = {
		{1., 3., 4.}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({1, 1}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
		BOOST_TEST( c[0][0] == 26. );
	}
	BOOST_TEST( herk(a) == gemm(a, hermitized(a)) );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size1){
	multi::array<complex, 2> const a = {
		{1. + 4.*I, 3. + 2.*I, 4. - 1.*I}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({1, 1}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
		BOOST_TEST( c[0][0] == 47. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size0){
	multi::array<complex, 2> const a;
	using namespace multi::blas;
	{
		multi::array<complex, 2> c;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_TEST( c[0][0] == 47. );
	}
}


BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_ordering_and_symmetrization_real_case){

	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., hermitized(a)); // c†=c=a†a
		BOOST_REQUIRE( size(hermitized(a))==3 );		
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, a); // c†=c=a†a
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, hermitized(a)); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
}


BOOST_AUTO_TEST_CASE(multi_blas_herk_real_automatic_ordering_and_symmetrization_real_case){

	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({3, 3}, 9999.);
		using multi::blas::hermitized;
		using multi::blas::herk;
		using multi::blas::filling;
	//	herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
	//	BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		using multi::blas::filling;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		using multi::blas::filling;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<double, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		multi::array<complex, 2> c = herk(a); // c†=c=a†a
		BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		multi::array<complex, 2> c = herk(hermitized(a)); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real_case){
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	using multi::blas::filling;
	{
		static_assert( not boost::multi::blas::is_complex_array<multi::array<double, 2>>{} , "!");
		multi::array<double, 2> c({2, 2}, 9999.);
		syrk(filling::lower, 1., a, 0., c);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., a, 0., c);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
	{
		static_assert( not boost::multi::blas::is_complex_array<multi::array<double, 2>>{} , "!");
		multi::array<double, 2> c = herk(filling::upper, a);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_real_case_1d){
	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::transposed;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		print(c);
		BOOST_REQUIRE( c[2][1]==complex(12.,0.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(2., hermitized(a), c);//c†=c=a†a=(a†a)†, `c` in lower triangular

		BOOST_REQUIRE( c[2][1]==complex(24.,0.) );
		BOOST_REQUIRE( c[1][2]==complex(24.,0.) );
		multi::array<complex, 2> c_gemm({3, 3});
	//	gemm(2., hermitized(a), a, c_gemm);
	}
}
#endif

#if 0
BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_timing){
	multi::array<complex, 2> const a({4000, 4000}); std::iota(data_elements(a), data_elements(a) + num_elements(a), 0.2);
	multi::array<complex, 2> c({4000, 4000}, 9999.);
	boost::timer::auto_cpu_timer t;
	using multi::blas::herk;
	using multi::blas::hermitized;
	using multi::blas::filling;
	herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
}
#endif
#endif
#endif

#endif
#endif

