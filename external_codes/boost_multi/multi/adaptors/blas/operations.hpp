#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&nvcc -x cu --expt-relaxed-constexpr`#c++ -Wall -Wextra -Wpedantic` -D_TEST_MULTI_ADAPTORS_BLAS_OPERATIONS $0.cpp -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019

#ifndef MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include    "../blas/core.hpp"
#include    "../blas/asum.hpp"
#include    "../blas/numeric.hpp"

#include "../../array_ref.hpp"

//#include<experimental/functional> // std::identity

namespace boost{
namespace multi{namespace blas{

template<class M> decltype(auto) transposed(M const& m){return rotated(m);}
//template<class M> decltype(auto) transposed(M&       m){return rotated(m);}

template<class T, typename = decltype(std::declval<typename T::element>().imag())>
std::true_type is_complex_array_aux(T const&);
std::false_type is_complex_array_aux(...);

template <typename T> struct is_complex_array: decltype(is_complex_array_aux(std::declval<T const&>())){};

template<class ComplexPtr> std::true_type is_conjugated_aux(multi::blas::detail::conjugater<ComplexPtr> const&);
template<class T> std::false_type is_conjugated_aux(T const&);

template<class A>
struct is_conjugated_t : decltype(is_conjugated_aux(typename std::decay_t<A>::element_ptr{})){};

template<class A> constexpr bool is_conjugated(A const&){return is_conjugated_t<A>{};}

template<class A> constexpr bool is_not_conjugated(A const& a){return is_complex_array<A>{} and not is_conjugated(a);}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element, class C=detail::conjugater<typename D::element_ptr>, typename = std::enable_if_t<not is_conjugated_t<std::decay_t<A>>()> >
decltype(auto) conjugated(A&& a, void* = 0){
	return multi::static_array_cast<E, C>(std::forward<A>(a));
}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element, class C=typename D::element_ptr::underlying_type, typename = std::enable_if_t<is_conjugated_t<std::decay_t<A>>{}> >
decltype(auto) conjugated(A&& a){
	return multi::static_array_cast<E, C>(std::forward<A>(a));
}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element>
decltype(auto) conjugated_transposed(A&& a){
	return transposed(conjugated(a));
}

template<class ComplexPtr> std::true_type is_hermitized_aux(multi::blas::detail::conjugater<ComplexPtr> const&);
template<class T> std::false_type is_hermitized_aux(T const&);

template<class A>
struct is_hermitized : std::decay_t<decltype(is_hermitized_aux(typename std::decay_t<A>::element_ptr{}))>{};

template<class A> decltype(auto) identity(A&& a){return std::forward<A>(a);}


template<class A>
decltype(auto) hermitized(A&& a, std::true_type){
	return conjugated_transposed(std::forward<A>(a));
}

template<class A>
decltype(auto) hermitized(A&& a, std::false_type){
	return transposed(std::forward<A>(a));
}

template<class A>
decltype(auto) hermitized(A&& a){
#if __cpp_if_constexpr>=201606
	if constexpr(is_complex_array<std::decay_t<A>>{}){
		return conjugated_transposed(std::forward<A>(a));
	}else{
		return transposed(std::forward<A>(a));
	}
#else
	return hermitized(std::forward<A>(a), is_complex_array<std::decay_t<A>>{});
#endif
}

template<class A>
decltype(auto) transposed(A&& a){return rotated(std::forward<A>(a));}

}}

}

#if _TEST_MULTI_ADAPTORS_BLAS_OPERATIONS

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"
#include "../blas/nrm2.hpp"
#include "../blas/gemm.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;

template<class M> 
decltype(auto) print(M const& C){
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<"---"<<std::endl;
}

namespace multi = boost::multi;
using complex = std::complex<double>;
auto const I = complex(0., 1.);

template<class T> void what();

BOOST_AUTO_TEST_CASE(m){
	multi::array<complex, 2> const A = {
		{1. - 3.*I, 6.  + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	using multi::blas::gemm;

	using multi::blas::hermitized;
	BOOST_REQUIRE( hermitized(A)[0][1] == conj(A[1][0]) );
	BOOST_REQUIRE( gemm(A, hermitized(A))[2][1] == 20. - 14.*I );

	using multi::blas::transposed;
	BOOST_REQUIRE( transposed(A)[0][1] == A[1][0] );
	BOOST_REQUIRE( gemm(A, transposed(A))[2][1] == 16. + 2.*I );

	static_assert( multi::blas::is_conjugated_t<decltype(hermitized(A))>{} , "!" );
	static_assert( not multi::blas::is_conjugated_t<std::decay_t<decltype( conjugated(hermitized(A)) )>>{}, "!");
	static_assert( not multi::blas::is_hermitized<std::decay_t<decltype( conjugated(hermitized(A)) )>>{}, "!");
}

BOOST_AUTO_TEST_CASE(is_complex_array_test){
	static_assert(multi::blas::is_complex_array<multi::array<std::complex<double>, 2>>{}, "!");
}

#if 0
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_operations_enums){
	BOOST_REQUIRE( multi::blas::operation::identity == multi::blas::real_operation::identity );
	BOOST_REQUIRE( multi::blas::operation::transposition == multi::blas::real_operation::transposition );
	BOOST_REQUIRE( multi::blas::operation::hermitian == multi::blas::complex_operation::hermitian );
	BOOST_REQUIRE( multi::blas::operation::identity == multi::blas::complex_operation::identity );

	BOOST_REQUIRE( multi::blas::operation{multi::blas::real_operation::identity} == multi::blas::real_operation::identity );
	BOOST_REQUIRE( multi::blas::operation{multi::blas::real_operation::transposition} == multi::blas::real_operation::transposition );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_operations){

	multi::array<complex, 2> const A = {
		{1. - 3.*I, 6.  + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	print(A);
	print(multi::blas::conjugated(A));

	auto&& Aconjd = multi::blas::conjugated(A);
	assert( Aconjd[1][2] == conj(A[1][2]) );
	multi::array<complex, 2> Aconj = multi::blas::conjugated(A);
	assert( Aconj[1][2] == conj(A[1][2]) );
	assert( Aconjd == Aconj );

	auto&& Aconjdconjd = multi::blas::conjugated(Aconjd);
	assert( Aconjdconjd[1][2] == A[1][2] );
	assert( &Aconjdconjd[1][2] == &A[1][2] );

	auto&& Atranspd = multi::blas::transposed(A);
	assert( Atranspd[1][2] == A[2][1] );
	multi::array<complex, 2> Atransp = multi::blas::transposed(A);
	assert( Atransp[1][2] == A[2][1] );
	assert( Atransp == Atranspd );

	auto&& Aconjdtranspd = multi::blas::conjugated_transposed(A); (void)Aconjdtranspd;
	assert( Aconjdtranspd[1][2] == conj(A[2][1]) );
	auto Aconjtransp = multi::blas::conjugated_transposed(A).decay();
	
	assert( Aconjtransp[1][2] == conj(A[2][1]) );
	assert( Aconjdtranspd == Aconjtransp );

	
{
	multi::array<complex, 2> const A = {
		{1. - 3.*I, 6.  + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	using multi::blas::hermitized;
	assert( hermitized(A)[0][1] == conj(A[1][0]) );
//	[]{}(hermitized(A));
	static_assert( multi::blas::is_conjugated<decltype(hermitized(A))>{} , "!");

	using multi::blas::conjugated;
//	[]{}(conjugated(conjugated(A)));

	using multi::blas::hermitized;
	[]{}(hermitized(hermitized(A)));

//	static_assert( not multi::blas::is_conjugated<decltype(hermitized(hermitized(A)))>{} , "!");

//	[]{}(hermitized(hermitized(A)));
//	[]{}(conjugated(conjugated(A)));

	static_assert( multi::blas::is_complex_array<std::decay_t<decltype(A)>>{} , "!");
//	auto&& AH = multi::blas::hermitized(A);
//	auto c = AH[0][0].imag();
//	static_assert( multi::blas::is_complex_array<std::decay_t<decltype(AH)>>{} , "!");

//	auto&& Aconjd = multi::blas::conjugated(A);
//	assert( Aconjd[1][2] == conj(A[1][2]) );
//	multi::array<complex, 2> Aconj = multi::blas::conjugated(A);
//	assert( Aconj[1][2] == conj(A[1][2]) );
//	assert( Aconjd == Aconj );

	auto&& Atranspd = multi::blas::T(A);
	assert( Atranspd[1][2] == A[2][1] );
	multi::array<complex, 2> Atransp = multi::blas::transposed(A);
	assert( Atransp[1][2] == A[2][1] );
	assert( Atransp == Atranspd );

	auto&& Aconjdtranspd = multi::blas::C(A); (void)Aconjdtranspd;
	assert( Aconjdtranspd[1][2] == conj(A[2][1]) );
	multi::array<complex, 2> Aconjtransp = multi::blas::conjugated_transposed(A);
	assert( Aconjtransp[1][2] == conj(A[2][1]) );
	assert( Aconjdtranspd == Aconjtransp );

}
	
}
#endif
#endif
#endif

