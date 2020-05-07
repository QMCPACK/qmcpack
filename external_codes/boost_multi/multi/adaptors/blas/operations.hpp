#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include    "../blas/numeric.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class M> decltype(auto) transposed(M const& m){return rotated(m);}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
decltype(auto) conjugated_transposed(A&& a){
	return transposed(blas::conj(std::forward<A>(a)));
}

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
decltype(auto) hermitized(A&& a){return conjugated_transposed(std::forward<A>(a));}

template<class A>
decltype(auto) transposed(A&& a){return rotated(std::forward<A>(a));}

template<class A> [[deprecated("use blas::H instead of blas::C for hermitized to avoid confusions")]]
decltype(auto) C(A&& a){return hermitized(std::forward<A>(a));}
template<class A> decltype(auto) H(A&& a){return hermitized(std::forward<A>(a));}
template<class A> decltype(auto) T(A&& a){return transposed(std::forward<A>(a));}
template<class A> decltype(auto) N(A&& a){return identity  (std::forward<A>(a));}

}}

}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_OPERATIONS

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi blas operations"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"

using std::cout;
template<class M> decltype(auto) print(M const& C){
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<"---"<<std::endl;
}

namespace multi = boost::multi;
using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(m){
	namespace blas = multi::blas;
	multi::array<complex, 2> const A = {
		{1. - 3.*I, 6.  + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};
	using blas::hermitized;
	BOOST_REQUIRE( hermitized(A)[0][1] == conj(A[1][0]) );

	static_assert( blas::is_conjugated<decltype(blas::H(A))>{}, "!" );
	BOOST_REQUIRE( blas::H(A)[0][1] == conj(A[1][0]) );

	using blas::transposed;
	BOOST_REQUIRE( transposed(A)[0][1] == A[1][0] );

	static_assert( not blas::is_conjugated<decltype(blas::T(A))>{}, "!" );
	BOOST_REQUIRE( blas::T(A)[0][1] == A[1][0] );
	
//	static_assert( multi::blas::is_conjugated<decltype(T(A))>{}, "!" );

/*	using multi::blas::gemm;



	BOOST_REQUIRE( gemm(A, hermitized(A))[2][1] == 20. - 14.*I );


	BOOST_REQUIRE( gemm(A, transposed(A))[2][1] == 16. + 2.*I );

	static_assert( multi::blas::is_conjugated_t<decltype(hermitized(A))>{} , "!" );
	static_assert( not multi::blas::is_conjugated_t<std::decay_t<decltype( conjugated(hermitized(A)) )>>{}, "!");
	static_assert( not multi::blas::is_hermitized<std::decay_t<decltype( conjugated(hermitized(A)) )>>{}, "!");
*/
}

BOOST_AUTO_TEST_CASE(is_complex_array_test){
	static_assert(multi::blas::is_complex<multi::array<std::complex<double>, 2>>{}, "!");
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

