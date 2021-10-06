#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2021

#ifndef MULTI_ADAPTORS_BLAS_OPERATIONS_HPP
#define MULTI_ADAPTORS_BLAS_OPERATIONS_HPP

#include "../blas/numeric.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class M> auto transposed(M const& m) -> decltype(auto){return rotated(m);}

template<class A, typename D=std::decay_t<A>, typename E=typename D::element_type>
auto conjugated_transposed(A&& a) -> decltype(auto){
	return transposed(blas::conj(std::forward<A>(a)));
}

template<class A> auto identity(A&& a) -> decltype(auto){return std::forward<A>(a);}

template<class A>
auto hermitized(A&& a, std::true_type /*true */) -> decltype(auto){
	return conjugated_transposed(std::forward<A>(a));
}

template<class A>
auto hermitized(A&& a, std::false_type /*false*/) -> decltype(auto){
	return transposed(std::forward<A>(a));
}

template<class A>
auto hermitized(A&& a) -> decltype(auto){return conjugated_transposed(std::forward<A>(a));}

template<class A>
auto transposed(A&& a) -> decltype(auto){return rotated(std::forward<A>(a));}

namespace operators{

MAYBE_UNUSED constexpr static struct {

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0>
	auto operator()(A&& a) const -> decltype(auto){return hermitized(std::forward<A>(a));}

	template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0>
	[[deprecated("use blas::C instead of blas::H for conjugated vectors to avoid confusions")]]
	auto operator()(A&& a) const -> decltype(auto){return blas::conj(std::forward<A>(a));}

} H;

template<class A, class Op>
auto operator^(A&& a, Op op)
->decltype(op(std::forward<A>(a))){
	return op(std::forward<A>(a));}

} // end namespace operators

using operators::H;

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 1, int> =0> 
auto C(A&& a) -> decltype(auto){return blas::conj(std::forward<A>(a));}

template<class A, std::enable_if_t<std::decay_t<A>::rank_v == 2, int> =0> 
auto C(A&& a) -> decltype(auto){return hermitized(std::forward<A>(a));}

namespace operators{

	template<class A>
	auto operator*(A&& a)
	->decltype(blas::conj(std::forward<A>(a))){
		return blas::conj(std::forward<A>(a));}

} // end namespace operators

template<class A> auto T(A&& a) -> decltype(auto){return transposed(std::forward<A>(a));}
template<class A> auto N(A&& a) -> decltype(auto){return identity  (std::forward<A>(a));}

} // end namespace blas
} // end namespace multi
} // end namespace boost

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

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

BOOST_AUTO_TEST_CASE(m){
	using complex = std::complex<double>; constexpr complex I{0., 1.};

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

	using namespace blas::operators;
	BOOST_REQUIRE( (*~A)[0][1] == conj(A[1][0]) );
	BOOST_REQUIRE( (~*A)[0][1] == conj(A[1][0]) );
	BOOST_REQUIRE( ( ~A)[0][1] ==      A[1][0]  );
	BOOST_REQUIRE( ( *A)[0][1] == conj(A[0][1]) );

}

BOOST_AUTO_TEST_CASE(is_complex_array_test){
	static_assert(multi::blas::is_complex_array<multi::array<std::complex<double>, 2>>{}, "!");
}

#endif
#endif

