#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -D_TEST_MULTI_ADAPTORS_LAPACK_POTRF $0.cpp -o$0x `pkg-config --libs blas lapack` -lboost_unit_test_framework&&valgrind $0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_LAPACK_POTRF_HPP
#define MULTI_ADAPTORS_LAPACK_POTRF_HPP

#include "../../array.hpp"
#include "../../config/NODISCARD.hpp"

#include "../lapack/core.hpp"
#include "../blas/numeric.hpp"

#include "../blas/filling.hpp"

#include<cassert>

namespace boost{namespace multi{namespace lapack{

using blas::filling;

namespace{

using ::core::potrf;

template<class Iterator>
auto potrf(filling t, Iterator first, Iterator last)
->decltype(potrf(static_cast<char>(t), typename std::iterator_traits<Iterator>::difference_type{}, base(first), stride(first), std::declval<int&>()), Iterator{})
{
	assert( stride(first) == stride(last) );
	assert( first->stride() == 1 );
	auto n = std::distance(first, last);
//	auto lda = stride(first);
	int info;
	potrf(static_cast<char>(t), n, base(first), stride(first), info);
	assert( info >= 0 );
	return info==0?last:first + info;
}
}

template<class A2D>
auto potrf(filling t, A2D&& A)
->decltype(potrf(t, begin(A), end(A)), A({0, 1}, {0, 1}))
{
	using blas::flip;
	if(stride(A)==1){
		auto last = potrf(flip(t), begin(rotated(A)), end(rotated(A)));
		using std::distance;
		return A({0, distance(begin(rotated(A)), last)}, {0, distance(begin(rotated(A)), last)});
	}
	auto last = potrf(t, begin(A), end(A));
	using std::distance;
	return A({0, distance(begin(A), last)}, {0, distance(begin(A), last)});
}

template<class A>
struct hermitic_t : private A{
	using underlying_type = A;
	underlying_type const& underlying()const &{return *this;}
	underlying_type& underlying()&{return *this;}
	underlying_type&& underlying()&&{return std::move(*this);}
	blas::filling side;
	hermitic_t(A const& a, blas::filling side) : A{a}, side{side}{}
	using A::size;
};

template<class A> hermitic_t<std::decay_t<decltype(std::declval<A>()())>> hermitic(blas::filling side, A&& a){
	return {a(), side};
}

template<class A2D>
NODISCARD("result is returned because third argument is const")
auto potrf(filling t, A2D const& A)
->decltype(potrf(t, decay(A)), decay(A)){
	auto ret = decay(A);
	auto last = potrf(t, ret); assert( size(last) == size(ret) );
	return ret;
}

template<class HA>
NODISCARD("result is returned because third argument is const")
decltype(auto) potrf(HA&& ha){
	return hermitic(ha.side, potrf(ha.side, std::forward<HA>(ha).underlying()));//static_cast<typename HA::underlying_type&>(ha)));
}

// orthonormalize rows
template<class A> auto onrm(A&& a, filling f = filling::upper)
->decltype(trsm(flip(f), hermitized(potrf(f, herk(f, a))), std::forward<A>(a))){assert(size(a) <= size(rotated(a)));
	return trsm(flip(f), hermitized(potrf(f, herk(f, a))), std::forward<A>(a));		
}

template<class A, class B> auto onrm(A&& a, B&& buffer, filling f = filling::upper)
->decltype(trsm(flip(f), hermitized(potrf(f, herk(f, a, buffer))), std::forward<A>(a))){assert(size(a) <= size(rotated(a)));
	return trsm(flip(f), hermitized(potrf(f, herk(f, a, buffer))), std::forward<A>(a));
}

//template<class A2D>
//decltype(auto) potrf(A2D&& A){return potrf(blas::detect_triangular(A), A);}

}}}

#if _TEST_MULTI_ADAPTORS_LAPACK_POTRF

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi lapack adaptor potrf"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<cmath> // std::isnan

namespace multi = boost::multi;
namespace lapack = multi::lapack;

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout << C[i][j] << ' ';
		cout << std::endl;
	}
	return cout << std::endl;
}

BOOST_AUTO_TEST_CASE(lapack_potrf, *boost::unit_test::tolerance(0.00001) ){
	using complex = std::complex<double>; complex const I{0, 1};
{
	multi::array<complex, 2> A = {
		{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		{NAN    , 167.381               , 126.746 + 0.0327519*I},
		{NAN    , NAN                   , 167.231              }
	};
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is hermitic in upper triangular (implicit below)
	BOOST_TEST( real(A[1][2]) == 3.78646 );
	BOOST_TEST( imag(A[1][2]) == 0.0170734 );
//	BOOST_TEST( std::isnan(norm(A[2][1])) );
}
{
	multi::array<complex, 2> A =
		{{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		 {NAN, 167.381, 126.746 + 0.0327519*I},
		 {NAN, NAN , 167.231}}
	;
	multi::array<complex, 2> At = rotated(A);
	auto&& Att = rotated(At);
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, Att); // A is hermitic in the upper triangular (implicit hermitic below)
	BOOST_TEST( real(Att[1][2]) == 3.78646 );
	BOOST_TEST( imag(Att[1][2]) == 0.0170734 );
//	BOOST_TEST( std::isnan(norm(Att[2][1])) );
}
{
	multi::array<complex, 2> A =
		{{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		 {NAN, 167.381, 126.746 + 0.0327519*I},
		 {NAN, NAN , 167.231}}
	;
	using lapack::potrf;
	using lapack::filling;
	potrf(filling::upper, A); // A is hermitic in the upper triangular (implicit hermitic below)
	BOOST_TEST( real(A[1][2]) == 3.78646 );
	BOOST_TEST( imag(A[1][2]) == 0.0170734 );
//	BOOST_TEST( std::isnan(A[2][1]) );
}
{
	multi::array<complex, 2> A =
		{{190., 126., 125.},
		 {NAN , 1110., 122.},
		 {NAN , NAN , 1350.}}
	;
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is the upper triangle (implicit hermitic/symmetric below), A becomes upper triangular with implicit zeros
	BOOST_TEST( real(A[1][2]) == 1.22058 );
//	BOOST_TEST( std::isnan(norm(A[2][1])) );
}
{
	multi::array<double, 2> A =
		{{190., 126., 125.},
		 {NAN , 1110., 122.},
		 {NAN , NAN , 1350.}}
	;
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is the upper triangle (implicit hermitic/symmetric below), A becomes upper triangular with implicit zeros
	BOOST_TEST( A[1][2] == 1.22058 );
//	BOOST_TEST( std::isnan(norm(A[2][1])) );
}
{
	multi::array<double, 2> A =
		{{190., 126., 125.},
		 {NAN , 1110., 122.},
		 {NAN , NAN , 1350.}}
	;
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::lower, rotated(A)); // A is the upper triangle (implicit hermitic/symmetric below), A becomes upper triangular with implicit symmetry
	print(A);
	BOOST_TEST( A[1][2] == 1.22058 );
//	BOOST_TEST( std::isnan(norm(A[2][1])) );
}
{
	multi::array<complex, 2> const A =
		{{190., 126., 125.},
		 {NAN , 1110., 122.},
		 {NAN , NAN , 1350.}}
	;
	using lapack::filling;
	using lapack::potrf;
	auto B = potrf(filling::upper, A);
	print(B);
	BOOST_TEST( real(B[1][2]) == 1.22058 );
}
}

#endif
#endif

