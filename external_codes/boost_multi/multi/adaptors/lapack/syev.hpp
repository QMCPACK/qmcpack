#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -D_TEST_MULTI_ADAPTORS_LAPACK_SYEV $0.cpp -o $0x `pkg-config --libs blas lapack` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_ADAPTORS_LAPACK_SYEV_HPP
#define MULTI_ADAPTORS_LAPACK_SYEV_HPP

#include "../lapack/core.hpp"
#include "../blas/filling.hpp"

#include "../../config/NODISCARD.hpp"

#include<cassert>

namespace boost{namespace multi{namespace lapack{

using blas::filling;

using ::core::syev;

template<class Array2D, class Array1D, class Array1DW>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w, Array1DW&& work)
->decltype(syev('V', uplo==blas::filling::upper?'L':'U', size(a), base(a), stride(a), base(w), base(work), size(work), std::declval<int&>()), a({0l, 1l}, {0l, 1l}))
{
	assert( size(work) >= std::max(1l, 3*size(a)-1l) );
	assert( size(a) == size(w) );
	assert( stride(w)==1 );
	assert( stride(work)==1 );
	if(size(a)==0) return a();
	int info = -1;
	     if(stride(rotated(a))==1) syev('V', uplo==blas::filling::upper?'L':'U', size(a), base(a), stride(        a ), base(w), base(work), size(work), info);
	else if(stride(        a )==1) syev('V', uplo==blas::filling::upper?'U':'L', size(a), base(a), stride(rotated(a)), base(w), base(work), size(work), info);
	else                           assert(0); // case not contemplated by lapack
	if(info < 0) assert(0); // bad argument
	return a({0, size(a)-info}, {0, size(a)-info});
}

template<class Array2D, class Array1D, class Array1DW = typename std::decay_t<Array1D>::decay_type>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w)
->decltype(syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1l, 3*size(a)-1l), get_allocator(w)))){
	return syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1l, 3*size(a)-1l), get_allocator(w)));}// TODO obtain automatic size from lapack info routine

template<class Array2D, class Array1D>
NODISCARD("because input array is const, output gives eigenvectors")
typename Array2D::decay_type syev(blas::filling uplo, Array2D const& a, Array1D&& w){
	auto ret = a.decay();
	auto l = syev(uplo, ret, std::forward<Array1D>(w));
	if(size(l) != size(a)) assert(0); // failed
	return ret;
}

template<class Array2D>
NODISCARD("because input array is const, output gives eigenvalues")
auto syev(blas::filling uplo, Array2D&& a){
	multi::array<typename std::decay_t<Array2D>::element_type, 1, decltype(get_allocator(a))> eigenvalues(size(a), get_allocator(a));
	syev(uplo, std::forward<Array2D>(a), eigenvalues);
	return eigenvalues;
}

template<class Array2D>
NODISCARD("because input array is const, output gives a structured binding of eigenvectors and eigenvactor")
auto syev(blas::filling uplo, Array2D const& a){
	struct{
		typename Array2D::decay_type eigenvectors;
		typename Array2D::value_type eigenvalues;
	} ret{a, {size(a), get_allocator(a)}};
	auto&& l = syev(uplo, ret.eigenvectors, ret.eigenvalues);
	assert( size(l) == size(a) );
	return ret;
}

}}}

#if _TEST_MULTI_ADAPTORS_LAPACK_SYEV

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi lapack adaptor syev"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"

#include<cmath> // std::isnan
#include<iostream>
#include<algorithm> // std::max

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

BOOST_AUTO_TEST_CASE(lapack_syev, *boost::unit_test::tolerance(0.00001) ){
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));
	multi::lapack::syev(multi::blas::filling::upper, A, W, WORK);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::upper, A, W);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::lower, rotated(A), W);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	namespace lapack = multi::lapack;
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto A_copy = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( A_copy[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 0.},
		{NAN    , 167.381, 0.},
		{NAN    , NAN    , 0.}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( W[0]==0. );
}
{
	multi::array<double, 2> A = {
		{1. , 1.,  1.},
		{NAN, 2 ,  1.},
		{NAN, NAN, 1.}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	print(A_ref);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( W[0]==0. );
}
{
	multi::array<double, 2> A = {{5.}};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[0][0] == 1. );
	BOOST_TEST( W[0]==5. );
}
{
	multi::array<double, 2> A;
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
}
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto sys = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( sys.eigenvectors[2][1] == -0.579092 );
	BOOST_TEST( sys.eigenvalues[1] == 42.2081 );
}
#if __cpp_structured_bindings
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto [eigenvecs, eigenvals] = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( eigenvecs[2][1] == -0.579092 );
	BOOST_TEST( eigenvals[1] == 42.2081 );
}
#endif

}
#endif
#endif


