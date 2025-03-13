#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x$OXX `pkg-config --libs blas lapack` -lboost_unit_test_framework&&$0x$OXX -x 0&&rm $0x$OXX;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_ADAPTORS_LAPACK_GEQRF_HPP
#define MULTI_ADAPTORS_LAPACK_GEQRF_HPP

#include "../lapack/core.hpp"
#include "../blas/filling.hpp"

#include "../../config/NODISCARD.hpp"

#include<cassert>

namespace boost{namespace multi{namespace lapack{

using blas::filling;

template<class Context, class A, class TAU, class WORK>
A&& geqrf(Context&& ctxt, A&& a, TAU&& tau, WORK&& work){
//  assert( stride(~a) == 1);
	assert( size(tau) == std::min(size(~a), size(a)) );
	int info = -1;
	geqrf_(std::forward<Context>(ctxt), size(~a), size(a), a.base(), stride(a), tau.base(), work.data(), work.size(), info);
	assert(info == 0);
	return std::forward<A>(a);
}

//using ::core::syev;
//using ::core::geqrf;

#if 0
template<class Array2D, class Array1D, class Array1DW>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w, Array1DW&& work)
->decltype(syev('V', uplo==blas::filling::upper?'L':'U', size(a), base(a), stride(a), base(w), base(work), size(work), std::declval<int&>()), a({0l, 1l}, {0l, 1l}))
{
	assert( size(work) >= std::max(1l, 3*size(a)-1l) );
	assert( size(a) == size(w) );
	assert( stride(w)==1 );
	assert( stride(work)==1 );
	if(size(a)==0) return std::forward<Array2D>(a)();
	int info = -1;
	     if(stride(rotated(a))==1) syev('V', uplo==blas::filling::upper?'L':'U', size(a), base(a), stride(        a ), base(w), base(work), size(work), info);
	else if(stride(        a )==1) syev('V', uplo==blas::filling::upper?'U':'L', size(a), base(a), stride(rotated(a)), base(w), base(work), size(work), info);
	else                           assert(0); // case not contemplated by lapack
	if(info < 0) assert(0); // bad argument
	return std::forward<Array2D>(a)({0, size(a)-info}, {0, size(a)-info});
}

template<class Array2D, class Array1D, class Array1DW = typename std::decay_t<Array1D>::decay_type>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w)
->decltype(syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1l, 3*size(a)-1l), get_allocator(w)))){
	return syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1l, 3*size(a)-1l), get_allocator(w)));}// TODO obtain automatic size from lapack info routine

template<class Array2D, class Array1D>
NODISCARD("because input array is const, output gives eigenvectors")
typename Array2D::decay_type syev(blas::filling uplo, Array2D const& a, Array1D&& w){
	auto ret = a.decay();
	if(syev(uplo, ret, std::forward<Array1D>(w)).size() != a.size()) assert(0); // failed
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
#endif

}}}
