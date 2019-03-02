#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_DOT -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_DOT_HPP
#define MULTI_ADAPTORS_BLAS_DOT_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class R, class It1, class It2>
auto dot(It1 first1, It1 last1, It2 first2){
	assert( stride(first1) == stride(first2) );
	auto d = std::distance(first1, last1);
	return dot<R>(d, base(first1), stride(first1), base(first2), stride(first2));
}

template<class It1, class It2>
auto dot(It1 first1, It1 last1, It2 first2){
	assert( stride(first1) == stride(first2) );
	auto d = std::distance(first1, last1);
	return dot(d, base(first1), stride(first1), base(first2), stride(first2));
}

template<class R, class X1D, class Y1D>
auto dot(X1D const& x, Y1D const& y){
	assert( size(x) == size(y) );
	assert( not offset(x) and not offset(y) );
	return dot<R>(begin(x), end(x), begin(y));
}

template<class X1D, class Y1D>
auto dot(X1D const& x, Y1D const& y){
	assert( size(x) == size(y) );
	assert( not offset(x) and not offset(y) );
	return dot(begin(x), end(x), begin(y));
}

template<class It1, class Size, class It2>
auto dotu(It1 first1, Size n, It2 first2){
	return dotu(n, base(first1), stride(first1), base(first2), stride(first2));
}

template<class It1, class It2>
auto dotu(It1 first1, It1 last1, It2 first2){
	assert( stride(first1) == stride(last1) );
	return dotu(first1, std::distance(first1, last1), first2);
}

template<class X1D, class Y1D>
auto dotu(X1D const& x, Y1D const& y){
	assert( size(x) == size(y) );
	assert( not offset(x) and not offset(y) );
	return dotu(begin(x), end(x), begin(y));
}

template<class It1, typename Size, class It2>
auto dotc_n(It1 first1, Size n, It2 first2){
	dotc(n, base(first1), stride(first1), base(first2), stride(first2));
	return first2 + n;
}

template<class It1, class It2>
auto dotc(It1 first1, It1 last1, It2 first2){
	assert( stride(first1) == stride(last1) );
	return dotc_n(first1, std::distance(first1, last1), first2);
}

template<class X1D, class Y1D>
auto dotc(X1D const& x, Y1D const& y){
	assert( size(x) == size(y) );
	assert( not offset(x) and not offset(y) );
	return dotc(begin(x), end(x), begin(y));
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_DOT

#include "../../array.hpp"
#include "../../utility.hpp"
#include "../blas/nrm2.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;
namespace multi = boost::multi;

int main(){
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	using multi::blas::dot;
	auto d = dot(cA[1], cA[2]);
	assert(d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.));

	using multi::blas::nrm2;
	assert( std::sqrt(dot(cA[1], cA[1])) == nrm2(cA[1]) );
}

#endif
#endif

