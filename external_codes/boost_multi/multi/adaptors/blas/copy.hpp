#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_COPY -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_COPY_HPP
#define MULTI_ADAPTORS_BLAS_COPY_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It, class Size, class OutIt>
OutIt copy_n(It first, Size n, OutIt d_first){
	copy(n, base(first), stride(first), base(d_first), stride(d_first));
	return d_first + n;
}

template<class It1, class OutIt>
OutIt copy(It1 first, It1 last, OutIt d_first){
	assert( stride(first) == stride(last) );
	return copy_n(first, std::distance(first, last), d_first);
}

template<class X1D, class Y1D>
Y1D&& copy(X1D const& x, Y1D&& y){
	assert( size(x) == size(y) );
	assert( offset(x) == 0 and offset(y) == 0 );
	auto it = copy(begin(x), end(x), begin(y));
	assert( it == end(y));
	return std::forward<Y1D>(y);
}

template<class X1D>
multi::array<typename X1D::value_type, 1> copy(X1D const& x){
	assert( not offset(x) );
	return copy(x, multi::array<typename X1D::value_type, 1>({0, size(x)}, 0.));
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_COPY

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

#include "../blas/dot.hpp"

using std::cout;
namespace multi = boost::multi;

int main(){
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	multi::array<double, 1> B = {5., 6., 7., 8.};
	using multi::blas::copy;
	copy(A, B);
	assert( B == A );
}

#endif
#endif

