#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_SWAP -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_SWAP_HPP
#define MULTI_ADAPTORS_BLAS_SWAP_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It1, class It2>
It2 swap(It1 first, It2 last, It2 first2){
	assert(stride(first) == stride(last));
	using std::distance;
	auto d = distance(first, last);
	swap(d, base(first), stride(first), base(first2), stride(first2));
	return first2 + d;
}

template<class X1D, class Y1D>
Y1D&& swap(X1D&& x, Y1D&& y){
	assert( size(x) == size(y) );
	assert( offset(x) == 0 and offset(y) == 0 );
	swap( begin(x), end(x), begin(y) );
	return std::move(y);
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_SWAP

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
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
//	using multi::blas::swap;
	multi::blas::swap(rotated(A)[1], rotated(A)[3]); // can ambiguate with (friend) multi::swap
}

#endif
#endif

