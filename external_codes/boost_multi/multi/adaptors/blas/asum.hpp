#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_SCAL -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_DOT_HPP
#define MULTI_ADAPTORS_BLAS_DOT_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It>
auto asum(It first, It last){
	assert( stride(first) == stride(last) );
	return asum(std::distance(first, last), base(first), stride(first));
}

template<class X1D> 
auto asum(X1D const& x){
	assert( not offset(x) );
//	using blas::asum;
	return asum(size(x), origin(x), stride(x));
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_SCAL

#include "../../array.hpp"
#include "../../utility.hpp"

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
	multi::array<std::complex<double>, 2> A = cA;
	using multi::blas::asum;
	assert(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a + std::abs(real(b)) + std::abs(imag(b));}));

}

#endif
#endif

