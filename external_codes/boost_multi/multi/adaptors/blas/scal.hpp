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

template<class T, class It>
void scal(T a, It first, It last){
	assert(stride(first) == stride(last));
	scal(std::distance(first, last), a, base(first), stride(first));
}

template<class T, class X1D>
X1D&& scal(T a, X1D&& m){
	assert( offset(m) == 0 );
	scal(a, begin(m), end(m));
	return std::forward<X1D>(m);
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
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};

	using multi::blas::scal;
	auto&& S = scal(2., A.rotated(1)[1]);

	assert( A[2][1] == 20. );
	assert( S[0] == 4. );

}

#endif
#endif

