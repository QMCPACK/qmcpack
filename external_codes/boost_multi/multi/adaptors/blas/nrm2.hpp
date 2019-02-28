#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_NRM2 -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_NRM2_HPP
#define MULTI_ADAPTORS_BLAS_NRM2_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It, class Size>
auto nrm2_n(It first, Size n)
->decltype(nrm2(n, base(first), stride(first)))
{
	return nrm2(n, base(first), stride(first));
}

template<class It>
auto nrm2(It first, It last)
->decltype(nrm2_n(first, std::distance(first, last)))
{
	assert( stride(first) == stride(last) );
	return nrm2_n(first, std::distance(first, last));
}

template<class X1D> 
auto nrm2(X1D const& x)
->decltype(nrm2(begin(x), end(x)))
{ 
	assert( not offset(x) );
	return nrm2(begin(x), end(x));
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_NRM2

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
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	using multi::blas::nrm2;
	using multi::blas::dot;
	assert( nrm2(cA[1]) == std::sqrt(dot(cA[1], cA[1])) );
}

#endif
#endif

