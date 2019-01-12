#ifdef COMPILATION_INSTRUCTIONS
c++ -std=c++14 -Wall -Wextra -Wpedantic -lblas -DADD_ $0 -o $0x.x  && time $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../blas.hpp"

#include "../../../array.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;

int main(){
	{
		using Z = std::complex<double>; Z const I(0.,1.);
		multi::array<Z, 2> const A = {
			{1. + 2.*I,  2.,  3.,  4.},
			{5.,  6. + 3.*I,  7.,  8.},
			{9., 10., 11.+ 4.*I, 12.}
		};
		using multi::blas::iamax; // izmax, uses zero-based indexing
		assert(iamax(A[1]) == std::distance(begin(A[1]), std::max_element(begin(A[1]), end(A[1]), [](auto&& a, auto&& b){return std::abs(real(a)) + std::abs(imag(a)) < std::abs(real(b)) + std::abs(imag(b));})));
	}
}

