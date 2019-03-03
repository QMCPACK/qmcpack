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
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		assert( A[0][2] == 3. and A[2][2] == 11. );
		multi::blas::swap(A[0], A[2]); // blas swap
		assert( A[0][2] == 11. and A[2][2] == 3. );
		             swap(A[0], A[2]); // built-in swap
		assert( A[0][2] == 3. and A[2][2] == 11. );
	}
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		assert( A[0][0] == 1. and A[0][3] == 4. );
		multi::blas::swap(rotated(A)[0], rotated(A)[3]); // blas swap (deep)
		assert( A[0][0] == 4. and A[0][3] == 1. );
		             swap(rotated(A)[0], rotated(A)[3]); // built-in swap (deep)
		assert( A[0][0] == 1. and A[0][3] == 4. );
	}
	{
		std::complex<double> const I(0., 1.);
		multi::array<std::complex<double>, 2> A = {
			{1.+ 2.*I,  2.,  3.,  4. + 3.*I},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		assert( A[0][0] == 1.+ 2.*I and A[0][3] == 4. + 3.*I );
		multi::blas::swap(rotated(A)[0], rotated(A)[3]); // blas swap (deep)
		assert( A[0][0] == 4. + 3.*I and A[0][3] == 1.+ 2.*I );
		             swap(rotated(A)[0], rotated(A)[3]); // built-in swap (deep)
		assert( A[0][0] == 1.+ 2.*I and A[0][3] == 4. + 3.*I );
	}
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		assert( A[0][2] == 3. and A[2][2] == 11. );
		auto it = multi::blas::swap(begin(A[0]), end(A[0]) - 1, begin(A[2])); // blas swap
		assert( it == end(A[2]) - 1 );
		assert( A[0][2] == 11. and A[2][2] == 3. );
		using std::swap_ranges;
		      swap_ranges(begin(A[0]), end(A[0]), begin(A[2])); // built-in swap
		assert( A[0][2] == 3. and A[2][2] == 11. );
	}
}

