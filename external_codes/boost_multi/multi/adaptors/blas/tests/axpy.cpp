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
		auto const AC = A;
		multi::array<double, 1> const B = A[2];
		using multi::blas::axpy;
		axpy(2., B, A[1]); // daxpy
		assert( A[1][2] == 2.*B[2] + AC[1][2] );
	}
	{
		using Z = std::complex<double>;
		multi::array<Z, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		auto const AC = A;
		multi::array<Z, 1> const B = A[2];
		using multi::blas::axpy;
		axpy(2., B, A[1]); // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		assert( A[1][2] == 2.*B[2] + AC[1][2] );
	}
}

