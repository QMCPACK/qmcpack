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
		multi::blas::copy(A[0], A[2]); // dcopy
		assert( A[0][2] == 3. and A[2][2] == 3. );
		multi::blas::copy(begin(A[1]), end(A[1]), begin(A[2])); // dcopy
		assert( A[1][3] == 8. and A[2][3] == 8. );
		auto AR3 = multi::blas::copy(rotated(A)[3]); // dcopy
		assert( AR3[1] == A[1][3] );
	}
}

