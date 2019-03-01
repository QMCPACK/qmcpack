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
		multi::blas::scal(2., A[2]); // dscal
		assert( A[0][2] == 3. and A[2][2] == 11.*2. );
	}
	{
		multi::array<std::complex<double>, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		assert( A[0][2] == 3. and A[2][2] == 11. );
		multi::blas::scal(std::complex<double>(2.), A[2]); // zscal
		assert( A[0][2] == 3. and A[2][2] == 11.*2. );
		multi::blas::scal(1./2, A[2]); // zdscal
		assert( A[0][2] == 3. and A[2][1] == 10. and A[2][2] == 11. );
		multi::blas::scal(2., begin(A[2]), begin(A[2]) + 2);
		assert( A[0][2] == 3. and A[2][1] == 20. and A[2][2] == 11. );
	}
}

