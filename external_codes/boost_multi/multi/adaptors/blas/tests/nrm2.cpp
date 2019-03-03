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
		multi::array<double, 2> const A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		using multi::blas::dot;
		using multi::blas::nrm2;
		assert( nrm2(A[1]) == std::sqrt(dot(A[1], A[1])) );
	}
	{
		using Z = std::complex<double>; Z const I(0.,1.);
		multi::array<Z, 1> A = {I, 2.*I, 3.*I};
		using multi::blas::dot;
		using multi::blas::dotc;
		using multi::blas::dotu;
		using multi::blas::nrm2;
		auto n = nrm2(A);
		static_assert( std::is_same<decltype(n), double>{}, "!");
		assert( std::abs(n*n - real(dot(A, A))) < 1e-14 );
		assert( std::abs(n*n - real(dotc(A, A))) < 1e-14 );
		assert( std::abs(n*n - real(dotu(A, A))) > 1e-14 );
	}
}

