#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic $0 -o $0.x && $0.x $@ && rm -f $0.x;exit
#endif

#include<iostream>
#include "../array.hpp"

namespace multi = boost::multi;
using std::cout;

int main(){

	multi::array<double, 3> A 
		#if __INTEL_COMPILER
		= (double[3][2][2])
		#endif
		{
			{{ 1.2,  1.1}, { 2.4, 1.}},
			{{11.2,  3.0}, {34.4, 4.}},
			{{ 1.2,  1.1}, { 2.4, 1.}}
		};
	
	assert( A[0] < A[1] );
	cout << A[0][0][1] <<std::endl;
	cout << A[1][0][1] <<std::endl;
	swap( A[0], A[1] );
	cout << A[0][0][1] <<std::endl;
	cout << A[1][0][1] <<std::endl;
	swap( A[0], A[1] );
	
	multi::array_ref<double, 3> AR(data(A), extensions(A));
	multi::array_cref<double, 3> AC(data(A), extensions(A));

	assert( A == A and AR == A and AR == AC );
	assert( A[0] == A[2] and AR[0]==A[2] and AR[0]==AC[2] );
	assert( not (A[0] != A[2]) and not (AR[0] != AR[2]) );
	assert( A[0] <= A[1] and AR[0] <= A[1] and AC[0] <= AC[1] );
	assert( A[0][0] <= A[0][1] and AR[0][0] <= A[0][1] );

	cout<< A[1][0][0] <<std::endl;
	assert( A[1][0][0] == 11.2 );
	assert( AR[1][0][0] == 11.2 );
	assert( AC[1][0][0] == 11.2 );
	assert( A[0][0][0] == 1.2 and AR[0][0][0] == 1.2 and AC[0][0][0] == 1.2);
	swap(AR[0], AR[1]);

	using std::begin;
	using std::end;

	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
//	assert( crbegin(A) < crend(A) );
//	assert( crend(A) > crbegin(A) );
	assert( end(A) - begin(A) == size(A));
	assert( rend(A) - rbegin(A) == size(A));

}

