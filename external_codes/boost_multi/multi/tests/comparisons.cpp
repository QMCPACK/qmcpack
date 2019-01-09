#ifdef COMPILATION_INSTRUCTIONS
c++ -O3 -std=c++14 -Wall -Wextra -Wpedantic -DNDEBUG $0 -o $0.x && $0.x $@ && rm -f $0.x;exit
#endif

#include<iostream>
#include "../../multi/array.hpp"

namespace multi = boost::multi;
using std::cout;

int main(){

	multi::array<double, 3> A {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};

	multi::array_ref<double, 3> AR(extensions(A), A.data());
	multi::array_cref<double, 3> AC(extensions(A), A.data());

	assert( A == A and AR == A and AR == AC );
	assert( A[0] == A[2] and AR[0]==A[2] and AR[0]==AC[2] );
	assert( not (A[0] != A[2]) and not (AR[0] != AR[2]) );
	assert( A[0] <= A[1] and AR[0] <= A[1] and AC[0] <= AC[1] );
	assert( A[0][0] <= A[0][1] and AR[0][0] <= A[0][1] );
	assert( A[1][0][0] == 1.2 and AR[1][0][0] == 1.2 and AC[1][0][0] == 1.2 );
	assert( A[0][0][0] == 11.2 and AR[0][0][0] == 11.2 and AC[0][0][0] == 11.2);
	swap(AR[0], AR[1]);
	swap(A[0], A[1]);
	assert( A[0] > A[1] );
	A[1] = A[0];
	assert( A[0] == A[1] );

	using std::begin;
	using std::end;

	assert( begin(A) < end(A) );
	assert( cbegin(A) < cend(A) );
	assert( crbegin(A) < crend(A) );
	assert( crend(A) > crbegin(A) );
	assert( end(A) - begin(A) == size(A));
	assert( rend(A) - rbegin(A) == size(A));

}

