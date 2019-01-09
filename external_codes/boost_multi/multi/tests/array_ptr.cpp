#ifdef COMPILATION_INSTRUCTIONS
c++ -O3 -std=c++14 `#-DNDEBUG` -Wall -Wextra $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>
#include<vector>
#include "../array.hpp"

namespace multi = boost::multi;

int main(){
	multi::array<double, 3> A {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};
	
	auto p = &A[1];
	auto p2 = p + 1;
	assert( &(*p )[0][0] == &A[1][0][0] );
	assert( &(*p2)[0][0] == &A[2][0][0] ); // this is true only because A is contiguous
	p2 = p;
	assert( p2 == p );
	auto p3 = &A[2][1];
	assert( &(*p3)[1] == &A[2][1][1] );
	assert( &p3->operator[](1) == &A[2][1][1] );
}

template<class It>
class F{
	It first;
	operator(It){ return std::distance() ; }
};

for(auto it = first; it != last; ++it){
	f(*it);
}

for(auto it = first; it != last; ++it){
	f(it);
}



for(auto it = first; it != last; ++it){
	*it = f();
}

