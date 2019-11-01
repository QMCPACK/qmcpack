#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>
#include<cassert>
#include<vector>

#include "../array.hpp"

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){
	std::vector<double> buffer(100);
	multi::array_ref<double, 2, std::vector<double>::iterator> A(buffer.begin(), {10, 10});
	A[1][1] = 9;
	assert(A[1][1] == 9);
	assert(buffer[11]==9);

	A[2]; // requires operator+ 
	A[1][1]; // requires operator*
}

