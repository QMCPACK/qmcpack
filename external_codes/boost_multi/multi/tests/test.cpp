#ifdef COMPILATION_INSTRUCTIONS
${CXX} -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && time $0.x $@ && rm -f $0.x; exit
#endif

#include "../array_ref.hpp"
#include "../array.hpp"

#include<iostream>

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){
	multi::array<int, 2> A({2, 2});
        std::cout<<(0 != A[0][0]) <<std::endl;
}

