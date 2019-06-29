#ifdef COMPILATION_INSTRUCTIONS
$CXX -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>
#include<cassert>

#include "../array.hpp"


using std::cout; using std::cerr;
namespace multi = boost::multi;

namespace min{
template<class T> struct ptr{ // minimalistic pointer
	T* impl_;
	T& operator*() const{return *impl_;}
	auto operator+(std::ptrdiff_t n) const{return ptr{impl_ + n};}
//	T& operator[](std::ptrdiff_t n) const{return impl_[n];} // optional
};
}

int main(){
	double* buffer = new double[100];
	multi::array_ref<double, 2, min::ptr<double> > CC(min::ptr<double>{buffer}, {10, 10});
	CC[2]; // requires operator+ 
	CC[1][1]; // requires operator*
	CC[1][1] = 9;
	assert(CC[1][1] == 9);
}

