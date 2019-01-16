#ifdef COMPILATION_INSTRUCTIONS
clang++ -O3 -std=c++17 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0.x && time $0.x $@ && rm -f $0.x; exit
#endif

#include<iostream>

////////////////////////////////////////////////////////////////////////////////
// minimal fancy ptr
////////////////////////////////////////////////////////////////////////////////
namespace fancy{

template<class T> struct ptr{
	static T value;
	using difference_type = std::ptrdiff_t;
	using value_type = T;
	using pointer = T*;
	using reference = T&;
	using iterator_category = std::random_access_iterator_tag;
	reference operator*() const{return value;}
	ptr operator+(difference_type) const{return *this;}
	ptr& operator+=(difference_type){return *this;}
	friend difference_type operator-(ptr const&, ptr const&){return 0;}
};
template<> double ptr<double>::value = 42.;

template<class T>
ptr<T> copy(ptr<T>, ptr<T>, ptr<T>){
	std::cerr << "called Pointer-based copy(ptr, ptr, ptr)" << std::endl; 
	return {};
}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// multi-fancy glue, where should this be? In boost/multi/adaptors/fancy.hpp if anything, or in user code if it is very special
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../array_ref.hpp"

namespace boost::multi{

template<class T>
void copy(
	multi::array_iterator<T, 1, fancy::ptr<T>> first, 
	multi::array_iterator<T, 1, fancy::ptr<T>> last, 
	multi::array_iterator<T, 1, fancy::ptr<T>> dest
){
	assert( stride(first) == stride(last) );
	std::cerr << "called Iterator-based copy(it1D, it1D, it1D) with strides " << stride(first) << " " << stride(dest) << std::endl;
}

};

////////////////////////////////////////////////////////////////////////////////
// user code
////////////////////////////////////////////////////////////////////////////////

using std::cout; using std::cerr;
namespace multi = boost::multi;

int main(){

	fancy::ptr<double> p1;
	multi::array_ref<double, 2, fancy::ptr<double>> A(p1, multi::iextensions<2>{5, 5});

	fancy::ptr<double> p2;
	multi::array_ref<double, 2, fancy::ptr<double>> B(p2, multi::iextensions<2>{5, 5});

	rotated(A)[1] = B[1];

}

