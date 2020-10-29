#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#include<iostream>
#include<cassert>

#include "../array.hpp"

using std::cout; using std::cerr;
namespace multi = boost::multi;

namespace minimalistic{
template<class T> struct ptr : std::iterator_traits<T*>{ // minimalistic pointer
	T* impl_;
	ptr(T* impl) : impl_{impl}{}
	typename ptr::reference operator*() const{return *impl_;}
	auto operator+(typename ptr::difference_type n) const{return ptr{impl_ + n};}
//	T& operator[](std::ptrdiff_t n) const{return impl_[n];} // optional
};

template<class T> struct ptr2 : std::iterator_traits<T*>{ // minimalistic pointer
	T* impl_;
	ptr2(T* impl) : impl_{impl}{}
	explicit ptr2(ptr<T> p) : impl_{p.impl_}{} 
	typename ptr2::reference operator*() const{return *impl_;}
	auto operator+(typename ptr2::difference_type n) const{return ptr2{impl_ + n};}
//	T& operator[](std::ptrdiff_t n) const{return impl_[n];} // optional
};

}

struct X{
	int a1;
	double a2;
	double b;
};

int main(){

	int X::* s = &X::a1;
	X x{1, 2.5, 1.2};
	assert( x.*s == x.a1 );

//	X X::*ss = &X::X; 

	double* buffer = new double[100];
	multi::array_ref<double, 2, minimalistic::ptr<double> > CC(minimalistic::ptr<double>{buffer}, {10, 10});
	CC[2]; // requires operator+ 
	CC[1][1]; // requires operator*
	CC[1][1] = 9;
	assert(CC[1][1] == 9);

	using multi::static_array_cast;
	auto&& CC2 = static_array_cast<double, minimalistic::ptr2<double>>(CC);
	assert( &CC2[1][1] == &CC[1][1] );

}

