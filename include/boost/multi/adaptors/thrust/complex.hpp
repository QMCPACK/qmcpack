#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2020

#ifndef MULTI_ADAPTORS_THRUST_COMPLEX_HPP
#define MULTI_ADAPTORS_THRUST_COMPLEX_HPP

#include<type_traits>
#include<thrust/complex.h>

namespace std{
	template<class T> struct is_trivially_copy_constructible<thrust::complex<T>                    > : std::true_type{};
	template<class T> struct is_trivially_constructible     <thrust::complex<T>, thrust::complex<T>> : std::true_type{};

	template<class T> struct is_trivially_assignable<thrust::complex<T>&, thrust::complex<T>>        : std::true_type{};
	template<class T> struct is_trivially_assignable<thrust::complex<T>&, thrust::complex<T> const > : std::true_type{};
	template<class T> struct is_trivially_assignable<thrust::complex<T>&, thrust::complex<T> const&> : std::true_type{};
	template<class T> struct is_trivially_assignable<thrust::complex<T>&, thrust::complex<T>      &> : std::true_type{};			

// this one is controversional because it changes the behavior of initialization of values in multi (even in the cpu)
	template<class T> struct is_trivially_default_constructible<thrust::complex<T>> : is_trivially_default_constructible<T>{};
}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_COMPLEX

#include "../../array.hpp"

#include<cassert>
#include<complex>

namespace multi = boost::multi;

template<class T> void what(T&&)=delete;

int main(){
	using complex = thrust::complex<double>;

	static_assert( std::is_trivially_default_constructible<complex>{}, "!");
	static_assert( std::is_trivially_copy_constructible<complex>{}   , "!");
	static_assert( std::is_trivially_assignable<complex&, complex>{} , "!");

	multi::array<complex, 2> A = {
		{ { 1., 2.}, { 3., 4.} },
		{ {22.,33.}, { 5., 9.} }
	};
}

#endif
#endif


