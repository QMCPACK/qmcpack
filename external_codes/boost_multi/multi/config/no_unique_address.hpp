#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -D_TEST_MULTI_CONFIG_NO_UNIQUE_ADDRESS $0.cpp -o$0x&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP
#define MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if __has_cpp_attribute(no_unique_address) >=201803
	#define NO_UNIQUE_ADDRESS [[no_unique_address]]
	#define no_unique_address_ no_unique_address
#else
	#define NO_UNIQUE_ADDRESS
#endif

////////////////////////////////////////////////////////////////////////////////
#ifdef _TEST_MULTI_CONFIG_NO_UNIQUE_ADDRESS

class A{};

class B{
	NO_UNIQUE_ADDRESS A x;
	double y;
};

int main(){
	static_assert( sizeof(B) == sizeof(double) , "may fail with no unique feauture"); // for example fails with clang++-8
}
#endif
#endif


