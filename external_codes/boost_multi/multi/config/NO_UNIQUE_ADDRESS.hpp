// © Alfredo A. Correa 2019-2021

#ifndef MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP
#define MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if __has_cpp_attribute(no_unique_address) >=201803 and not defined(__NVCC__) and not defined(__PGI)
	#define MULTI_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
	#define MULTI_NO_UNIQUE_ADDRESS
#endif

////////////////////////////////////////////////////////////////////////////////
#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

class A{};

struct B{
	MULTI_NO_UNIQUE_ADDRESS A x;
	double y;
};

int main(){
#if not defined(__INTEL_COMPILER) and not defined(__NVCC__)
	static_assert( sizeof(B) == sizeof(double) , "!");
#endif
	B b;
	double& by = b.y; (void)by;
}
#endif
#endif

