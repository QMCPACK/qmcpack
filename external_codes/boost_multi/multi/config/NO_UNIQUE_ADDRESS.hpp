#ifdef COMPILATION_INSTRUCTIONS
echo $CXXX $CXXFLAGS
$CXXX $CXXFLAGS $0 -o $0.$X&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP
#define MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if __has_cpp_attribute(no_unique_address) >=201803 and not defined(__NVCC__)
	#define MULTI_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
	#define MULTI_NO_UNIQUE_ADDRESS
#endif

////////////////////////////////////////////////////////////////////////////////
#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__ // _TEST_MULTI_CONFIG_NO_UNIQUE_ADDRESS

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

