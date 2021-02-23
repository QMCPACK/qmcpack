#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS $0 -o $0.$X &&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MPI3_CONFIG_NODISCARD_HPP
#define MPI3_CONFIG_NODISCARD_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#ifndef NODISCARD
#define MPI3_GET_MACRO(_1,NAME,...) NAME

#if   (__has_cpp_attribute(nodiscard)>=201603L) and not defined(__NVCC__)
	#define NODISCARD0()    [[nodiscard]]
	#define NODISCARD1(MsG) [[nodiscard]]
#elif (__has_cpp_attribute(nodiscard)>=201907L) and not defined(__NVCC__)
	#define NODISCARD0()    [[nodiscard]]
	#define NODISCARD1(MsG) [[nodiscard(MsG)]]
#elif (__has_cpp_attribute(warn_unused_result))
	#define NODISCARD0()    __attribute__((warn_unused_result))
	#define NODISCARD1(MsG) __attribute__((warn_unused_result))
#else
	#if __INTEL_COMPILER
		#pragma warning( push,    1028 )
		#pragma warning( disable: 1028 )
	#endif
	#pragma warning "compiler doesn't support no discard warning?"
	#if __INTEL_COMPILER
		#pragma warning( pop )
	#endif
	#define NODISCARD0()
	#define NODISCARD1(MsG)
#endif

#define NODISCARD(...) MPI3_GET_MACRO(__VA_ARGS__, NODISCARD1, NODISCARD0)(__VA_ARGS__)
#endif

#if defined(__INCLUDE_LEVEL__)  and not __INCLUDE_LEVEL__

NODISCARD("because...") int f(){return 5;}
//[[nodiscard]] int g(){return 5;} // ok in g++ -std=c++14

int main(){
	int i; 
	i = f(); // ok
//	f();  // warning
	i += 1;
}
#endif
#endif


