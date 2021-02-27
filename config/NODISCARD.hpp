#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
echo $X
$CXXX $CXXFLAGS $0 -o $0.$X &&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_CONFIG_NODISCARD_HPP
#define MULTI_CONFIG_NODISCARD_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#ifndef NODISCARD
#if defined(__NVCC__)
	#define NODISCARD(MsG)
#elif (__has_cpp_attribute(nodiscard) and (__cplusplus>=201703L))
	#if (__has_cpp_attribute(nodiscard)>=201907) and (__cplusplus>201703L)
		#define NODISCARD(MsG) [[nodiscard(MsG)]]
	#else
		#define NODISCARD(MsG) [[nodiscard]]
	#endif
#elif __has_cpp_attribute(gnu::warn_unused_result)
	#define NODISCARD(MsG) [[gnu::warn_unused_result]]
#else
	#define NODISCARD(MsG)
#endif
#endif

#ifndef NODISCARD_CLASS
	#if(__has_cpp_attribute(nodiscard) and not defined(__NVCC__) and (not defined(__clang__) or (defined(__clang__) and (__cplusplus >= 202002L))))
		#if (__has_cpp_attribute(nodiscard)>=201907)
			#define NODISCARD_CLASS(MsG) [[nodiscard_(MsG)]]
		#else
			#define NODISCARD_CLASS(MsG) [[nodiscard]]
		#endif
	#else
		#define NODISCARD_CLASS(MsG)
	#endif
#endif

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include "../config/MAYBE_UNUSED.hpp"

NODISCARD("because...") int f(){return 5;}

struct A{
	NODISCARD("because...")
	friend
	int ff(){return 5.;}
};

struct NODISCARD_CLASS("because...") B{};

B create_B(){return B{};}

int main(){
	int i; 
	i = f(); // ok
//	f();  // warning
	++i;
	(void)i;
	
	MAYBE_UNUSED auto b = create_B();
}
#endif
#endif


