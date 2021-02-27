#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_CONFIG_MAYBE_UNUSED_HPP
#define MULTI_CONFIG_MAYBE_UNUSED_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if (__has_cpp_attribute(maybe_unused)) and (__cplusplus>=201703L)
	#define MULTI_MAYBE_UNUSED [[maybe_unused]]
#elif __has_cpp_attribute(gnu::unused)
	#define MULTI_MAYBE_UNUSED [[gnu::unused]]
#elif __has_cpp_attribute(__attribute__((unused)))
	#define MULTI_MAYBE_UNUSED __attribute__((unused))
#else
	#define MULTI_MAYBE_UNUSED
#endif

#ifndef MAYBE_UNUSED
#define MAYBE_UNUSED MULTI_MAYBE_UNUSED
#endif

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_CONFIG_NODISCARD

int main(){
	MAYBE_UNUSED int i; 
}
#endif
#endif


