#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS $0 -o $0x &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_CONFIG_DEPRECATED_HPP
#define MULTI_CONFIG_DEPRECATED_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#ifdef __NVCC__
	#define DEPRECATED(MsG) __attribute__((deprecated))
#else
	#if __has_cpp_attribute(deprecated)
		#define	DEPRECATED(MsG) [[deprecated(MsG)]]
	#else
		#define DEPRECATED(MsG)
	#endif
#endif

#if not defined(__INTEL_COMPILER)
#define BEGIN_NO_DEPRECATED \
\
_Pragma("GCC diagnostic push") \
_Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"") \
\

#else
#define BEGIN_NO_DEPRECATED \
_Pragma("warning push") \
_Pragma("warning disable 1786") \

#endif

#if not defined(__INTEL_COMPILER)
#define END_NO_DEPRECATED \
\
_Pragma("GCC diagnostic pop") \
\

#else
#define END_NO_DEPRECATED \
\
_Pragma("warning pop") \
\

#endif

#define BEGIN_CUDA_SLOW BEGIN_NO_DEPRECATED
#define END_CUDA_SLOW   END_NO_DEPRECATED

#define NO_DEPRECATED(ExpR) \
	BEGIN_NO_DEPRECATED \
	ExpR \
	END_NO_DEPRECATED

#if not __INCLUDE_LEVEL__

DEPRECATED("because...") int f(){return 5;}

int main(){
	int i;
//	f(); // warning
	NO_DEPRECATED( i = f(); ) // ok
	++i;
	(void)i;
}
#endif
#endif


