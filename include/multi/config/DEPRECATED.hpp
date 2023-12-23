// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#ifndef MULTI_CONFIG_DEPRECATED_HPP_
#define MULTI_CONFIG_DEPRECATED_HPP_
// Copyright 2019-2022 Alfredo A. Correa

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#define DEPRECATED(MsG) [[deprecated]]

//  #ifdef __NVCC__
//  	#define DEPRECATED(MsG) __attribute__((deprecated))
//  #else
//  	#if __has_cpp_attribute(deprecated)
//  		#define	DEPRECATED(MsG) [[deprecated(MsG)]]
//  	#else
//  		#define DEPRECATED(MsG)
//  	#endif
//  #endif

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

#endif  // MULTI_CONFIG_DEPRECATED_HPP_

