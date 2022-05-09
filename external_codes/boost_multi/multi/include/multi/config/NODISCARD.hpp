// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

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
		#define NODISCARD(MsG) [[nodiscard]]  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is needed in C++17
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

#endif
