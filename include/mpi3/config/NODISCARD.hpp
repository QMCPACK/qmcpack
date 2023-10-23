// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// © Alfredo A. Correa 2019-2021
#ifndef MPI3_CONFIG_NODISCARD_HPP
#define MPI3_CONFIG_NODISCARD_HPP

//#ifndef __has_cpp_attribute
//#define __has_cpp_attribute(name) 0
//#endif

//#ifndef BMPI3_NODISCARD
//#if defined(__NVCC__)
//	#define BMPI3_NODISCARD(MsG)
//#elif (__has_cpp_attribute(nodiscard) and (__cplusplus>=201703L))
//	#if (__has_cpp_attribute(nodiscard)>=201907) and (__cplusplus>201703L)
//		#define BMPI3_NODISCARD(MsG) [[nodiscard(MsG)]]
//	#else
//		#define BMPI3_NODISCARD(MsG) [[nodiscard]]
//	#endif
//#elif __has_cpp_attribute(gnu::warn_unused_result)
//	#define BMPI3_NODISCARD(MsG) [[gnu::warn_unused_result]]  // NOLINT(cppcoreguidelines-macro-usage) : replaced by `nodiscard` in C++17
//#else
//	#define BMPI3_NODISCARD(MsG)
//#endif
//#endif

#endif

