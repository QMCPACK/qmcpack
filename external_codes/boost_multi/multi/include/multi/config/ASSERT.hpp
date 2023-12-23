// Copyright 2019-2023 Alfredo A. Correa

#ifndef MULTI_CONFIG_ASSERT_HPP_
#define MULTI_CONFIG_ASSERT_HPP_

#include<cassert>

#if defined(MULTI_ACCESS_NDEBUG) or defined(__CUDACC__)
	#define MULTI_ACCESS_ASSERT(Expr)  // NOLINT(cppcoreguidelines-macro-usage
#else
	// #include<stacktrace>
	// // NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this is for very inefficient asserts
	// #if defined(__cpp_lib_stacktrace) && (__cpp_lib_stacktrace >= 202011L)
	// #define MULTI_ACCESS_ASSERT(Expr) assert((std::cerr<<std::stacktrace()<<std::endl) && (Expr))
	// #else
	#define MULTI_ACCESS_ASSERT(Expr) assert(Expr)  // NOLINT(cppcoreguidelines-macro-usage)
	// #endif
#endif

#endif  // MULTI_CONFIG_ASSERT_HPP_
