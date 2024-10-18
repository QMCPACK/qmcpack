// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_CONFIG_ASSERT_HPP
#define BOOST_MULTI_DETAIL_CONFIG_ASSERT_HPP

#include<cassert>

#if defined(BOOST_MULTI_ACCESS_NDEBUG) || defined(__CUDACC__)
	#define BOOST_MULTI_ACCESS_ASSERT(Expr)  // NOLINT(cppcoreguidelines-macro-usage
#else
	// #include<stacktrace>
	// // NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this is for very inefficient asserts
	// #if defined(__cpp_lib_stacktrace) && (__cpp_lib_stacktrace >= 202011L)
	// #define BOOST_MULTI_ACCESS_ASSERT(Expr) assert((std::cerr<<std::stacktrace()<<std::endl) && (Expr))
	// #else
	#define BOOST_MULTI_ACCESS_ASSERT(Expr) assert(Expr)  // NOLINT(cppcoreguidelines-macro-usage)
	// #endif
#endif

#endif  // BOOST_MULTI_DETAIL_CONFIG_ASSERT_HPP
