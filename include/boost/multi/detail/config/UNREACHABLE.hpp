// Copyright 2024-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_CONFIG_UNREACHABLE_HPP
#define BOOST_MULTI_DETAIL_CONFIG_UNREACHABLE_HPP

// BOOST_MULTI_UNREACHABLE(): marks a point the control flow can never reach.
// Prefer the standardized `std::unreachable()` (C++23), fall back to compiler
// intrinsics, and finally to a hard trap so the fallback stays defined behavior.

#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif

#if defined(__cpp_lib_unreachable) && (__cpp_lib_unreachable >= 202202L)
#  include <utility>
#  define BOOST_MULTI_UNREACHABLE() ::std::unreachable()
#elif defined(__GNUC__) || defined(__clang__)  // includes nvcc/clang-cuda host paths
#  define BOOST_MULTI_UNREACHABLE() __builtin_unreachable()
#elif defined(__EDG__)  // e.g. nvcc's front end in MSVC-host mode
#  define BOOST_MULTI_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#  define BOOST_MULTI_UNREACHABLE() __assume(false)
#else
#  include <cstdlib>
#  define BOOST_MULTI_UNREACHABLE() ::std::abort()
#endif

#endif  // BOOST_MULTI_DETAIL_CONFIG_UNREACHABLE_HPP
