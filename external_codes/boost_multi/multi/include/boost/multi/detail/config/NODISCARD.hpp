// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_CONFIG_NODISCARD_HPP
#define BOOST_MULTI_DETAIL_CONFIG_NODISCARD_HPP

// clang-format off
#ifdef __has_cpp_attribute

// No discard first
#  ifdef __NVCC__
#    define BOOST_MULTI_NODISCARD(MsG)
#  elif __has_cpp_attribute(nodiscard) && (__cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L))
#    if (__has_cpp_attribute(nodiscard) >= 201907L) && (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
#      define BOOST_MULTI_NODISCARD(MsG) [[nodiscard]]  // [[nodiscard(MsG)]] in c++20 empty message is not allowed with paren
#    else
#      define BOOST_MULTI_NODISCARD(MsG) [[nodiscard]]  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is needed in C++17
#    endif
#  elif __has_cpp_attribute(gnu::warn_unused_result)
#    define BOOST_MULTI_NODISCARD(MsG) [[gnu::warn_unused_result]]
#  endif 

// No discard class
#  if(__has_cpp_attribute(nodiscard) && !defined(__NVCC__) && (!defined(__clang__) || (defined(__clang__) && (__cplusplus >= 202002L)))) && (__cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L))
#    if (__has_cpp_attribute(nodiscard) >= 201907L) && (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
#      define BOOST_MULTI_NODISCARD_CLASS(MsG) [[nodiscard_(MsG)]]
#    else
#      define BOOST_MULTI_NODISCARD_CLASS(MsG) [[nodiscard]]
#    endif
#  endif

#endif

#ifndef BOOST_MULTI_NODISCARD
#  define BOOST_MULTI_NODISCARD(MsG)
#endif

#ifndef BOOST_MULTI_NODISCARD_CLASS
#  define BOOST_MULTI_NODISCARD_CLASS(MsG)
#endif

// clang-format on

#endif  // BOOST_MULTI_DETAIL_CONFIG_NODISCARD_HPP
