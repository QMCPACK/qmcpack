// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_CONFIG_NO_UNIQUE_ADDRESS_HPP
#define BOOST_MULTI_DETAIL_CONFIG_NO_UNIQUE_ADDRESS_HPP

#ifdef __has_cpp_attribute
	#if __has_cpp_attribute(no_unique_address) >= 201803L && !defined(__NVCC__) && !defined(__PGI) && (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
		// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this macro will be needed until C++20
		#define BOOST_MULTI_NO_UNIQUE_ADDRESS [[no_unique_address]]
	#endif
#endif

#ifndef BOOST_MULTI_NO_UNIQUE_ADDRESS
	#ifdef _MSC_VER
		#define BOOST_MULTI_NO_UNIQUE_ADDRESS  // [[msvc::no_unique_address]]
	#else
		#define BOOST_MULTI_NO_UNIQUE_ADDRESS
	#endif
#endif

#endif  // BOOST_MULTI_DETAIL_CONFIG_NO_UNIQUE_ADDRESS_HPP
