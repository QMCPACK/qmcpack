// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP_
#define MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP_

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if __has_cpp_attribute(no_unique_address) >=201803 and not defined(__NVCC__) and not defined(__PGI)
	// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this macro will be needed until C++20
	#define MULTI_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
	// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this macro will be needed until C++20
	#define MULTI_NO_UNIQUE_ADDRESS
#endif

#endif  // MULTI_CONFIG_NO_UNIQUE_ADDRESS_HPP_
