// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MULTI_CONFIG_MAYBE_UNUSED_HPP
#define MULTI_CONFIG_MAYBE_UNUSED_HPP

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if (__has_cpp_attribute(maybe_unused)) and (__cplusplus>=201703L)
	#define MULTI_MAYBE_UNUSED [[maybe_unused]]  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is really necessary in C++17
#elif __has_cpp_attribute(gnu::unused)
	#define MULTI_MAYBE_UNUSED [[gnu::unused]]  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is really necessary in C++17
#elif __has_cpp_attribute(__attribute__((unused)))
	#define MULTI_MAYBE_UNUSED __attribute__((unused))  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is really necessary in C++17
#else
	#define MULTI_MAYBE_UNUSED  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is really necessary in C++17
#endif

#ifndef MAYBE_UNUSED
	#define MAYBE_UNUSED MULTI_MAYBE_UNUSED  // NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) check if this is really necessary in C++17
#endif

#endif

