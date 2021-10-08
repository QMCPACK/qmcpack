#ifndef MULTI_CONFIG_MAYBE_UNUSED_HPP // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#define MULTI_CONFIG_MAYBE_UNUSED_HPP
// © Alfredo A. Correa 2020-2021

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(name) 0
#endif

#if (__has_cpp_attribute(maybe_unused)) and (__cplusplus>=201703L)
	#define MULTI_MAYBE_UNUSED [[maybe_unused]]
#elif __has_cpp_attribute(gnu::unused)
	#define MULTI_MAYBE_UNUSED [[gnu::unused]]
#elif __has_cpp_attribute(__attribute__((unused)))
	#define MULTI_MAYBE_UNUSED __attribute__((unused))
#else
	#define MULTI_MAYBE_UNUSED
#endif

#ifndef MAYBE_UNUSED
#define MAYBE_UNUSED MULTI_MAYBE_UNUSED
#endif

// TODO(correaa): add MAYBE_UNUSED test
// TODO(correaa): remove header when library is C++17 or greater only
#endif

