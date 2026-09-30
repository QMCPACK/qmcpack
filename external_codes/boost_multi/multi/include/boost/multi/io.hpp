// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_IO_HPP
#define BOOST_MULTI_IO_HPP
// #pragma once

#include "boost/multi/utility.hpp"

#include <iostream>

// #if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
// #if __has_include(<format>)
// #include <boost/multi/io/format.hpp>
// #endif
// #endif

namespace boost::multi {

namespace detail {

template<class Array, std::enable_if_t<!has_dimensionality<Array>::value, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
void print(std::ostream& os, Array const& arr, std::string_view /*open*/, std::string_view /*sep*/, std::string_view /*close*/, std::string_view /*tag*/, int /*indent*/) {
	os << arr;
}

template<class Array, std::enable_if_t<has_dimensionality<Array>::value && Array::dimensionality == 0, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
void print(std::ostream& os, Array const& arr, std::string_view /*open*/, std::string_view /*sep*/, std::string_view /*close*/, std::string_view /*tag*/, int /*indent*/) {
	assert(!arr.empty());
	os << static_cast<typename Array::element_cref>(arr);
}

template<class Array, std::enable_if_t<has_dimensionality<Array>::value && (Array::dimensionality > 0), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
void print(std::ostream& os, Array const& arr, std::string_view open, std::string_view sep, std::string_view close, std::string_view tab, int indent) {
	for(auto count = 0; count != indent; ++count) {  // NOLINT(altera-unroll-loops)
		os << tab;
	}
	os << open[0];
	if constexpr(Array::dimensionality > 1) {
		os << '\n';
	}
	for(auto idx : arr.extent()) {  // NOLINT(altera-unroll-loops) TODO(correaa) use an algorithm
		multi::detail::print(os, arr[idx], open.size() == 1 ? open : open.substr(1), sep.size() == 1 ? sep : sep.substr(1), close.size() == 1 ? close : close.substr(1), tab.size() == 1 ? tab : tab.substr(1), indent + 1);
		if(idx != arr.extent().back()) {
			os << sep[0];
			if constexpr(Array::dimensionality > 1) {
				os << '\n';
			} else {
				os << ' ';
			}
		}
	}

	if constexpr(Array::dimensionality > 1) {
		os << sep[0] << ' ' << '\n';
		for(auto count = 0; count != indent; ++count) {  // NOLINT(altera-unroll-loops) TODO(correaa) use an algorithm
			os << tab;
		}
	}

	os << close[0];
}

}  // namespace detail

template<class Array, std::enable_if_t<Array::dimensionality >= 0, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
auto operator<<(std::ostream& os, Array const& arr) -> std::ostream& {
	multi::detail::print(os, arr, "{", ",", "}", "\t", 0);
	return os;
}

template<typename Integer>
auto operator<<(std::ostream& os, extent_t<Integer> const& ext) -> std::ostream& {
	if(ext.empty()) {
		return os << "[)";
	}
	// if(ext.front() != 0) {
	// 	return os << "[" << ext.front() << ", " << ext.back() + 1 << ")";
	// }
	return os << "[" << ext.front() << ", " << ext.back() + 1 << ")";
}

}  // namespace boost::multi

#endif  // BOOST_MULTI_IO_HPP
