// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_IO_FORMAT_HPP
#define BOOST_MULTI_IO_FORMAT_HPP
// #pragma once

#include <boost/multi/utility.hpp>

#include <iostream>
#if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
#if __has_include(<format>)
#include <format>
#endif
#endif

#if defined(__cpp_lib_format) && (__cpp_lib_format >= 202110L)
template<typename T, boost::multi::dimensionality_type D>
struct std::formatter<boost::multi::array<T, D>, char> {  // NOLINT(cert-dcl58-cpp) this is the way

	constexpr auto parse(std::format_parse_context& ctx) {
		return ctx.begin();
	}

	template<class FormatContext>
	auto format(boost::multi::array<T, D> const& arr, FormatContext& ctx) const {
		auto out = ctx.out();
		*out++   = '[';

		bool first = true;
		for(auto const& elem : arr) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
			if(!first) {
				*out++ = ',';
				*out++ = ' ';
			}
			out   = std::format_to(out, "{}", elem);
			first = false;
		}

		*out++ = ']';
		return out;
	}
};

template<typename T, boost::multi::dimensionality_type D, typename P, class L>
struct std::formatter<::boost::multi::const_subarray<T, D, P, L>, char> {  // NOLINT(cert-dcl58-cpp) this is the way

	constexpr auto parse(std::format_parse_context& ctx) {
		return ctx.begin();
	}

	template<class FormatContext>
	auto format(::boost::multi::const_subarray<T, D, P, L> const& arr, FormatContext& ctx) const {
		auto out = ctx.out();
		*out++   = '[';

		bool first = true;
		for(auto const& elem : arr) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
			if(!first) {
				*out++ = ',';
				*out++ = ' ';
			}
			out   = std::format_to(out, "{}", elem);
			first = false;
		}

		*out++ = ']';
		return out;
	}
};

#endif

#endif  // BOOST_MULTI_IO_FORMAT_HPP
