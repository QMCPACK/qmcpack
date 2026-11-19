// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_EXTENSIONS_HPP
#define BOOST_MULTI_DETAIL_EXTENSIONS_HPP
// #pragma once

#include "boost/multi/detail/index_range.hpp"

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && !defined(_MSC_VER) && __has_include(<format>)
#include <format>
#endif

#include <tuple>

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost::multi::detail {

template<class... Exts>
class extensions;

template<>
class extensions<> {
 public:
	extensions() = default;
};

template<class Ex, class... Exts>
class extensions<Ex, Exts...> : private Ex {
	extensions<Exts...> rest_;

 public:
	extensions() = default;

	BOOST_MULTI_HD explicit constexpr extensions(Ex ex, Exts... rest) : Ex{ex}, rest_{rest...} {}

	template<std::size_t I>
	BOOST_MULTI_HD constexpr auto get() const {
		if constexpr(I == 0) {
			return static_cast<Ex const&>(*this);
		} else {
			return rest_.template get<I - 1>();
		}
	}

	using extension_type = Ex;
	using sub_type       = extensions<Exts...>;

	BOOST_MULTI_HD constexpr auto extension() const { return static_cast<Ex const&>(*this); }
	[[nodiscard]] BOOST_MULTI_HD constexpr auto extent() const { return static_cast<Ex const&>(*this); }

	BOOST_MULTI_HD constexpr auto sub() const { return rest_; }
};

template<class... Exts> extensions(Exts...) -> extensions<Exts...>;

template<std::size_t I, class... Exts>
auto get(::boost::multi::detail::extensions<Exts...> const& exts)
	-> decltype(exts.template get<I>()) {
	return exts.template get<I>();
}

template<class Self>
struct tyid {
	using type = Self;
};

}  // end namespace boost::multi::detail

template<class... Exts>
struct std::tuple_size<::boost::multi::detail::extensions<Exts...>> {  // NOLINT(cert-dcl58-cpp) structured binding
	static constexpr std::size_t value = sizeof...(Exts);
};

template<std::size_t I, class... Exts>
struct std::tuple_element<I, ::boost::multi::detail::extensions<Exts...>> {  // NOLINT(cert-dcl58-cpp) structured binding
	using type = typename std::conditional_t<  // NOLINT(modernize-type-traits) bug in clang-tidy
		I == 0,
		::boost::multi::detail::tyid<
			typename ::boost::multi::detail::extensions<Exts...>::extension_type
		>,
		::std::tuple_element<
			I - 1, typename ::boost::multi::detail::extensions<Exts...>::sub_type
		>
	>::type;
};

#if defined(__cpp_lib_format) && (__cpp_lib_format >= 202106L) && !defined(_MSC_VER)

template<class... Exts>
struct std::formatter<::boost::multi::detail::extensions<Exts...> > {  // NOLINT(cert-dcl58-cpp) it's the way
	constexpr auto parse(std::format_parse_context& /*ctx*/) { return /* */; }

	auto format(::boost::multi::detail::extensions<Exts...> const& obj, std::format_context& ctx) const {
		using std::get;
		if constexpr(sizeof...(Exts) == 1) {
			return std::format_to(ctx.out(), "({})", get<0>(obj));
		}
		if constexpr(sizeof...(Exts) == 2) {
			return std::format_to(ctx.out(), "({} x {})", get<0>(obj), get<1>(obj));
		}
		if constexpr(sizeof...(Exts) == 3) {
			return std::format_to(ctx.out(), "({} x {})", get<0>(obj), get<1>(obj), get<2>(obj));
		}
		if constexpr(sizeof...(Exts) == 4) {
			return std::format_to(ctx.out(), "({} x {} x {} x {})", get<0>(obj), get<1>(obj), get<2>(obj), get<3>(obj));
		}
	}
};

#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_EXTENSIONS_HPP
