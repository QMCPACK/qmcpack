// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #pragma once
#ifndef BOOST_MULTI_DETAIL_OUTER_HPP
#define BOOST_MULTI_DETAIL_OUTER_HPP

#include "boost/multi/detail/index_range.hpp"  // IWYU pragma: export  // for index_extension, extension_t, tuple, intersection, range, operator!=, operator==
#include "boost/multi/detail/interfaces.hpp"
#include "boost/multi/detail/operators.hpp"      // IWYU pragma: export  // for equality_comparable
#include "boost/multi/detail/serialization.hpp"  // IWYU pragma: export  // for archive_traits
#include "boost/multi/detail/tuple_zip.hpp"      // IWYU pragma: export  // for get, tuple, tuple_prepend, tail, tuple_prepend_t, ht_tuple
#include "boost/multi/detail/types.hpp"          // IWYU pragma: export  // for dimensionality_type, index, size_type, difference_type, size_t

// #include <algorithm>         // for max
// #include <array>             // for array
// #include <cassert>           // for assert

// #ifdef __HIP_PLATFORM_AMD__
// #include <hip/hip_runtime.h>  // it seems that AMD, HIP, ROCM 6.4, clang 21 needs this to have a working assert in host device functions
// #endif

// #include <cstddef>           // for size_t, ptrdiff_t, __GLIBCXX__
// #include <cstdlib>           // for abs
// #include <initializer_list>  // for initializer_list
// #include <iterator>
// #include <memory>       // for swap
#include <tuple>  // for tuple_element, tuple, tuple_size, tie, make_index_sequence, index_sequence
// #include <type_traits>  // for enable_if_t, integral_constant, decay_t, declval, make_signed_t, common_type_t
// #include <utility>      // for forward

// #if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
// #if !defined(__clang_major__) || !(__clang_major__ == 16)
// #include <ranges>    // IWYU pragma: keep
// #endif
// #endif

// // clang-format off
// namespace boost::multi { template <boost::multi::dimensionality_type D, typename SSize = multi::size_type> struct layout_t; }
// namespace boost::multi::detail { template <class ...Ts> class tuple; }
// // clang-format on

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

// #ifdef _MSC_VER
// #pragma warning(push)
// #pragma warning(disable : 4514)  // inline function removed, in MSVC C++17 mode
// #pragma warning(disable : 5045)  // Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified
// #endif

namespace boost::multi {

// template<class... Args> using tuple = ::std::tuple<Args...>;  // TODO(correaa) use libcuda++ in the future

namespace stdx {
template<class Tuple> auto tail(Tuple const& tup) noexcept {
	return std::apply([](auto const& /*head*/, auto const&... rest) noexcept -> auto { return std::make_tuple(rest...); }, tup);
}
template<class Tuple> auto head(Tuple const& tup) noexcept {
	return std::apply([](auto const& head, auto const&... /*rest*/) noexcept -> decltype(auto) { return head; }, tup);
}

template<class T, class Tuple> auto append_front(T const& value, Tuple const& tup) noexcept {
	return std::apply([&](auto const&... rest) noexcept -> auto { return std::make_tuple(value, rest...); }, tup);
}

template<class Tuple> struct tuple_tail;
template<class T, class... Ts> struct tuple_tail<std::tuple<T, Ts...>> {
	using type = std::tuple<Ts...>;
};

template<class Tuple> using tuple_type_t = typename tuple_tail<Tuple>::type;
}  // end namespace stdx

namespace detail {

template<class... Exts>
class outer_t;

template<class Outer>
class outer_elements_t;  // lazy random-access range over the flattened coordinate tuples (defined below)

template<>
class outer_t<> : public std::tuple<> {  // base case: zero-dimensional extents
 public:
	static constexpr dimensionality_type dimensionality = 0;

	using index        = multi::index;
	using indices_type = std::tuple<>;
	using element      = std::tuple<>;

	outer_t() = default;

	constexpr auto num_elements() const noexcept -> multi::ssize_t { return 1; }  // NOLINT(readability-convert-member-functions-to-static)

	constexpr auto from_linear(index /*n*/) const noexcept -> indices_type { return {}; }  // NOLINT(readability-convert-member-functions-to-static)
	constexpr auto sizes() const noexcept -> std::tuple<> { return {}; }                    // NOLINT(readability-convert-member-functions-to-static)

	friend constexpr auto operator==(outer_t const& /*self*/, outer_t const& /*other*/) noexcept -> bool { return true; }
	friend constexpr auto operator!=(outer_t const& /*self*/, outer_t const& /*other*/) noexcept -> bool { return false; }
};

// A lazy, random-access, sized range over the flattened coordinate tuples of an `outer_t`:
// `outer.elements()[n] == outer.from_linear(n)`. This is what makes `outer` behave like an
// array of (coordinate) tuples -- the cartesian product materialized lazily, element by element.
template<class Outer>
class outer_elements_t {
	Outer outer_;

 public:
	using value_type      = typename Outer::indices_type;
	using size_type       = multi::ssize_t;
	using difference_type = std::ptrdiff_t;

	explicit constexpr outer_elements_t(Outer outer) : outer_{std::move(outer)} {}

	constexpr auto size() const -> size_type { return outer_.num_elements(); }
	constexpr auto operator[](difference_type n) const -> value_type { return outer_.from_linear(static_cast<typename Outer::index>(n)); }

	class iterator {
		Outer          outer_;
		std::ptrdiff_t n_ = 0;

	 public:
		using value_type        = typename Outer::indices_type;
		using reference         = value_type;  // proxy: dereference yields a fresh coordinate tuple by value
		using pointer           = void;
		using difference_type   = std::ptrdiff_t;
		using iterator_category = std::random_access_iterator_tag;

		iterator() = default;
		constexpr iterator(Outer outer, difference_type n) : outer_{std::move(outer)}, n_{n} {}

		constexpr auto operator*() const -> value_type { return outer_.from_linear(static_cast<typename Outer::index>(n_)); }
		constexpr auto operator[](difference_type d) const -> value_type { return *(*this + d); }

		constexpr auto operator++() -> iterator& { ++n_; return *this; }
		constexpr auto operator--() -> iterator& { --n_; return *this; }
		constexpr auto operator++(int) -> iterator { auto ret = *this; ++n_; return ret; }  // NOLINT(cert-dcl21-cpp)
		constexpr auto operator--(int) -> iterator { auto ret = *this; --n_; return ret; }  // NOLINT(cert-dcl21-cpp)

		constexpr auto operator+=(difference_type d) -> iterator& { n_ += d; return *this; }
		constexpr auto operator-=(difference_type d) -> iterator& { n_ -= d; return *this; }

		friend constexpr auto operator+(iterator it, difference_type d) -> iterator { it += d; return it; }
		friend constexpr auto operator+(difference_type d, iterator it) -> iterator { it += d; return it; }
		friend constexpr auto operator-(iterator it, difference_type d) -> iterator { it -= d; return it; }
		friend constexpr auto operator-(iterator const& lhs, iterator const& rhs) -> difference_type { return lhs.n_ - rhs.n_; }

		friend constexpr auto operator==(iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ == rhs.n_; }
		friend constexpr auto operator!=(iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ != rhs.n_; }
		friend constexpr auto operator< (iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ <  rhs.n_; }
		friend constexpr auto operator> (iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ >  rhs.n_; }
		friend constexpr auto operator<=(iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ <= rhs.n_; }
		friend constexpr auto operator>=(iterator const& lhs, iterator const& rhs) -> bool { return lhs.n_ >= rhs.n_; }
	};

	constexpr auto begin() const -> iterator { return iterator{outer_, 0}; }
	constexpr auto end()   const -> iterator { return iterator{outer_, size()}; }
};

template<class Ext, class... Exts>
class outer_t<Ext, Exts...> : public std::tuple<Ext, Exts...> {  // TODO(correaa) use libcuda++ in the future https://github.com/boostorg/math/blob/develop/include/boost/math/tools/cstdint.hpp
	using base_ = std::tuple<Ext, Exts...>;

 public:
	static constexpr dimensionality_type dimensionality = sizeof...(Exts);
	// static constexpr static dimensionality_type rank_v = sizeof...(Exts);

	using element = std::tuple<Ext, Exts...>;
	// using element = std::tuple<typename Exts::value_type...>;

	using extension_type  = Ext;  // std::tuple_element_t<0, std::tuple<Exts...>>;
	using difference_type = typename extension_type::difference_type;
	using size_type       = typename extension_type::size_type;
	using index           = typename extension_type::index;

	outer_t() = default;
	// using std::tuple<Exts...>::tuple;
	explicit outer_t(Ext ext, Exts... xs) : std::tuple<Ext, Exts...>{std::move(ext), std::move(xs)...} {}

	auto extension() const -> extension_type {
		using std::get;
		return get<0>(static_cast<std::tuple<Ext, Exts...> const&>(*this));
	}
	[[nodiscard]] constexpr auto extent() const -> extension_type {
		using std::get;
		return get<0>(static_cast<std::tuple<Ext, Exts...> const&>(*this));
	}

	auto size() const -> size_type { return extension().size(); }

	auto num_elements() const { return size() * sub().num_elements(); }

	using sub_type = outer_t<Exts...>;

	constexpr auto sub() const -> sub_type {
		return std::apply(
			[](auto... xs) noexcept -> sub_type { return sub_type(xs...); },
			stdx::tail(static_cast<std::tuple<Ext, Exts...> const&>(*this))
		);
	}

	// the coordinate type produced when fully indexing: a tuple of `dimensionality + 1` indices
	using indices_type = decltype(std::tuple_cat(std::declval<std::tuple<index>>(), std::declval<typename sub_type::indices_type>()));

	// row-major decomposition of a linear index `n` into the full coordinate tuple
	constexpr auto from_linear(index n) const -> indices_type {
		auto const sub_num = sub().num_elements();
		return std::tuple_cat(
			std::make_tuple(static_cast<index>(extension().first() + (n / sub_num))),
			sub().from_linear(static_cast<index>(n % sub_num))
		);
	}

	// the structured factors of the cartesian product (the D extensions, as a tuple)
	constexpr auto extensions() const -> element { return static_cast<base_ const&>(*this); }

	// the per-dimension sizes
	constexpr auto sizes() const { return std::tuple_cat(std::make_tuple(size()), sub().sizes()); }

	// lazy random-access array of the `num_elements()` coordinate tuples (cartesian product, flattened)
	constexpr auto elements() const -> outer_elements_t<outer_t> { return outer_elements_t<outer_t>{*this}; }

	class iterator : public detail::forward_iterator_facade<iterator> {
		friend class outer_t;

		typename extension_type::iterator idx_;
		outer_t<Exts...>                  rest_;

		constexpr explicit iterator(typename extension_type::iterator idx, outer_t<Exts...> rest)
		: idx_{idx}, rest_{std::move(rest)} {}

	 public:
		iterator() = default;

		constexpr auto operator++() -> iterator& {
			++idx_;
			return *this;
		}

		constexpr auto operator--() -> iterator& {
			--idx_;
			return *this;
		}

		auto operator==(iterator const& other) const -> bool {
			assert( rest_ == other.rest_ );
			return idx_ == other.idx_;
		}

		using difference_type   = typename extension_type::iterator::difference_type;
		using value_type        = void;
		using reference_type    = value_type;
		using pointer           = void;
		using iterator_category = std::random_access_iterator_tag;

		constexpr auto operator+=(difference_type d) -> iterator& { idx_ += d; return *this; }
		constexpr auto operator+(difference_type d) const {
			iterator ret{*this};  // mull-ignore: cxx_init_const
			ret += d;
			return ret;
		}
		// constexpr auto operator-(difference_type d) const { return iterator{idx_ - d, rest_}; }

		constexpr auto operator*() const {
			return stdx::append_front(*idx_, static_cast<std::tuple<Exts...> const&>(rest_));
		}

		constexpr auto operator[](index idx) const { return *((*this) + idx); }

		// 	 public:


		// 		friend constexpr auto operator-(iterator const& self, iterator const& other) -> difference_type { return self.idx_ - other.idx_; }
		// 		friend constexpr auto operator+(difference_type n, iterator const& self) { return self + n; }

		// 		constexpr auto operator+=(difference_type d) -> iterator& { idx_ += d; return *this; }
		// 		constexpr auto operator-=(difference_type d) -> iterator& { idx_ -= d; return *this; }

		// 		constexpr auto operator--() -> iterator& { --idx_; return *this; }

		// 		constexpr auto operator++(int) -> iterator { iterator ret{*this}; operator++(); return ret; }  // NOLINT(cert-dcl21-cpp)
		// 		constexpr auto operator--(int) -> iterator { iterator ret{*this}; operator--(); return ret; }  // NOLINT(cert-dcl21-cpp)


		// 		friend constexpr auto operator==(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ == other.idx_; }
		// 		friend constexpr auto operator!=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ != other.idx_; }

		// 		friend constexpr auto operator<(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ < other.idx_; }
		// 		friend constexpr auto operator>(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ > other.idx_; }

		// 		friend constexpr auto operator<=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ <= other.idx_; }
		// 		friend constexpr auto operator>=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ >= other.idx_; }
		// 	};
	};

	constexpr auto begin() const noexcept {
		return iterator{
			stdx::head(static_cast<std::tuple<Ext, Exts...> const&>(*this)).begin(),
			std::apply([](auto... xs) noexcept -> auto { return outer_t<Exts...>(xs...); }, stdx::tail(static_cast<std::tuple<Ext, Exts...> const&>(*this)))
		};
	}

	template<std::size_t I>
	friend auto get(outer_t const& self) {
		using std::get;
		return get<I>(static_cast<std::tuple<Ext, Exts...> const&>(self));
	}
};

template<class... Exts> outer_t(Exts...) -> outer_t<decltype(multi::extent_t(std::declval<Exts>()))...>;

}  // end namespace detail

}  // end namespace boost::multi

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmismatched-tags"  // some std libs declare tuple_size as a class
#endif
template<class... Exts>
struct std::tuple_size<::boost::multi::detail::outer_t<Exts...>> {  // NOLINT(cert-dcl58-cpp) structured binding
	static constexpr std::size_t value = sizeof...(Exts);
};

template<std::size_t I, class... Exts>
struct std::tuple_element<I, ::boost::multi::detail::outer_t<Exts...>> {  // NOLINT(cert-dcl58-cpp) structured binding
	using type = std::tuple_element_t<I, std::tuple<Exts...>>;
};
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_OUTER_HPP
