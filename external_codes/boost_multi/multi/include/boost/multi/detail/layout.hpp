// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_LAYOUT_HPP
#define BOOST_MULTI_DETAIL_LAYOUT_HPP
// #pragma once

#include "boost/multi/detail/config/NODISCARD.hpp"
#include "boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp"

#include "boost/multi/detail/index_range.hpp"    // IWYU pragma: export  // for index_extension, extension_t, tuple, intersection, range, operator!=, operator==
#include "boost/multi/detail/operators.hpp"      // IWYU pragma: export  // for equality_comparable
#include "boost/multi/detail/serialization.hpp"  // IWYU pragma: export  // for archive_traits
#include "boost/multi/detail/tuple_zip.hpp"      // IWYU pragma: export  // for get, tuple, tuple_prepend, tail, tuple_prepend_t, ht_tuple
#include "boost/multi/detail/types.hpp"          // IWYU pragma: export  // for dimensionality_type, index, size_type, difference_type, size_t

#include <algorithm>         // for max
#include <array>             // for array
#include <cassert>           // for assert

#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>  // it seems that AMD, HIP, ROCM 6.4, clang 21 needs this to have a working assert in host device functions
#endif

#include <cstddef>           // for size_t, ptrdiff_t, __GLIBCXX__
#include <cstdlib>           // for abs
#include <initializer_list>  // for initializer_list
#include <iterator>
#include <memory>       // for swap
#include <tuple>        // for tuple_element, tuple, tuple_size, tie, make_index_sequence, index_sequence
#include <type_traits>  // for enable_if_t, integral_constant, decay_t, declval, make_signed_t, common_type_t
#include <utility>      // for forward

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#if !defined(__clang_major__) || !(__clang_major__ == 16)
#include <ranges>    // IWYU pragma: keep
#endif
#endif

// clang-format off
namespace boost::multi { template <boost::multi::dimensionality_type D, typename SSize = multi::ssize_t> struct layout_t; }
namespace boost::multi::detail { template <class ...Ts> class tuple; }
// clang-format on

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4514)  // inline function removed, in MSVC C++17 mode
#pragma warning(disable : 5045)  // Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified
#endif

namespace boost::multi {

template<typename Stride>
struct stride_traits;

template<>
struct stride_traits<std::ptrdiff_t> {
	using category = std::random_access_iterator_tag;
};

template<typename Stride>
struct stride_traits {
	using category = typename Stride::category;
};

template<typename Integer>
struct stride_traits<std::integral_constant<Integer, 1>> {
#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && (!defined(__clang__) || __clang_major__ != 10)
	using category = std::contiguous_iterator_tag;
#else
	using category = std::random_access_iterator_tag;
#endif
};

namespace detail {

// template<class DeviceFun>
// class device {
// 	DeviceFun fun_;

//  public:
// 	explicit device(DeviceFun fun) : fun_{std::move(fun)} {}

// 	template<class... Args>
// #ifdef __NVCC__
// 	__device__
// #endif
// 	constexpr auto operator()(Args&&... args) const
// 	->decltype(fun_(std::forward<Args>(args)...)) {
// 		return fun_(std::forward<Args>(args)...); }
// };

template<class Tuple, std::size_t... Ns>
constexpr auto tuple_tail_impl(Tuple&& tup, std::index_sequence<Ns...> /*012*/) {
	(void)tup;  // workaround bug warning in nvcc
	using boost::multi::detail::get;
	return boost::multi::detail::tuple{std::forward<decltype(get<Ns + 1U>(std::forward<Tuple>(tup)))>(get<Ns + 1U>(std::forward<Tuple>(tup)))...};
}

template<class Tuple>
constexpr auto tuple_tail(Tuple&& t)  // NOLINT(readability-identifier-length) std naming
	-> decltype(tuple_tail_impl(std::forward<Tuple>(t), std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>> - 1U>())) {
	return tuple_tail_impl(std::forward<Tuple>(t), std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>> - 1U>());
}

}  // end namespace detail

// template<dimensionality_type D, typename SSize=multi::size_type> struct layout_t;

template<dimensionality_type D>
struct extents_t;

template<dimensionality_type D>
using sizes_t = typename extents_t<D>::sizes_type;

/// A multidimensional array value
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam Alloc Allocator type
template<typename T, dimensionality_type D, class Alloc = std::allocator<T> > struct array;  // TODO(correaa) why the declaration is in this header

/// A multidimensional array value
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam Alloc Allocator type
template<typename T, dimensionality_type D, class Alloc = std::allocator<T> > struct dynamic_array;  // TODO(correaa) why the declaration is in this header

namespace detail {
struct non_copyable_base {
	non_copyable_base(non_copyable_base const&) = delete;
	non_copyable_base(non_copyable_base&&) = default;

	non_copyable_base() = default;

	auto operator=(non_copyable_base const&) -> non_copyable_base& = default;
	auto operator=(non_copyable_base&&) -> non_copyable_base& = default;

	~non_copyable_base() = default;
};

struct copyable_base {
	copyable_base(copyable_base const&) = default;
	copyable_base(copyable_base&&) = default;

	copyable_base() = default;

	auto operator=(copyable_base const&) -> copyable_base& = default;
	auto operator=(copyable_base&&) -> copyable_base& = default;

	~copyable_base() = default;
};
}  // end namespace detail

template<dimensionality_type D>
struct extents_t : boost::multi::detail::tuple_prepend_t<index_extension, typename extents_t<D - 1>::base_> {
	using base_ = boost::multi::detail::tuple_prepend_t<index_extension, typename extents_t<D - 1>::base_>;

 public:
	static constexpr dimensionality_type dimensionality = D;
	constexpr static dimensionality_type rank_v = D;

	using difference_type = index_extension::difference_type;
	using nelems_type = multi::index;
	using size_type = index_extension::size_type;

	using element = boost::multi::detail::tuple_prepend_t<index_extension::value_type, typename extents_t<D - 1>::element>;

	extents_t() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) intentionally trivial; default-init by design

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(multi::ssize_t size1)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : allow terse syntax
	: extents_t{index_extension{size1}} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: base_{ext1} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 2, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2)
	: base_{ext1, ext2} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 3, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3)
	: base_{ext1, ext2, ext3} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 4, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4) noexcept
	: base_{ext1, ext2, ext3, ext4} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 5, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4, index_extension ext5)
	: base_{ext1, ext2, ext3, ext4, ext5} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 6, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4, index_extension ext5, index_extension ext6)
	: base_{ext1, ext2, ext3, ext4, ext5, ext6} {}

	template<class... Ts, std::enable_if_t<sizeof...(Ts) == static_cast<std::size_t>(D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; allow terse syntax // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(detail::tuple<Ts...> const& exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: extents_t(exts, std::make_index_sequence<static_cast<std::size_t>(D)>()) {}

	template<class... Ts, std::enable_if_t<sizeof...(Ts) == static_cast<std::size_t>(D), int> = 0, class = decltype(base_{std::declval<::std::tuple<Ts...> >()})>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; allow terse syntax // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(::std::tuple<Ts...> exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_{std::move(exts)} {}

	template<
		class... Exts,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
			(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
			&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
			&& !std::conjunction_v<std::is_integral<Exts>...>,  // NOLINT(modernize-type-traits) for C++20
			int> = 0
	>
	BOOST_MULTI_HD constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: base_{index_extension(exts)...} {}

	template<
		class... Exts,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
			(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
			&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
			&& std::conjunction_v<std::is_integral<Exts>...>
			&& std::conjunction_v<std::is_unsigned<Exts>...>,
			int> = 0
	>
	BOOST_MULTI_HD explicit constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: base_{index_extension(static_cast<index_extension::size_type>(exts))...} {}

	template<
		class... Exts,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
			(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
			&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
			&& std::conjunction_v<std::is_integral<Exts>...>
			&& std::conjunction_v<std::is_signed<Exts>...>,  // NOLINT(modernize-type-traits) for C++20
			int> = 0
	>
	BOOST_MULTI_HD constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: base_{index_extension(exts)...} {}

	// template<
	// 	class... Exts,
	// 	std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
	// 		(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
	// 		&& std::conjunction_v<std::is_convertible<Exts, index_extension::size_type>...>
	// 		&& std::conjunction_v<multi::detail::is_implicitly_convertible<Exts, index_extension::size_type>...>,  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
	// 		int> = 0
	// >
	// BOOST_MULTI_HD /*implicit*/ constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	// : base_{static_cast<index_extension>(static_cast<index_extension::size_type>(exts))...} {}

	// template<
	// 	class... Exts,
	// 	std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
	// 		(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
	// 		&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>
	// 		&& std::conjunction_v<multi::detail::is_implicitly_convertible<Exts, index_extension>...>,  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
	// 		int> = 0
	// >
	// BOOST_MULTI_HD /*implicit*/ constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	// : base_{static_cast<index_extension>(static_cast<typename index_extension::index>(exts))...} {}

	// template<
	// 	class... Exts,
	// 	std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
	// 		(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
	// 		&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>
	// 		&& !std::conjunction_v<multi::detail::is_implicitly_convertible<Exts, index_extension>...>,  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
	// 		int> = 0
	// >
	// BOOST_MULTI_HD explicit constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	// : base_{index_extension(exts)...} {}

	// template<class OtherExtensions,
	// 	decltype( multi::detail::implicit_cast<index_extension>(OtherExtensions{}.extent()) )* = nullptr,
	// 	decltype( multi::detail::implicit_cast<typename layout_t<D - 1>::extents_type>(OtherExtensions{}.sub()) )* = nullptr
	// >
	// // cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	// BOOST_MULTI_HD constexpr extents_t(OtherExtensions const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	// : extents_t(other.extent(), other.sub()) {}

	BOOST_MULTI_HD constexpr extents_t(index_extension const& ext, typename layout_t<D - 1>::extents_type const& other)
	: extents_t(multi::detail::ht_tuple(ext, other.base())) {}

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return *this; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return *this; }

	friend constexpr auto operator*(index_extension const& ext, extents_t const& self) -> extents_t<D + 1> {
		// return extents_t<D + 1>(tuple(extension, self.base()));
		return extents_t<D + 1>(ext, self);
	}

	friend BOOST_MULTI_HD auto operator==(extents_t const& self, extents_t const& other) { return self.base() == other.base(); }
	friend BOOST_MULTI_HD auto operator!=(extents_t const& self, extents_t const& other) { return self.base() != other.base(); }

	using index        = multi::index;
	using indices_type = multi::detail::tuple_prepend_t<index, typename extents_t<D - 1>::indices_type>;

	// template<class Func>
	// friend BOOST_MULTI_HD constexpr auto operator^(Func fun, extents_t const& xs) {
	// 	return restriction<D, Func>(xs, std::move(fun));
	// }
	// template<class Func>
	// friend constexpr auto operator->*(extents_t const& xs, Func fun) {
	// 	return restriction<D, Func>(xs, std::move(fun));
	// }

	BOOST_MULTI_HD constexpr auto sub() const {
		return extents_t<D - 1>{static_cast<base_ const&>(*this).tail()};
	}

	[[nodiscard]]
	BOOST_MULTI_HD constexpr auto from_linear(nelems_type const& n) const -> indices_type {
		auto const sub_num_elements = sub().num_elements();
		#if !(defined(__NVCC__) || defined(__HIP_PLATFORM_NVIDIA__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__))
		assert(sub_num_elements != 0);  // clang hip doesn't allow assert in host device functions
		#endif
		return multi::detail::ht_tuple(n / sub_num_elements, sub().from_linear(n % sub_num_elements));
	}

	friend constexpr auto operator%(nelems_type idx, extents_t const& exts) { return exts.from_linear(idx); }

	constexpr explicit operator bool() const { return !layout_t<D>{*this}.empty(); }

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto to_linear(index const& idx, Indices const&... rest) const {
		auto const sub_extensions = extents_t<D - 1>{this->base().tail()};
		return (idx * sub_extensions.num_elements()) + sub_extensions.to_linear(rest...);
	}

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto operator()(index idx, Indices... rest) const { return to_linear(idx, rest...); }

	template<class Before, dimensionality_type DD>
	class cursor_t {
		Before bef_;
		// missing start indices information
		template<class, dimensionality_type> friend class cursor_t;
		friend extents_t;

	 public:
		cursor_t() = default;
		explicit cursor_t(Before const& bef) : bef_{bef} {}
		
		static constexpr dimensionality_type dimensionality = DD;

		constexpr auto operator[](difference_type n) const {
			using std::apply;
			if constexpr(DD != 1) {
				return cursor_t<typename multi::layout_t<std::tuple_size_v<Before> + 1>::indexes, DD - 1> (
					apply([n] (auto... idxs) -> auto {return detail::mk_tuple(idxs..., n);}, bef_)
				);
			} else {
				return apply([n] (auto... idxs) -> auto {return detail::mk_tuple(idxs..., n);}, bef_);
			}
		}
	};

	/// Returns a cursor to the home (e.g. top-left) element
	static auto home() { return cursor_t<tuple<>, D>{}; }

	class iterator {  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) constructor does not initialize these fields: idx_
		index idx_;
		extents_t<D - 1> rest_;
		friend extents_t;
	
		constexpr iterator(index idx, extents_t<D - 1> rest) : idx_{idx}, rest_{rest} {}

	 public:
		iterator() = default;

		using difference_type = index;
		using value_type = decltype(ht_tuple(std::declval<index>(), std::declval<extents_t<D - 1>>().base()));
		using pointer = void;
		using reference = value_type;
		using iterator_category = std::random_access_iterator_tag;

		constexpr auto operator+=(difference_type n) -> iterator& { idx_ += n; return *this; }
		constexpr auto operator-=(difference_type n) -> iterator& { idx_ -= n; return *this; }

		constexpr auto operator+(difference_type n) const { return iterator{idx_ + n, rest_}; }
		constexpr auto operator-(difference_type n) const { return iterator{idx_ - n, rest_}; }

		friend constexpr auto operator-(iterator const& self, iterator const& other) -> difference_type { assert( self.rest_ == other.rest_ ); return self.idx_ - other.idx_; }

		friend constexpr auto operator+(difference_type n, iterator const& self) { return self + n; }

		constexpr auto operator++() -> auto& { ++idx_; return *this; }
		constexpr auto operator--() -> auto& { --idx_; return *this; }

		constexpr auto operator++(int) -> iterator { iterator ret{*this}; ++idx_; return ret; }
		constexpr auto operator--(int) -> iterator { iterator ret{*this}; --idx_; return ret; }

		constexpr auto operator*() const {
			// multi::detail::what(rest_);
			return ht_tuple(idx_, rest_.base());
		}

		BOOST_MULTI_HD constexpr auto operator[](difference_type const& n) const -> reference { return *((*this) + n); }

		friend constexpr auto operator==(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ == other.idx_; }
		friend constexpr auto operator!=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ != other.idx_; }

		friend constexpr auto operator<(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ < other.idx_; }
		friend constexpr auto operator>(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ > other.idx_; }

		friend constexpr auto operator<=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ <= other.idx_; }
		friend constexpr auto operator>=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ >= other.idx_; }
	};

	constexpr auto begin() const { return iterator{this->base().head().first(), this->base().tail()}; }
	constexpr auto end()   const { return iterator{this->base().head().last() , this->base().tail()}; }

	BOOST_MULTI_HD constexpr auto operator[](index idx) const {
		return static_cast<base_ const&>(*this)[idx];
	}

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto next_canonical(index& idx, Indices&... rest) const -> bool {  // NOLINT(google-runtime-references) idx is mutated
		if(extents_t<D - 1>{this->base().tail()}.next_canonical(rest...)) {
			++idx;
		}
		if(idx == this->base().head().last()) {
			idx = this->base().head().first();
			return true;
		}
		return false;
	}
	template<class... Indices>
	constexpr auto prev_canonical(index& idx, Indices&... rest) const -> bool {  // NOLINT(google-runtime-references) idx is mutated
		if(extents_t<D - 1>{this->base().tail()}.prev_canonical(rest...)) {
			--idx;
		}
		if(idx < static_cast<index>(this->base().head().first())) {
			idx = static_cast<index>(this->base().head().back());
			return true;
		}
		return false;
	}

	class elements_t {
		extents_t xs_;
		explicit constexpr elements_t(extents_t const& exts) : xs_{exts} {}

		friend struct extents_t;

	 public:
		using difference_type = extents_t::difference_type;

		class iterator {  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) TODO(correaa) investigate
			index_extension::iterator curr_;

			static_assert( std::is_default_constructible_v<index_extension::iterator> );

			typename extents_t<D - 1>::elements_t::iterator rest_it_;
			typename extents_t<D - 1>::elements_t::iterator rest_begin_;
			typename extents_t<D - 1>::elements_t::iterator rest_end_;

			BOOST_MULTI_HD constexpr iterator(
				index_extension::iterator curr,
				typename extents_t<D - 1>::elements_t::iterator rest_it,
				typename extents_t<D - 1>::elements_t::iterator rest_begin,
				typename extents_t<D - 1>::elements_t::iterator rest_end
			)
			: curr_{curr}, rest_it_{rest_it}, rest_begin_{rest_begin}, rest_end_{rest_end} {}

			friend class elements_t;

		 public:		
			using difference_type   = elements_t::difference_type;
			using value_type        = indices_type;
			using pointer           = void;
			using reference         = value_type;
			using iterator_category = std::random_access_iterator_tag;

			iterator() = default;

			template<class CUT>
			class mk_tup {
				CUT cu_;

			 public:
				constexpr explicit mk_tup(CUT current) : cu_{current} {}
				template<class... Ts>
				constexpr auto operator()(Ts... idxs) const { return detail::mk_tuple(cu_, idxs...); }
			};

			BOOST_MULTI_HD constexpr auto operator*() const {
				// printf("op* %ld ...\n", *curr_);
				using std::apply;
				return apply(mk_tup<decltype(*curr_)>{*curr_}, *rest_it_);
				// return apply([cu = *curr_] BOOST_MULTI_HD (auto... es) {return detail::mk_tuple(cu, es...);}, *rest_it_); 
			}

			BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> iterator& {
				auto len = rest_end_ - rest_begin_;
				auto off = rest_it_ - rest_begin_;
				auto tot = off + n;

				auto quo = tot / len;
				auto res = tot % len;

				if(res < 0) {
					res += len;
					--quo;
				}

				curr_ += quo;
				rest_it_ = rest_begin_ + res;

				// if(n >= 0) {
				// 	curr_ += (rest_it_ - rest_begin_ + n) / (rest_end_ - rest_begin_);
				// 	rest_it_ = rest_begin_ + ((rest_it_ - rest_begin_ + n) % (rest_end_ - rest_begin_));
				// } else {
				// 	curr_ -= (rest_end_ - rest_it_ - n) / (rest_end_ - rest_begin_);
				// 	rest_it_ = rest_end_ - ((rest_end_ - rest_it_ - n) % (rest_end_ - rest_begin_));
				// 	if(rest_it_ == rest_end_) {
				// 		rest_it_ = rest_begin_;
				// 		++curr_;
				// 	}
				// }
				return *this;
			}

			BOOST_MULTI_HD constexpr auto operator-=(difference_type n) -> iterator& {
				if(n > 0) {  // TODO(correaa) I don't know how to overcome this mutation:  // mull-ignore: cxx_gt_to_ge
					curr_ -= (rest_end_ - rest_it_ + n) / (rest_end_ - rest_begin_);
					rest_it_ = rest_end_ - ((rest_end_ - rest_it_ + n) % (rest_end_ - rest_begin_));
					if(rest_it_ == rest_end_) {
						rest_it_ = rest_begin_;
						++curr_;
					}
				} else if(n < 0) {  // TODO(correaa) I don't know how to overcome this mutation:  // mull-ignore: cxx_lt_to_le
					curr_ += (rest_it_ - rest_begin_ - n) / (rest_end_ - rest_begin_);
					rest_it_ = rest_begin_ + ((rest_it_ - rest_begin_ - n) % (rest_end_ - rest_begin_));
				}
				return *this;
			}

			friend BOOST_MULTI_HD constexpr auto operator-(iterator const& self, iterator const& other) -> difference_type {
				return ((self.curr_ - other.curr_) * (self.rest_end_ - self.rest_begin_)) + (self.rest_it_ - self.rest_begin_) - (other.rest_it_ - other.rest_begin_);
			}

			BOOST_MULTI_HD constexpr auto operator-(difference_type n) const {
				return iterator{*this} -= n;
			}

			BOOST_MULTI_HD constexpr auto operator+(difference_type n) const {
				return iterator{*this} += n;
			}
			friend BOOST_MULTI_HD constexpr auto operator+(difference_type n, iterator const& self) -> iterator { return self + n; }  // `n + it` form, required by std::random_access_iterator

			BOOST_MULTI_HD constexpr auto operator++() -> auto& {
				++rest_it_;
				if( rest_it_ == rest_end_ ) {
					rest_it_ = rest_begin_;
					++curr_;
				}
				return *this;
			}
			BOOST_MULTI_HD constexpr auto operator++(int) -> iterator { iterator ret{*this}; ++(*this); return ret; }  // NOLINT(cert-dcl21-cpp) required by std::weakly_incrementable

			BOOST_MULTI_HD constexpr auto operator--() -> auto& {
				if( rest_it_ == rest_begin_ ) {
					rest_it_ = rest_end_;
					--curr_;
				}
				--rest_it_;
				return *this;
			}
			BOOST_MULTI_HD constexpr auto operator--(int) -> iterator { iterator ret{*this}; --(*this); return ret; }  // NOLINT(cert-dcl21-cpp)

			BOOST_MULTI_HD constexpr auto operator[](difference_type n) const { return *((*this) + n); }

			friend BOOST_MULTI_HD constexpr auto operator==(iterator const& self, iterator const& other) { return (self.curr_ == other.curr_) && (self.rest_it_ == other.rest_it_); }
			friend BOOST_MULTI_HD constexpr auto operator!=(iterator const& self, iterator const& other) { return (self.curr_ != other.curr_) || (self.rest_it_ != other.rest_it_); }

			friend BOOST_MULTI_HD constexpr auto operator< (iterator const& self, iterator const& other) { return (self.curr_ <  other.curr_) || ((self.curr_ == other.curr_) && (self.rest_it_ < other.rest_it_)); }
			friend BOOST_MULTI_HD constexpr auto operator<=(iterator const& self, iterator const& other) { return (self < other) || (self == other); }
			friend BOOST_MULTI_HD constexpr auto operator> (iterator const& self, iterator const& other) { return  other <  self; }  // for std::totally_ordered
			friend BOOST_MULTI_HD constexpr auto operator>=(iterator const& self, iterator const& other) { return !(self  <  other); }
		};

		constexpr auto begin() const {
			return iterator{
				xs_.head().begin(),
				extents_t<D - 1>{xs_.tail()}.elements().begin(),
				extents_t<D - 1>{xs_.tail()}.elements().begin(),
				extents_t<D - 1>{xs_.tail()}.elements().end(),
			};
		}

		constexpr auto end() const {
			return iterator{
				xs_.head().end(),
				extents_t<D - 1>{xs_.tail()}.elements().begin(),
				extents_t<D - 1>{xs_.tail()}.elements().begin(),
				extents_t<D - 1>{xs_.tail()}.elements().end(),
			};
		}

		BOOST_MULTI_HD constexpr auto operator[](index idx) const { return begin()[idx]; }

		BOOST_MULTI_HD constexpr auto  size() const noexcept { return xs_.num_elements(); }
		BOOST_MULTI_HD constexpr auto ssize() const noexcept { return this->size(); }
		BOOST_MULTI_HD constexpr auto usize() const noexcept { return static_cast<std::size_t>(xs_.num_elements()); }
	};

	constexpr auto elements() const { return elements_t{*this}; }

	template<class Func>
	BOOST_MULTI_HD constexpr auto element_transformed(Func fun) const { return [fun](auto const&... idxs) -> decltype(auto) { return fun(detail::mk_tuple(idxs...)); } ^(*this); }

	BOOST_MULTI_HD constexpr auto               extension() const { return this->get<0>(); }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0
	[[nodiscard]] BOOST_MULTI_HD constexpr auto extent() const { return this->get<0>(); }     // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0

	BOOST_MULTI_HD constexpr auto size() const noexcept { return this->get<0>().size(); }
	BOOST_MULTI_HD constexpr auto sizes() const {
		return this->apply([](auto const&... exts) -> auto { return multi::detail::mk_tuple(exts.size()...); });
	}

	constexpr auto rotate() const {
		this->apply([](auto const& head, auto const&... rest) -> extents_t { return extents_t(rest..., head); });
	}

	constexpr auto unrotate() const {
		this->apply([](auto const&... rest, auto const& tail) -> extents_t { return extents_t(tail, rest...); });
	}

	constexpr auto transpose() const {
		return this->apply([](auto const& head1, auto const& head2, auto const&... rest) -> extents_t { return extents_t(head2, head1, rest...); });
	}

	[[deprecated]] BOOST_MULTI_HD constexpr auto extensions() const {
		using std::apply;
		return apply([](auto... sizes) -> extents_t { return extents_t(sizes...); }, sizes());
	}
	BOOST_MULTI_HD constexpr auto extents() const {
		using std::apply;
		return apply([](auto... sizes) -> extents_t { return extents_t(sizes...); }, sizes());
	}

	using sizes_type = boost::multi::detail::tuple_prepend_t<ssize_t, typename extents_t<D - 1>::sizes_type>;

 private:
	template<class Archive, std::size_t... I>
	void serialize_impl_(Archive& arxiv, std::index_sequence<I...> /*unused012*/) {
		using boost::multi::detail::get;
		(void)std::initializer_list<unsigned>{(arxiv & multi::archive_traits<Archive>::make_nvp("extent", get<I>(this->base())), 0U)...};
	}

 public:
	template<class Archive>
	void serialize(Archive& arxiv, unsigned int const /*version*/) {
		serialize_impl_(arxiv, std::make_index_sequence<static_cast<std::size_t>(D)>());
	}

 private:
	template<class Array, std::size_t... I, typename = decltype(base_{boost::multi::detail::get<I>(std::declval<Array const&>())...})>
	BOOST_MULTI_HD constexpr extents_t(Array const& tup, std::index_sequence<I...> /*unused012*/)
	: base_{boost::multi::detail::get<I>(tup)...} {}

	static BOOST_MULTI_HD constexpr auto multiply_fold_() -> multi::ssize_t { return static_cast<multi::ssize_t>(1U); }
	static BOOST_MULTI_HD constexpr auto multiply_fold_(multi::ssize_t const& size) -> multi::ssize_t { return size; }
	template<class... As>
	static BOOST_MULTI_HD constexpr auto multiply_fold_(multi::ssize_t const& size, As const&... rest) -> multi::ssize_t { return size * static_cast<multi::ssize_t>(multiply_fold_(rest...)); }

	template<std::size_t... I>
	BOOST_MULTI_HD constexpr auto num_elements_impl_(std::index_sequence<I...> /*unused012*/) const -> multi::ssize_t {
		using boost::multi::detail::get;
		return static_cast<multi::ssize_t>(multiply_fold_(static_cast<multi::ssize_t>(get<I>(this->base()).size())...));
	}

 public:
	BOOST_MULTI_HD constexpr auto num_elements() const -> multi::ssize_t {
		return static_cast<multi::ssize_t>(num_elements_impl_(std::make_index_sequence<static_cast<std::size_t>(D)>()));
	}

	friend constexpr auto intersection(extents_t const& self, extents_t const& other) -> extents_t {
		using boost::multi::detail::get;
		return extents_t(
			multi::detail::ht_tuple(
				index_extension(intersection(get<0>(self.base()), get<0>(other.base()))),
				intersection(
					extents_t<D - 1>(self.base().tail()),
					extents_t<D - 1>(other.base().tail())
				).base()
			)
		);
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend constexpr auto get(extents_t const& self) -> typename std::tuple_element_t<Index, base_> {
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	constexpr auto get() const -> std::tuple_element_t<Index, base_> {
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}

	template<class F>
	constexpr auto apply(F&& fun) const -> decltype(auto) {
		return std::apply(std::forward<F>(fun), this->base());
	}
};

template<> struct extents_t<0> : tuple<> {
	using base_ = tuple<>;

 private:
	// base_ impl_;

 public:
	static constexpr dimensionality_type dimensionality = 0;  // TODO(correaa): consider deprecation

	using rank = std::integral_constant<dimensionality_type, 0>;
	using element = tuple<>;

	using index = multi::index;

	using nelems_type = index;
	using difference_type = index;
	using size_type = index_extension::size_type;  // TODO(correaa) or void?

	explicit BOOST_MULTI_HD constexpr extents_t(tuple<> const& tup)
	: base_{tup} {}

	extents_t() = default;

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return *this; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return *this; }

	template<class Archive> static void serialize(Archive& /*ar*/, unsigned /*version*/) { /*noop*/ }

	static BOOST_MULTI_HD constexpr auto num_elements() /*const*/ -> multi::ssize_t { return 1; }

	using indices_type = tuple<>;  // TODO(correaa) or boost::multi::detail::tuple<>; ?

	[[nodiscard]] static constexpr auto from_linear(nelems_type const& n) /*const*/ -> indices_type {
		assert(n == 0);
		(void)n;  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : constexpr function
		return indices_type{};
	}

	friend constexpr auto operator%(nelems_type const& n, extents_t const& /*s*/) -> tuple<> { return /*s.*/ from_linear(n); }

	static BOOST_MULTI_HD constexpr auto to_linear() /*const*/ -> difference_type { return 0; }
	BOOST_MULTI_HD constexpr auto        operator()() const { return to_linear(); }

	constexpr auto operator[](index) const -> element = delete;

	static BOOST_MULTI_HD constexpr auto next_canonical() /*const*/ -> bool { return true; }
	static BOOST_MULTI_HD constexpr auto prev_canonical() /*const*/ -> bool { return true; }

	friend constexpr auto intersection(extents_t const& /*x1*/, extents_t const& /*x2*/) -> extents_t { return {}; }

	constexpr BOOST_MULTI_HD auto operator==(extents_t const& /*other*/) const { return true; }
	constexpr BOOST_MULTI_HD auto operator!=(extents_t const& /*other*/) const { return false; }

	template<std::size_t Index>  // TODO(correaa) = detele ?
	friend constexpr auto get(extents_t const& self) -> typename std::tuple_element_t<Index, base_> {
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}

	template<std::size_t Index>  // TODO(correaa) = detele ?
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto get() const -> typename std::tuple_element_t<Index, base_> {
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}
};

template<> struct extents_t<1> : tuple<multi::index_extension> {
	using base_ = tuple<multi::index_extension>;

	static constexpr auto dimensionality = 1;  // TODO(correaa): consider deprecation

	constexpr static dimensionality_type rank_v = 1;

	using size_type = multi::index_extension::size_type;
	using difference_type = multi::index_extension::difference_type;
	using element = tuple<multi::index_extension::value_type>;
	using index = multi::index;
	using sizes_type = tuple<size_type>;

	[[deprecated("use .extent()")]] constexpr auto extension() const { using std::get; return get<0>(static_cast<base_ const&>(*this)); }
	[[nodiscard]] constexpr auto extent() const { using std::get; return get<0>(static_cast<base_ const&>(*this)); }

	constexpr auto sizes() const { return sizes_type{this->size()}; }  // using std::get; return get<0>(static_cast<base_ const&>(*this)); }

	constexpr auto sub() const { return extents_t<0>{this->base().tail()}; }

	class cursor_t {
		index idx_;
		extents_t<0> rest_;
		friend extents_t;

		constexpr cursor_t(index idx, extents_t<0> rest) : idx_{idx}, rest_{rest} {}

	 public:
		cursor_t() = default;
		using value_type = decltype(ht_tuple(std::declval<index>(), std::declval<extents_t<0>>().base()));
		using reference = value_type;

		BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> reference { return ht_tuple(idx_ + n, rest_.base()); }
	};

	auto home() const -> cursor_t {
		return cursor_t{this->base().head().first(), extents_t<0>{this->base().tail()}};
	}

	class iterator {  // : public weakly_incrementable<iterator> {
		index idx_;
		extents_t<0> rest_;
		friend extents_t;
	
		constexpr iterator(index idx, extents_t<0> rest) : idx_{idx}, rest_{rest} {}

	 public:
		iterator() = default;

		using difference_type = index;
		using value_type = decltype(ht_tuple(std::declval<index>(), std::declval<extents_t<0>>().base()));
		using pointer = void;
		using reference = value_type;
		using iterator_category = std::random_access_iterator_tag;

		constexpr auto operator+(difference_type n) const { return iterator{idx_ + n, rest_}; }
		constexpr auto operator-(difference_type n) const { return iterator{idx_ - n, rest_}; }

		friend BOOST_MULTI_HD constexpr auto operator-(iterator const& self, iterator const& other) -> difference_type { return self.idx_ - other.idx_; }
		friend BOOST_MULTI_HD constexpr auto operator+(difference_type n, iterator const& self) { return self + n; }

		constexpr auto operator+=(difference_type n) -> iterator& { idx_ += n; return *this; }
		constexpr auto operator-=(difference_type n) -> iterator& { idx_ -= n; return *this; }

		constexpr auto operator++() -> iterator& { ++idx_; return *this; }
		constexpr auto operator--() -> iterator& { --idx_; return *this; }

		constexpr auto operator++(int) -> iterator { iterator ret{*this}; operator++(); return ret; }  // NOLINT(cert-dcl21-cpp)
		constexpr auto operator--(int) -> iterator { iterator ret{*this}; operator--(); return ret; }  // NOLINT(cert-dcl21-cpp)

		constexpr auto operator*() const {
			// multi::detail::what(rest_);
			return ht_tuple(idx_, rest_.base());
			}

		BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> reference { return *(*this + n); }  // NOLINT(readability-redundant-parentheses) bug in clang-tidy trunk

		friend constexpr auto operator==(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ == other.idx_; }
		friend constexpr auto operator!=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ != other.idx_; }

		friend constexpr auto operator<(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ < other.idx_; }
		friend constexpr auto operator>(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ > other.idx_; }

		friend constexpr auto operator<=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ <= other.idx_; }
		friend constexpr auto operator>=(iterator const& self, iterator const& other) { assert( self.rest_ == other.rest_ ); return self.idx_ >= other.idx_; }
	};

	constexpr auto begin() const { return iterator{this->base().head().first(), extents_t<0>{this->base().tail()}}; }
	constexpr auto end()   const { return iterator{this->base().head().last() , extents_t<0>{this->base().tail()}}; }

	class elements_t {
		multi::index_range rng_;

	 public:
		class iterator : multi::index_range::iterator {
			friend class elements_t;  // enclosing class is friend automatically?
			BOOST_MULTI_HD constexpr explicit iterator(multi::index_range::iterator other)
			: multi::index_range::iterator{other} {}

			BOOST_MULTI_HD constexpr auto base_() const -> multi::index_range::iterator const& { return *this; }
			BOOST_MULTI_HD constexpr auto base_() -> multi::index_range::iterator& { return *this; }

		 public:
			using value_type      = std::tuple<multi::index_range::iterator::value_type>;
			using multi::index_range::iterator::difference_type;  // using difference_type = multi::index_range::iterator::difference_type;
			using reference = value_type;
			using pointer = void;

			iterator() = default;

			BOOST_MULTI_HD constexpr auto operator*() const -> reference { return *base_(); }

			BOOST_MULTI_HD constexpr auto operator++() -> iterator& {
				++base_();
				return *this;
			}

			BOOST_MULTI_HD constexpr auto operator--() -> iterator& {
				--base_();
				return *this;
			}

			BOOST_MULTI_HD constexpr auto operator++(int) { iterator ret{*this}; ++(*this); return ret; }
			BOOST_MULTI_HD constexpr auto operator--(int) { iterator ret{*this}; --(*this); return ret; }

			BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> iterator& {
				base_() += n;
				return *this;
			}

			BOOST_MULTI_HD constexpr auto operator-=(difference_type n) -> iterator& {
				base_() -= n;
				return *this;
			}

			BOOST_MULTI_HD constexpr auto operator+(difference_type n) const -> iterator { iterator ret{*this}; return ret += n; }  // mull-ignore: cxx_init_const
			BOOST_MULTI_HD constexpr auto operator-(difference_type n) const -> iterator { iterator ret{*this}; return ret -= n; }  // mull-ignore: cxx_init_const

			friend BOOST_MULTI_HD constexpr auto operator-(iterator const& self, iterator const& other) -> difference_type {
				return self.base_() - other.base_();
			}

			BOOST_MULTI_HD constexpr auto operator==(iterator const& other) const { return base_() == other.base_(); }
			BOOST_MULTI_HD constexpr auto operator!=(iterator const& other) const { return base_() != other.base_(); }

			BOOST_MULTI_HD constexpr auto operator<(iterator const& other) const { return base_() < other.base_(); }
			BOOST_MULTI_HD constexpr auto operator<=(iterator const& other) const { return base_() <= other.base_(); }

			BOOST_MULTI_HD auto operator[](difference_type n) const { return *(*this + n); }  // NOLINT(readability-redundant-parentheses) bug in clang-tidy trunk
		};

		BOOST_MULTI_HD constexpr auto begin() const noexcept -> iterator { return iterator{rng_.begin()}; }
		BOOST_MULTI_HD constexpr auto end() const noexcept -> iterator { return iterator{rng_.end()}; }

		using size_type = multi::index_extension::size_type;
		using difference_type = multi::index_extension::difference_type;
		using value_type      = iterator::value_type;
		using reference       = iterator::reference;

		BOOST_MULTI_HD constexpr auto operator[](difference_type n) const noexcept(noexcept(*(std::declval<iterator>() + n))) -> reference { return *(begin() + n); }  // NOLINT(readability-redundant-parentheses) bug in clang-tidy

		BOOST_MULTI_HD constexpr auto size() const -> size_type { return end() - begin(); }

		BOOST_MULTI_HD constexpr explicit elements_t(multi::index_range rng)
		: rng_{rng} {}
	};

	auto elements() const {
		using std::get;
		// auto rng = get<0>(static_cast<tuple<multi::index_extension> const&>(*this));
		return elements_t{get<0>(static_cast<tuple<multi::index_extension> const&>(*this))};
	}

	using nelems_type = index;

	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax (compatible with std::vector(int) constructor
	BOOST_MULTI_HD constexpr extents_t(multi::ssize_t size)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_(multi::index_extension{0, size}) {}

	template<class T1>
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int>  // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(tuple<T1> extensions)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_{static_cast<multi::index_extension>(extensions.head())} {}

	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(multi::index_extension const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_{other} {}

	BOOST_MULTI_HD constexpr explicit extents_t(base_ tup)
	: base_{tup} {}

	template<class OtherExtents,
		decltype( multi::detail::implicit_cast<multi::index_extension>(OtherExtents{}.extent()) )* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(OtherExtents const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_{other.extent()} {}

	extents_t() = default;

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return *this; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return *this; }

	BOOST_MULTI_HD constexpr auto operator==(extents_t const& other) const { return base() == other.base(); }
	BOOST_MULTI_HD constexpr auto operator!=(extents_t const& other) const { return base() != other.base(); }

	BOOST_MULTI_HD constexpr auto size() const noexcept -> size_type { return this->base().head().size(); }

	BOOST_MULTI_HD constexpr auto num_elements() const { return size(); }

	using indices_type = multi::detail::tuple<multi::index>;

	[[nodiscard]] BOOST_MULTI_HD constexpr auto from_linear(nelems_type const& n) const -> indices_type {  // NOLINT(readability-convert-member-functions-to-static) TODO(correaa)
		return indices_type{n};
	}

	friend constexpr auto operator%(nelems_type idx, extents_t const& extensions)
		-> multi::detail::tuple<multi::index> {
		return extensions.from_linear(idx);
	}

	static BOOST_MULTI_HD constexpr auto to_linear(index const& idx) -> difference_type { return idx; }

	BOOST_MULTI_HD constexpr auto operator[](index idx) const {
		using std::get;
		return multi::detail::tuple<multi::index>{get<0>(this->base())[idx]};
	}
	BOOST_MULTI_HD constexpr auto operator()(index idx) const { return idx; }

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto next_canonical(index& idx) const -> bool {  // NOLINT(google-runtime-references) idx is mutated
		using boost::multi::detail::get;
		// if(idx == ::boost::multi::detail::get<0>(this->base()).back()) {
		// 	idx = ::boost::multi::detail::get<0>(this->base()).first();
		// 	return true;
		// }
		++idx;
		if(idx == get<0>(this->base()).last()) {
			idx = get<0>(this->base()).first();
			return true;
		}
		return false;
	}
	constexpr auto prev_canonical(index& idx) const -> bool {  // NOLINT(google-runtime-references) idx is mutated
		using boost::multi::detail::get;
		if(idx == get<0>(this->base()).first()) {
			// idx = 42;  // TODO(correaa) implement and test
			idx = get<0>(this->base()).back();
			return true;
		}
		--idx;
		return false;
	}

	friend auto intersection(extents_t const& self, extents_t const& other) {
		return extents_t{
			intersection(
				boost::multi::detail::get<0>(self.base()),
				boost::multi::detail::get<0>(other.base())
			)
		};
	}
	template<class Archive>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		using boost::multi::detail::get;
		auto&  extension_ = get<0>(this->base());
		arxiv& multi::archive_traits<Archive>::make_nvp("extent", extension_);
	}

	template<std::size_t Index, std::enable_if_t<(Index < 1), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto get() const -> std::tuple_element_t<Index, base_> {  // by value, to match the other extents_t<D> specializations (structured bindings)
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < 1), int> = 0>                       // NOLINT(modernize-use-constraints) TODO(correaa)
	friend constexpr auto get(extents_t const& self) -> std::tuple_element_t<Index, base_> {  // by value, to match the other extents_t<D> specializations (structured bindings)
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}
};

template<dimensionality_type D> using iextensions = extents_t<D>;

template<dimensionality_type D> using extensions_t = extents_t<D>;  // backward-compatibility alias for the former name of extents_t

template<boost::multi::dimensionality_type D>
constexpr auto array_size_impl(boost::multi::extents_t<D> const&)
	-> std::integral_constant<std::size_t, static_cast<std::size_t>(D)>;

extents_t(multi::ssize_t) -> extents_t<1>;
extents_t(multi::ssize_t, multi::ssize_t) -> extents_t<2>;
extents_t(multi::ssize_t, multi::ssize_t, multi::ssize_t) -> extents_t<3>;
extents_t(multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t) -> extents_t<4>;
extents_t(multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t) -> extents_t<5>;
extents_t(multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t) -> extents_t<6>;
extents_t(multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t, multi::ssize_t) -> extents_t<7>;

}  // end namespace boost::multi

// Some versions of Clang throw warnings that stl uses class std::tuple_size instead
// of struct std::tuple_size like it should be
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

template<boost::multi::dimensionality_type D>
struct std::tuple_size<boost::multi::extents_t<D>>  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) to implement structured binding
: std::integral_constant<std::size_t, static_cast<std::size_t>(D)> {};

template<>
struct std::tuple_element<0, boost::multi::extents_t<0>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) to implement structured binding
	using type = void;
};

template<std::size_t Index, boost::multi::dimensionality_type D>
struct std::tuple_element<Index, boost::multi::extents_t<D>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) to implement structured binding
	using type = typename std::tuple_element_t<Index, typename boost::multi::extents_t<D>::base_>;
};

namespace std {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification)

// clang wants tuple_size to be a class, not a struct with -Wmismatched-tags
#if !defined(__GLIBCXX__) || (__GLIBCXX__ <= 20190406)
template<> struct tuple_size<boost::multi::extents_t<0>> : std::integral_constant<boost::multi::dimensionality_type, 0> {};
template<> struct tuple_size<boost::multi::extents_t<1>> : std::integral_constant<boost::multi::dimensionality_type, 1> {};
template<> struct tuple_size<boost::multi::extents_t<2>> : std::integral_constant<boost::multi::dimensionality_type, 2> {};
template<> struct tuple_size<boost::multi::extents_t<3>> : std::integral_constant<boost::multi::dimensionality_type, 3> {};
template<> struct tuple_size<boost::multi::extents_t<4>> : std::integral_constant<boost::multi::dimensionality_type, 4> {};
template<> struct tuple_size<boost::multi::extents_t<5>> : std::integral_constant<boost::multi::dimensionality_type, 5> {};
#else
template<> class tuple_size<boost::multi::extents_t<0>> : public std::integral_constant<boost::multi::dimensionality_type, 0> {};
template<> class tuple_size<boost::multi::extents_t<1>> : public std::integral_constant<boost::multi::dimensionality_type, 1> {};
template<> class tuple_size<boost::multi::extents_t<2>> : public std::integral_constant<boost::multi::dimensionality_type, 2> {};
template<> class tuple_size<boost::multi::extents_t<3>> : public std::integral_constant<boost::multi::dimensionality_type, 3> {};
template<> class tuple_size<boost::multi::extents_t<4>> : public std::integral_constant<boost::multi::dimensionality_type, 4> {};
template<> class tuple_size<boost::multi::extents_t<5>> : public std::integral_constant<boost::multi::dimensionality_type, 5> {};
#endif

#if !defined(_MSC_VER) && (!defined(__GLIBCXX__) || (__GLIBCXX__ <= 20240707))
template<std::size_t N, ::boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D> const& tp)  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal idiom to defined tuple get, gcc workaround
	-> decltype(tp.template get<N>()) {
	return tp.template get<N>();
}

// template<std::size_t N>  // , boost::multi::dimensionality_type D>
// constexpr auto get(boost::multi::extents_t<2> const& tp)  // NOLINT(cert-dcl58-cpp) normal idiom to defined tuple get, gcc workaround
// // ->decltype(tp.template get<N>()) {
// -> decltype(auto) {
//  return tp.template get<N>(); }

template<std::size_t N, ::boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D>& tp)  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal idiom to defined tuple get, gcc workaround
	-> decltype(tp.template get<N>()) {
	return tp.template get<N>();
}

template<std::size_t N, boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D>&& tp)  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal idiom to defined tuple get, gcc workaround
	-> decltype(std::move(tp).template get<N>()) {
	return std::move(tp).template get<N>();
}
#endif

template<typename Fn, boost::multi::dimensionality_type D>
constexpr auto
apply(Fn&& fun, boost::multi::extents_t<D> const& exts) noexcept -> decltype(auto) {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) workaround
	return exts.apply(std::forward<Fn>(fun));
}

}  // end namespace std

namespace boost::multi {

struct monostate : equality_comparable<monostate> {
	friend BOOST_MULTI_HD constexpr auto operator==(monostate const& /*self*/, monostate const& /*other*/) { return true; }
};

template<typename SSize = multi::index>
class stride_t {
	difference_type stride_;

 public:
	BOOST_MULTI_HD constexpr auto operator()() const -> difference_type { return stride_; }

	template<class Ptr>
	BOOST_MULTI_HD constexpr auto operator()(Ptr ptr) const -> Ptr { return ptr + stride_; }

	using category = std::random_access_iterator_tag;
};

template<typename SSize = multi::index>
class contiguous_stride_t {
 public:
	// using difference_type = SSize;

	BOOST_MULTI_HD constexpr auto operator()() const -> SSize { return 1; }

	template<class Ptr>
	BOOST_MULTI_HD constexpr auto operator()(Ptr const& ptr) const -> Ptr { return ptr + 1; }

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
	using category = std::random_access_iterator_tag;  // std::contiguous_iterator_tag;
#else
	using category = std::random_access_iterator_tag;
#endif
};

// using multi::detail::tuple;

template<typename SSize = multi::index>
class contiguous_layout {

 public:
	using dimensionality_type                           = multi::dimensionality_t;
	using rank                                          = std::integral_constant<dimensionality_type, 1>;
	static constexpr auto                rank_v         = rank::value;
	static constexpr dimensionality_type dimensionality = rank_v;

	using size_type  = SSize;
	using sizes_type = typename boost::multi::detail::tuple<size_type>;

	using difference_type = SSize;

	using index           = size_type;
	using index_range     = multi::range<index>;
	using index_extension = multi::extent_t<index>;

	using indexes = tuple<index>;

	using extent_type = multi::extent_t<index>;
	using extension_type [[deprecated]] = multi::extent_t<index>;
	using extents_type = multi::extents_t<1>;
	using extensions_type [[deprecated("use extent_type")]] = multi::extents_t<1>;

	using stride_type  = std::integral_constant<int, 1>;
	using strides_type = boost::multi::detail::tuple<stride_type>;

	using offset_type = std::integral_constant<int, 0>;

	using nelems_type = SSize;

	using sub_type = layout_t<0, SSize>;

 private:
	// BOOST_MULTI_NO_UNIQUE_ADDRESS sub_type sub_;
	// BOOST_MULTI_NO_UNIQUE_ADDRESS stride_type stride_;
	size_type nelems_;

	template<std::size_t N, class Tup>
	static constexpr auto get_(Tup&& tup) {
		using std::get;
		return get<N>(std::forward<Tup>(tup));
	}

 public:
	constexpr explicit contiguous_layout(multi::extents_t<1> exts)
	: nelems_{get_<0>(exts).size()} {}

	BOOST_MULTI_HD constexpr contiguous_layout(
		sub_type /*sub*/,
		stride_type /*stride*/,
		offset_type /*offset*/,
		nelems_type nelems
	)
	: /*sub_{sub}, stride_{} offset_{},*/ nelems_{nelems} {}

 private:
	static constexpr auto at_aux_(index /*idx*/) {
		return sub_type{};  // sub_.sub_, sub_.stride_, sub_.offset_ + offset_ + (idx*stride_), sub_.nelems_}();
	}

 public:
	BOOST_MULTI_HD constexpr auto operator[](index idx) const { return at_aux_(idx); }

	template<typename... Indices>
	BOOST_MULTI_HD constexpr auto operator()(index idx, Indices... rest) const { return operator[](idx)(rest...); }
	BOOST_MULTI_HD constexpr auto operator()(index idx) const { return at_aux_(idx); }
	BOOST_MULTI_HD constexpr auto operator()() const { return *this; }

	static BOOST_MULTI_HD constexpr auto stride() { return std::integral_constant<int, 1>{}; }
	static BOOST_MULTI_HD constexpr auto offset() { return std::integral_constant<int, 0>{}; }

	BOOST_MULTI_HD constexpr auto extent() const { return extent_type{0, nelems_}; }
	[[deprecated]] BOOST_MULTI_HD constexpr auto extension() const { return extent_type{0, nelems_}; }

	BOOST_MULTI_HD constexpr auto num_elements() const { return nelems_; }

	BOOST_MULTI_HD constexpr auto size() const { return nelems_; }
	BOOST_MULTI_HD constexpr auto sizes() const { return sizes_type{size()}; }

	BOOST_MULTI_HD constexpr auto nelems() const noexcept { return nelems_; }

	[[deprecated]] BOOST_MULTI_HD constexpr auto extensions() const { return multi::extents_t<1>{extent()}; }
	BOOST_MULTI_HD constexpr auto extents() const { return multi::extents_t<1>{extent()}; }

	BOOST_MULTI_HD constexpr auto is_empty() const -> bool { return nelems_ == 0; }

	BOOST_MULTI_NODISCARD("empty checks for emptyness, it performs no action. Use `is_empty()` instead")
	BOOST_MULTI_HD constexpr auto empty() const { return is_empty(); }

	static constexpr auto sub() { return layout_t<0, SSize>{}; }

	static constexpr auto is_compact() { return std::true_type{}; }

	BOOST_MULTI_HD constexpr auto drop(difference_type count) const {
		assert(count <= this->size());

		return contiguous_layout(
			/*this->*/sub(),
			/*this->*/stride(),
			/*this->*/offset(),
			/*this->*/stride() * (this->size() - count)
		);
	}

	BOOST_MULTI_HD constexpr auto slice(index first, index last) const {
		return contiguous_layout(
			/*this->*/ sub(),
			/*this->*/ stride(),
			/*this->*/ offset(),
			this->is_empty() ? 0 : this->nelems() / this->size() * (last - first)
		);
	}
};

template<typename Stride1, typename Stride2, typename Size1, typename Pointer = void*>
class bistride {
 public:
	using stride1_type = Stride1;
	using size1_type   = Size1;
	using stride2_type = Stride2;
	using offset_type  = std::ptrdiff_t;

 private:
	stride1_type stride1_;
	stride2_type stride2_;
	multi::ssize_t nelems2_;
	Pointer ptr_;
	std::ptrdiff_t n_;

 public:
	auto stride1() const { return stride1_; }
	auto stride2() const { return stride2_; }

	auto nelems2() const { return nelems2_; }  // cppcheck-suppress functionStatic;  // bug in cppcheck 2,.20

	using category = std::random_access_iterator_tag;

	BOOST_MULTI_HD constexpr explicit bistride(stride1_type stride1, stride2_type stride2, multi::ssize_t size, Pointer ptr)  // NOLINT(bugprone-easily-swappable-parameters)
	: stride1_{stride1}, stride2_{stride2}, nelems2_{size}, ptr_{ptr}, n_{1} {}

	BOOST_MULTI_HD constexpr explicit bistride(stride1_type stride1, stride2_type stride2, multi::ssize_t size, Pointer ptr, std::ptrdiff_t n)  // NOLINT(bugprone-easily-swappable-parameters)
	: stride1_{stride1}, stride2_{stride2}, nelems2_{size}, ptr_{ptr}, n_{n} {}

	BOOST_MULTI_HD constexpr auto operator*(std::ptrdiff_t n) const {
		assert(n_ == 1);  // TODO(correaa) test n_ != 1
		return bistride{stride1_, stride2_, nelems2_, ptr_, n /**n_*/};
	}

	#if (defined(__clang__) && (__clang_major__ >= 16)) && !defined(__INTEL_LLVM_COMPILER)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif
	template<class Ptr>
	friend BOOST_MULTI_HD constexpr auto operator+=(Ptr& ptr, bistride const& self) -> Ptr& {
		ptr = ptr + self;
		return ptr;
	}

	template<class Ptr>
	friend BOOST_MULTI_HD constexpr auto operator+=(Ptr& ptr, bistride& self) -> Ptr& {
		assert(self.n_ == 1);
		// if(self.n_ == 1) {
		// 	ptr += self.stride2_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		// 	if(ptr == static_cast<Ptr>(self.ptr_) + self.nelems2_) {  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		// 		self.ptr_ = static_cast<Ptr>(self.ptr_) + self.stride1_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		// 		ptr = static_cast<Ptr>(self.ptr_);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		// 	}
		// } else {
			ptr = ptr + self;
		// }
		return ptr;
	}

	BOOST_MULTI_HD constexpr auto operator-(offset_type /*unused*/) const { return *this; }

	template<class Ptr>
	friend BOOST_MULTI_HD constexpr auto operator+(Ptr const& ptr, bistride const& self) {
		auto base = static_cast<Ptr>(self.ptr_);
		auto dist = ptr - base;
		auto outer = dist / self.stride1_;

		// vvv TODO(correaa) Survived: Replaced / with *
		auto inner = (dist % self.stride1_) / self.stride2_;  // mull-ignore: cxx_div_to_mul

		auto shift = inner + self.n_;
		auto size2 = self.nelems2_ / self.stride2_;

		auto new_outer = (shift / size2) + outer;
		auto new_inner = shift % size2;

		return base + (new_outer * self.stride1_) + (new_inner * self.stride2_);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}

	template<class Ptr>
	BOOST_MULTI_HD constexpr auto segment_base(Ptr const& ptr) const {
		auto base = static_cast<Ptr>(ptr_);
		auto dist = ptr - base;
		auto segment_index = dist / stride1_;
		auto ret = base + (segment_index * stride1_);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return ret;
	}
	#if (defined(__clang__) && (__clang_major__ >= 16)) && !defined(__INTEL_LLVM_COMPILER)
	#pragma clang diagnostic pop
	#endif
};

template<dimensionality_type D>
struct bilayout {
	using size_type       = multi::ssize_t;
	using difference_type = std::make_signed_t<size_type>;
	using index           = difference_type;

	using stride1_type = difference_type;
	using stride2_type = difference_type;
	// using bistride_type = std::pair<index, index>;
	using sub_type = layout_t<D - 1>;

	using dimensionality_type    = typename sub_type::dimensionality_type;
	using rank                   = std::integral_constant<dimensionality_type, sub_type::rank::value + 1>;
	constexpr static auto rank_v = rank::value;

	constexpr static auto dimensionality() { return rank_v; }

 private:
	stride1_type stride1_;
	size_type    nelems1_;
	stride2_type stride2_;
	size_type    nelems2_;
	sub_type     sub_;
	void* ptr_;

 public:
	bilayout(
		stride1_type stride1,  // NOLINT(bugprone-easily-swappable-parameters)
		size_type    nelems1,
		stride2_type stride2,  // NOLINT(bugprone-easily-swappable-parameters)
		size_type    nelems2,
		sub_type     sub,
		void* ptr
	)
	: stride1_{stride1}, nelems1_{nelems1}, stride2_{stride2}, nelems2_{nelems2}, sub_{std::move(sub)}, ptr_{ptr} {}

	using offset_type     = std::ptrdiff_t;
	// using stride_type     = void;  // std::pair<stride1_type, stride2_type>;

	using stride_type = bistride<stride1_type, stride2_type, size_type>;

	using index_range     = multi::range<index>;
	using extent_type  = void;
	using extension_type [[deprecated]] = void;

	using extents_type = extents_t<D>;
	using extensions_type [[deprecated("use extents_t")]] = extents_type;

	// using extents_type  = void;
	// using extensions_type [[deprecated]] = void;

	using sizes_type      = void;
	using indexes         = void;

	using strides_type = void;

	// auto stride() const = delete;
	BOOST_MULTI_HD constexpr auto stride() const {
		return stride_type{stride1_, stride2_, nelems2_, ptr_, 1};
	}
	auto num_elements() const = delete;

	BOOST_MULTI_HD static constexpr auto offset() { return offset_type{}; }
	BOOST_MULTI_HD constexpr auto size() const { return (nelems2_ / stride2_) * (nelems1_ / stride1_); }

	constexpr auto nelems() const { return nelems1_ - stride1_ + nelems2_; }  // span to one-past-end: (nsegs-1) outer strides + the last segment (NOT just nelems2_, which is one segment short)
	[[deprecated("use extent")]] void extension() const  = delete;
	constexpr auto extent() const { return multi::extent_t{0, size()}; }

	[[deprecated("use extents")]] auto extensions() const = delete;
	auto extents() const = delete;

	auto is_empty() const   = delete;
	auto empty() const      = delete;
	BOOST_MULTI_HD constexpr auto sub() const { return sub_; }
	auto sizes() const      = delete;

	auto is_compact() const = delete;

	using index_extension = multi::index_extension;
};

template<class Ptr>
class segmented_ptr {
	Ptr ptr_;
	Ptr first_;
	Ptr last_;
	std::ptrdiff_t stride_;

 public:
	segmented_ptr(Ptr ptr, Ptr first, Ptr last, std::ptrdiff_t stride) : ptr_{ptr}, first_{first}, last_{last}, stride_{stride} {}
	auto operator++() -> segmented_ptr& {
		++ptr_;
		if(ptr_ == last_) {
			first_ += stride_;
			last_ += stride_;
			ptr_ = first_;
		}
		return *this;
	}

	auto operator--() -> segmented_ptr& {
		if(ptr_ == first_) {
			first_ -= stride_;
			last_ -= stride_;
			ptr_ = last_;
		}
		--ptr_;
		return *this;
	}
};

template<dimensionality_type D, typename SSize>
struct layout_t
	: multi::equality_comparable<layout_t<D, SSize>> {
	template<class Ptr = void*>
	auto flatten(Ptr ptr) const {
		return bilayout<D - 1>(
			stride(),
			nelems(),
			sub().stride(),
			sub().nelems(),
			sub().sub(),
			ptr
		);
	}

	using dimensionality_type = multi::dimensionality_type;
	using rank                = std::integral_constant<dimensionality_type, D>;

	using sub_type        = layout_t<D - 1>;
	using size_type       = SSize;
	using difference_type = std::make_signed_t<size_type>;

	/// indexing type in the leading dimension (usually `std::ptrdiff_t`)
	using index           = difference_type;

	using index_extension = multi::index_extension;
	using index_range     = multi::range<index>;

	using stride_type = index;
	using offset_type = index;
	using nelems_type = index;

	using strides_type = typename boost::multi::detail::tuple_prepend<stride_type, typename sub_type::strides_type>::type;
	using offsets_type = typename boost::multi::detail::tuple_prepend<offset_type, typename sub_type::offsets_type>::type;
	using nelemss_type = typename boost::multi::detail::tuple_prepend<nelems_type, typename sub_type::nelemss_type>::type;

	using extent_type = index_extension;  // not index_range!
	using extension_type [[deprecated]] = index_extension;  // not index_range!

	using extensions_type [[deprecated]] = extents_t<rank::value>;
	using extents_type = extents_t<rank::value>;

	using sizes_type      = typename boost::multi::detail::tuple_prepend<size_type, typename sub_type::sizes_type>::type;

	using indexes = typename boost::multi::detail::tuple_prepend<index, typename sub_type::indexes>::type;

	static constexpr dimensionality_type rank_v         = rank::value;
	static constexpr dimensionality_type dimensionality = rank_v;  // TODO(correaa): consider deprecation

	[[deprecated("for compatibility with Boost.MultiArray, use static `dimensionality` instead")]]
	static constexpr auto num_dimensions() { return dimensionality; }  // NOSONAR(cpp:S1133)

	friend constexpr auto dimensionality(layout_t const& /*self*/) { return rank_v; }

 private:
	sub_type    sub_;
	stride_type stride_;  // =  1;  // or std::numeric_limits<stride_type>::max()?
	offset_type offset_;
	nelems_type nelems_;

	template<dimensionality_type, typename> friend struct layout_t;

 public:
	layout_t() = default;

	template<
		class OtherLayout,
		class = decltype(sub_type{std::declval<OtherLayout const&>().sub()}),
		class = decltype(stride_type{std::declval<OtherLayout const&>().stride()}),
		class = decltype(offset_type{std::declval<OtherLayout const&>().offset()}),
		class = decltype(nelems_type{std::declval<OtherLayout const&>().nelems()})>
	BOOST_MULTI_HD constexpr explicit layout_t(OtherLayout const& other)
	: sub_{other.sub()}, stride_{other.stride()}, offset_{other.offset()}, nelems_{other.nelems()} {}

 private:
	template<class Fun, class Tup>
	static BOOST_MULTI_HD constexpr auto apply_(Fun&& fun, Tup&& tup) -> decltype(auto) {  // this is workaround for icc 2021
		using std::apply;
		return apply(std::forward<Fun>(fun), std::forward<Tup>(tup));
	}

 public:
	#ifdef __NVCC__
	#pragma nv_diagnostic push
	#pragma nv_diag_suppress = 20013  // TODO(correa) use multi::apply  // calling a constexpr __host__ function("apply") from a __host__ __device__ function("layout_t") is not allowed.
	#endif
 private:
	template<class... Args>
	static BOOST_MULTI_HD constexpr auto std_apply_(Args&&... args) ->decltype(auto) { using std::apply; return apply(std::forward<Args>(args)...); }

 public:
	BOOST_MULTI_HD constexpr explicit layout_t(extents_type const& exts)
	: sub_{apply_        ([](auto const&... subexts) -> auto { return multi::extents_t<D - 1>{subexts...}; }, detail::tail(exts.base()))}
	// : sub_{/*std::*/apply([](auto const&... subexts) { return multi::extents_t<D - 1>{subexts...}; }, detail::tail(exts.base()))}
	, stride_{sub_.num_elements() ? sub_.num_elements() : 1}
	, offset_{boost::multi::detail::get<0>(exts.base()).first() * stride_}
	, nelems_{boost::multi::detail::get<0>(exts.base()).size() * sub().num_elements()} {}

	BOOST_MULTI_HD constexpr explicit layout_t(extents_type const& exts, strides_type const& strds)
	: sub_{std::apply([](auto const&... subexts) -> auto { return multi::extents_t<D - 1>{subexts...}; }, detail::tail(exts.base())), detail::tail(strds)}, stride_{boost::multi::detail::get<0>(strds)}, offset_{boost::multi::detail::get<0>(exts.base()).first() * stride_}, nelems_{boost::multi::detail::get<0>(exts.base()).size() * sub().num_elements()} {}
	#ifdef __NVCC__
	#pragma nv_diagnostic pop
	#endif

	BOOST_MULTI_HD constexpr explicit layout_t(sub_type const& sub, stride_type str, offset_type off, nelems_type nlms)  // NOLINT(bugprone-easily-swappable-parameters)
	: sub_{sub}, stride_{str}, offset_{off}, nelems_{nlms} {}

	BOOST_MULTI_HD constexpr explicit layout_t(sub_type const& sub, stride_type str, offset_type off /*, nelems_type nelems*/)  // NOLINT(bugprone-easily-swappable-parameters)
	: sub_{sub}, stride_{str}, offset_{off} /*, nelems_{nelems}*/ {}                                                            // this leaves nelems_ uninitialized

	constexpr auto origin() const { return sub_.origin() - offset_; }

 private:
	#ifdef __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wlarge-by-value-copy"
	#endif

	BOOST_MULTI_HD constexpr auto at_aux_(index idx) const {
		return sub_type{sub_.sub_, sub_.stride_, sub_.offset_ + offset_ + (idx * stride_), sub_.nelems_}();
	}

 public:
	BOOST_MULTI_HD constexpr auto operator[](index idx) const { return at_aux_(idx); }

	template<typename... Indices>
	BOOST_MULTI_HD constexpr auto operator()(index idx, Indices... rest) const { return operator[](idx)(rest...); }
	BOOST_MULTI_HD constexpr auto operator()(index idx) const { return at_aux_(idx); }

	#ifdef __clang__
	#pragma clang diagnostic pop
	#endif

	#ifdef __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wlarge-by-value-copy"  // TODO(correaa) can it be returned by reference?
	#endif

	BOOST_MULTI_HD constexpr auto operator()() const { return *this; }

	#ifdef __clang__
	#pragma clang diagnostic pop
	#endif

	BOOST_MULTI_HD constexpr auto        sub() & -> sub_type& { return sub_; }
	BOOST_MULTI_HD constexpr auto        sub() const& -> sub_type const& { return sub_; }
	friend BOOST_MULTI_HD constexpr auto sub(layout_t const& self) -> sub_type const& { return self.sub(); }

	BOOST_MULTI_HD constexpr auto        nelems() & -> nelems_type& { return nelems_; }
	BOOST_MULTI_HD constexpr auto        nelems() const& -> nelems_type const& { return nelems_; }
	// friend BOOST_MULTI_HD constexpr auto nelems(layout_t const& self) -> nelems_type const& { return self.nelems(); }

	// constexpr BOOST_MULTI_HD auto nelems(dimensionality_type dim) const { return (dim != 0) ? sub_.nelems(dim - 1) : nelems_; }

	friend BOOST_MULTI_HD constexpr auto operator==(layout_t const& self, layout_t const& other) -> bool {
		return self.sub_ == other.sub_ && self.stride_ == other.stride_ && self.offset_ == other.offset_ && self.nelems_ == other.nelems_;
		// return std::tie(self.sub_, self.stride_, self.offset_, self.nelems_) == std::tie(other.sub_, other.stride_, other.offset_, other.nelems_);
	}

	friend BOOST_MULTI_HD constexpr auto operator!=(layout_t const& self, layout_t const& other) -> bool {
		return !(self == other);
		// return std::tie(self.sub_, self.stride_, self.offset_, self.nelems_) != std::tie(other.sub_, other.stride_, other.offset_, other.nelems_);
	}

	constexpr BOOST_MULTI_HD auto operator<(layout_t const& other) const -> bool {
		return std::tie(sub_, stride_, offset_, nelems_) < std::tie(other.sub_, other.stride_, other.offset_, other.nelems_);
	}

	constexpr auto reindex() const { return *this; }
	constexpr auto reindex(index idx) const {
		return layout_t(
			sub(),
			stride(),
			idx * stride(),
			nelems()
		);
	}
	template<class... Indexes>
	constexpr auto reindexed(index first, Indexes... idxs) const {
		return reindexed(first).rotate().reindexed(idxs...).unrotate();
	}

	BOOST_MULTI_HD constexpr auto num_elements() const noexcept -> size_type { return size() * sub_.num_elements(); }  // TODO(correaa) investigate mutation * -> /

	BOOST_MULTI_HD constexpr auto is_empty() const noexcept { return nelems_ == 0; }  // mull-ignore: cxx_eq_to_ne

	BOOST_MULTI_NODISCARD("empty checks for emptyness, it performs no action. Use `is_empty()` for clarity instead")
	BOOST_MULTI_HD constexpr auto empty() const noexcept { return is_empty(); }

	BOOST_MULTI_HD constexpr  auto size() const noexcept -> size_type {
		if(nelems_ == 0) {
			return 0;
		}
		// BOOST_MULTI_ACCESS_ASSERT(stride_);
		// if(nelems_ != 0) {MULTI_ACCESS_ASSERT(stride_ != 0);}
		// return nelems_ == 0?0:nelems_/stride_;
		// assert(stride_ != 0);
		return nelems_ / stride_;
	}

	BOOST_MULTI_HD constexpr auto stride() -> stride_type& { return stride_; }
	BOOST_MULTI_HD constexpr auto stride() const -> stride_type const& { return stride_; }

	BOOST_MULTI_HD constexpr auto        strides() const -> strides_type { return strides_type{stride(), sub_.strides()}; }

	BOOST_MULTI_HD constexpr auto        offset() const -> index { return offset_; }

	constexpr BOOST_MULTI_HD auto        offsets() const { return boost::multi::detail::tuple{offset(), sub_.offsets()}; }
	constexpr BOOST_MULTI_HD auto        nelemss() const { return boost::multi::detail::tuple{nelems(), sub_.nelemss()}; }

 private:
	constexpr auto base_size_() const {
		return (std::max)(nelems_, sub_.base_size_());
	}

 public:
	constexpr auto is_compact() const { return base_size_() == num_elements(); }

	constexpr auto is_flattable() const & {
		return (this->size() <= 1) || (this->stride() == this->sub().nelems());
	}

	constexpr auto shape() const& -> decltype(auto) { return sizes(); }

	BOOST_MULTI_HD constexpr auto sizes() const noexcept { return multi::detail::ht_tuple(size(), sub_.sizes()); }

	[[nodiscard]] BOOST_MULTI_HD constexpr auto extent() const -> extent_type {
		if(nelems_ == 0) {
			return index_extension{};
		}
		// assert(stride_ != 0);
		assert(offset_ % stride_ == 0);
		assert(nelems_ % stride_ == 0);
		return index_extension{offset_ / stride_, (offset_ + nelems_) / stride_};
	}
	[[deprecated]] [[nodiscard]] BOOST_MULTI_HD constexpr auto extension() const -> extent_type {
		return extent();
	}

	BOOST_MULTI_HD constexpr auto        extents() const {
		// auto fa = extension();
		// auto sa = sub_.extensions().base();
		// auto ht_tuple = multi::detail::ht_tuple(fa, sa);
		// auto ret = extensions_type{ht_tuple};
		// return ret;
		return extents_type{multi::detail::ht_tuple(extent(), sub_.extents().base())};
	}

	[[deprecated]] BOOST_MULTI_HD constexpr auto        extensions() const {
		// auto fa = extension();
		// auto sa = sub_.extensions().base();
		// auto ht_tuple = multi::detail::ht_tuple(fa, sa);
		// auto ret = extensions_type{ht_tuple};
		// return ret;
		return extents_type{multi::detail::ht_tuple(extent(), sub_.extents().base())};
	}

	[[deprecated("use get<d>(m.extensions()")]]  // TODO(correaa) redeprecate, this is commented to give a smaller CI output
	constexpr auto
	extension(dimensionality_type dim) const {
		return std::apply([](auto... extensions) -> auto { return std::array<index_extension, static_cast<std::size_t>(D)>{{extensions...}}; }, extensions().base()).at(static_cast<std::size_t>(dim));
	}  // cppcheck-suppress syntaxError ; bug in cppcheck 2.14
	   //  [[deprecated("use get<d>(m.strides())  ")]]  // TODO(correaa) redeprecate, this is commented to give a smaller CI output
	// constexpr auto stride(dimensionality_type dim) const {
	// 	return std::apply([](auto... strides) -> auto { return std::array<stride_type, static_cast<std::size_t>(D)>{{strides...}}; }, strides()).at(static_cast<std::size_t>(dim));
	// }
	//  [[deprecated("use get<d>(m.sizes())    ")]]  // TODO(correaa) redeprecate, this is commented to give a smaller CI output
	//  constexpr auto size     (dimensionality_type dim) const {return std::apply([](auto... sizes     ) {return std::array<size_type      , static_cast<std::size_t>(D)>{sizes     ...};}, sizes     ()       ).at(static_cast<std::size_t>(dim));}

	BOOST_MULTI_HD constexpr auto sort() const {
		auto ret = layout_t(
			this->sub().sort(),
			this->stride(),
			this->offset(),
			this->nelems()
		);

		if constexpr(D > 1) {
			if(
				    ret.stride() < ret.sub().stride()  // if strides are equal, one of the sizes is 1  // mull-ignore: cxx_lt_to_le
				|| (
					   ret.stride() == ret.sub().stride()
					&& ret.sub().size() < ret.size()  // mull-ignore: cxx_lt_to_le
				)
			)
			{
				auto ret2 = ret.transpose();
				ret = layout_t(
					ret2.sub().sort(),
					ret2.stride(),
					ret2.offset(),
					ret2.nelems()
				);
			}
		}

		return ret;
	}

	BOOST_MULTI_HD constexpr auto drop(difference_type count) const {
		assert(count <= this->size());

		return layout_t{
			this->sub(),
			this->stride(),
			this->offset(),
			this->stride() * (this->size() - count),
		};
	}

	BOOST_MULTI_HD constexpr auto slice(index first, index last) const {
		return layout_t(
			this->sub(),
			this->stride(),
			this->offset(),
			this->is_empty() ? 0 : this->nelems() / this->size() * (last - first)
		);
	}

	constexpr auto partition(size_type n) const {
		assert(n != 0);
		// vvv TODO(correaa) should be size() here?
		assert((this->nelems() % n) == 0);  // if you get an assertion here it means that you are partitioning an array with an incommunsurate partition
		return multi::layout_t<D + 1>{
			multi::layout_t<D>{
				this->sub(),
				this->stride(),
				this->offset(),
				this->nelems() / n,  // mull-ignore: cxx_div_to_mul
			},
			this->nelems() / n,  // mull-ignore: cxx_div_to_mul
			0,
			this->nelems(),
		};
		// new_layout.sub().nelems() /= n;
	}

	// template<class TT>
	// constexpr static void ce_swap(TT& left, TT& right) {
	// 	TT tmp = std::move(left);
	// 	left   = std::move(right);
	// 	right  = tmp;
	// }

	BOOST_MULTI_HD constexpr auto transpose() const {
		return layout_t(
			sub_type(
				sub().sub(),
				stride(),
				offset(),
				nelems()
			),
			sub().stride(),
			sub().offset(),
			sub().nelems()
		);
	}

	constexpr auto reverse() const {
		auto ret = unrotate();
		return layout_t(
			ret.sub().reverse(),
			ret.stride(),
			ret.offset(),
			ret.nelems()
		);
	}

	BOOST_MULTI_HD constexpr auto rotate() const {
		if constexpr(D > 1) {
			auto const ret = transpose();
			return layout_t(
				ret.sub().rotate(),
				ret.stride(),
				ret.offset(),
				ret.nelems()
			);
		} else {
			return *this;
		}
	}

	BOOST_MULTI_HD constexpr auto unrotate() const {
		if constexpr(D > 1) {
			auto const ret = layout_t(
				sub().unrotate(),
				stride(),
				offset(),
				nelems()
			);
			return ret.transpose();
		} else {
			return *this;
		}
	}

	constexpr auto hull_size() const -> size_type {
		if(is_empty()) {
			return 0;
		}
		return std::abs(size() * stride()) > std::abs(sub_.hull_size()) ? size() * stride() : sub_.hull_size();
	}

	[[deprecated("use two arg version")]] constexpr auto scale(size_type factor) const {
		return layout_t{sub_.scale(factor), stride_ * factor, offset_ * factor, nelems_ * factor};
	}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wlarge-by-value-copy"  // TODO(correaa) use checked span
#endif

	BOOST_MULTI_HD constexpr auto take(size_type n) const {
		return layout_t(
			this->sub(),
			this->stride(),
			this->offset(),
			this->stride() * n
		);
	}

	BOOST_MULTI_HD constexpr auto halve() const {
		assert(this->size() % 2 == 0);
		return layout_t<D + 1>(
			this->take(this->size() / 2),
			this->nelems() / 2,
			0,
			this->nelems()
		);
	}

	constexpr auto scale(size_type num, size_type den) const {
		assert((stride_ * num) % den == 0);
		assert(offset_ == 0);  // TODO(correaa) implement ----------------vvv
		return layout_t{sub_.scale(num, den), stride_ * num / den, offset_ /* *num/den */, nelems_ * num / den};
	}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
};

template<typename SSize>
struct layout_t<0, SSize>
: multi::equality_comparable<layout_t<0, SSize>> {
	using dimensionality_type = multi::dimensionality_type;
	using rank                = std::integral_constant<dimensionality_type, 0>;

	using size_type       = SSize;
	using difference_type = std::make_signed_t<size_type>;
	using index           = difference_type;
	using index_extension = multi::index_extension;
	using index_range     = multi::range<index>;

	using sub_type    = monostate;
	using stride_type = monostate;
	using offset_type = index;
	using nelems_type = index;

	using strides_type = tuple<>;
	using offsets_type = tuple<>;
	using nelemss_type = tuple<>;

	using extent_type = void;
	using extension_type [[deprecated]] = void;

	using extents_type    = extents_t<rank::value>;
	using extensions_type [[deprecated]] = extents_t<rank::value>;
	using sizes_type      = tuple<>;
	using indexes         = tuple<>;

	static constexpr dimensionality_type rank_v         = rank::value;
	static constexpr dimensionality_type dimensionality = rank_v;  // TODO(correaa) : consider deprecation

	friend constexpr auto dimensionality(layout_t const& /*self*/) { return rank_v; }

 private:
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // '6' bytes padding added after data member
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif
	BOOST_MULTI_NO_UNIQUE_ADDRESS sub_type    sub_;
	BOOST_MULTI_NO_UNIQUE_ADDRESS stride_type stride_;  // TODO(correaa) padding struct 'boost::multi::layout_t<0>' with 1 byte to align 'stride_' [-Werror,-Wpadded]

	offset_type offset_;

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

	nelems_type nelems_;

	template<dimensionality_type, typename> friend struct layout_t;

 public:
	layout_t() = default;

	BOOST_MULTI_HD constexpr explicit layout_t(extents_type const& /*nil*/)
	: offset_{0}, nelems_{1} {}

	// BOOST_MULTI_HD constexpr explicit layout_t(extensions_type const& /*nil*/, strides_type const& /*nil*/) {}

	BOOST_MULTI_HD constexpr layout_t(sub_type sub, stride_type stride, offset_type offset, nelems_type nelems)  // NOLINT(bugprone-easily-swappable-parameters)
	: sub_{sub}, stride_{stride}, offset_{offset}, nelems_{nelems} {}

	[[nodiscard]] BOOST_MULTI_HD constexpr auto extensions() const { return extents_type{}; }  // cppcheck-suppress functionStatic
	// friend BOOST_MULTI_HD constexpr auto        extensions(layout_t const& self) { return self.extensions(); }
	[[nodiscard]] BOOST_MULTI_HD constexpr auto extents() const { return extents_type{}; }  // cppcheck-suppress functionStatic

	[[nodiscard]] BOOST_MULTI_HD constexpr auto num_elements() const { return nelems_; }
	// friend BOOST_MULTI_HD constexpr auto        num_elements(layout_t const& self) { return self.num_elements(); }

	[[nodiscard]] BOOST_MULTI_HD constexpr auto sizes() const { return tuple<>{}; }
	// friend BOOST_MULTI_HD constexpr auto        sizes(layout_t const& self) { return self.sizes(); }

	[[nodiscard]] BOOST_MULTI_HD constexpr auto strides() const { return strides_type{}; }  // cppcheck-suppress functionStatic
	[[nodiscard]] BOOST_MULTI_HD constexpr auto offsets() const { return offsets_type{}; }  // cppcheck-suppress functionStatic
	[[nodiscard]] BOOST_MULTI_HD constexpr auto nelemss() const { return nelemss_type{}; }  // cppcheck-suppress functionStatic

	BOOST_MULTI_HD constexpr auto operator()() const { return offset_; }
	// constexpr explicit operator offset_type() const {return offset_;}

	constexpr auto stride() const -> stride_type = delete;
	constexpr auto offset() const -> offset_type { return offset_; }
	constexpr auto nelems() const noexcept -> nelems_type { return nelems_; }
	constexpr auto sub() const -> sub_type = delete;

	constexpr auto sort() const noexcept { return *this; }

	constexpr auto size() const -> size_type           = delete;

	[[deprecated]] constexpr auto extension() const -> extent_type = delete;
	constexpr auto extent() const -> extent_type = delete;

	BOOST_MULTI_HD constexpr auto is_empty() const noexcept { return nelems_ == 0; }

	BOOST_MULTI_NODISCARD("empty checks for emptyness, it performs no action. Use `is_empty()` instead")
	constexpr auto empty() const noexcept { return nelems_ == 0; }

	friend constexpr auto empty(layout_t const& self) noexcept { return self.empty(); }

	[[deprecated("is going to be removed")]]
	constexpr auto is_compact() const -> bool = delete;

 private:
	static constexpr auto base_size_() -> size_type { return 0; }

 public:
	static constexpr auto origin() -> offset_type { return 0; }

	constexpr auto reverse() const { return *this; }
	// constexpr auto reverse()          -> layout_t& {return *this;}

	BOOST_MULTI_HD constexpr auto take(size_type /*n*/) const {
		return layout_t<0, SSize>{};
	}

	BOOST_MULTI_HD constexpr auto halve() const {
		return layout_t<1, SSize>(*this, 0, 0, 0);
	}

	// [[deprecated("use two arg version")]] constexpr auto scale(size_type /*size*/) const {return *this;}
	constexpr auto scale(size_type /*num*/, size_type /*den*/) const { return *this; }

	//  friend constexpr auto operator!=(layout_t const& self, layout_t const& other) {return not(self == other);}
	friend BOOST_MULTI_HD constexpr auto operator==(layout_t const& self, layout_t const& other) {
		return 
			self.sub_ == other.sub_ &&
			self.stride_ == other.stride_ &&
			self.nelems_ == other.nelems_
		;
		// return std::tie(self.sub_, self.stride_, self.offset_, self.nelems_) == std::tie(other.sub_, other.stride_, other.offset_, other.nelems_);
	}

	friend BOOST_MULTI_HD constexpr auto operator!=(layout_t const& self, layout_t const& other) {
		return !(self==other);
		// return std::tie(self.sub_, self.stride_, self.offset_, self.nelems_) != std::tie(other.sub_, other.stride_, other.offset_, other.nelems_);
	}

	constexpr auto operator<(layout_t const& other) const -> bool {
		return std::tie(offset_, nelems_) < std::tie(other.offset_, other.nelems_);
	}

	BOOST_MULTI_HD constexpr auto rotate() const { return *this; }
	BOOST_MULTI_HD constexpr auto unrotate() const { return *this; }

	constexpr auto hull_size() const -> size_type { return num_elements(); }  // not in bytes
};

BOOST_MULTI_HD constexpr auto
operator*(layout_t<0>::index_extension const& extensions_0d, layout_t<0>::extents_type const& /*zero*/)
	-> layout_t<1>::extents_type {
	return layout_t<1>::extents_type{tuple<layout_t<0>::index_extension>{extensions_0d}};
}

BOOST_MULTI_HD constexpr auto operator*(extents_t<1> const& extensions_1d, extents_t<1> const& self) {
	using boost::multi::detail::get;
	return extents_t<2>({get<0>(extensions_1d.base()), get<0>(self.base())});
}

}  // end namespace boost::multi

namespace boost::multi::detail {

template<class Tuple>
struct convertible_tuple : Tuple {
	using Tuple::Tuple;
	BOOST_MULTI_HD explicit convertible_tuple(Tuple const& other)
	: Tuple(other) {}

 public:
	using array_type = std::array<std::ptrdiff_t, std::tuple_size_v<Tuple>>;
	auto to_array() const noexcept {
		return std::apply([](auto... elems) noexcept -> auto {
			return std::array<std::common_type_t<decltype(elems)...>, sizeof...(elems)>{{static_cast<multi::ssize_t>(elems)...}};
		},
						  static_cast<Tuple const&>(*this));
	}

	/*explicit*/ operator array_type() const& noexcept { return to_array(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	/*explicit*/ operator array_type() && noexcept { return to_array(); }      // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-stack-address"
#endif
	[[deprecated("This is here for nominal compatiblity with Boost.MultiArray, this would be a dangling conversion")]]
	operator std::ptrdiff_t const*() const&&;  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
											   /*{ return to_array().data(); }*/
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	template<std::size_t Index, std::enable_if_t<(Index < std::tuple_size_v<Tuple>), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend BOOST_MULTI_HD constexpr auto get(convertible_tuple const& self) -> std::tuple_element_t<Index, Tuple> {
		using std::get;
		return get<Index>(static_cast<Tuple const&>(self));
	}
};

template<class Array>
struct decaying_array : Array {
	using Array::Array;
	explicit decaying_array(Array const& other)
	: Array(other) {}

	[[deprecated("possible dangling conversion, use `std::array<T, D> p` instead of `auto* p`")]]
	constexpr operator std::ptrdiff_t const*() const { return Array::data(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

	template<std::size_t Index, std::enable_if_t<(Index < std::tuple_size_v<Array>), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	friend constexpr auto get(decaying_array const& self) -> std::tuple_element_t<Index, Array> {
		using std::get;
		return get<Index>(static_cast<Array const&>(self));
	}
};
}  // end namespace boost::multi::detail

template<class Tuple> struct std::tuple_size<boost::multi::detail::convertible_tuple<Tuple>> : std::integral_constant<std::size_t, std::tuple_size_v<Tuple>> {};  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal to define tuple size
template<class Array> struct std::tuple_size<boost::multi::detail::decaying_array<Array>> : std::integral_constant<std::size_t, std::tuple_size_v<Array>> {};     // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) normal to define tuple size

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
namespace std::ranges {  // NOLINT(cert-dcl58-cpp) to enable borrowed, nvcc needs namespace
template<>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::extents_t<1>::elements_t> = true;  // NOLINT(misc-definitions-in-headers)
}  // end namespace std::ranges
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_LAYOUT_HPP
