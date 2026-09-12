// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_EXTENTS_HPP
#define BOOST_MULTI_DETAIL_EXTENTS_HPP
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
namespace boost::multi::detail {
	/// Layout for strided arrays
	template <boost::multi::dimensionality_type D, typename SSize = multi::ssize_t>
	struct layout_t; 
}  // end namespace boost::multi::detail

namespace boost::multi::detail { template <class ...Ts> class tuple; }
// clang-format on

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost::multi {

// template<dimensionality_type D, typename SSize=multi::size_type> struct layout_t;

template<dimensionality_type D>
class extents_t;

/// `D`-tuple to store sizes of arrays
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

template<dimensionality_type D> class extents_t;

/// A structured cartesian product of extensions, it can be decomposed as a tuple-like into its cartesian factors. Used to determine the extents of array.
template<dimensionality_type D>
class extents_t {
	using base_ = boost::multi::detail::tuple_prepend_t<index_extension, typename extents_t<D - 1>::base_>;

	base_ impl_;  // NOLINT(misc-non-private-member-variables-in-classes) make private

 public:
	template<::boost::multi::dimensionality_type DD>
	using projection_type = std::tuple_element_t<static_cast<std::size_t>(DD), base_>;

	template<dimensionality_type> friend class extents_t;

	static constexpr dimensionality_type dimensionality = D;
	constexpr static dimensionality_type rank_v = D;

	/// Signed integer type to represent difference between indices (usually std::ptrdiff_t)
	using difference_type = index_extension::difference_type;

 private:
	using nelems_type = multi::index;

 public:
	/// A type to hold the size of the Cartesian product in the leading dimension
	using size_type = index_extension::size_type;

	/// Element type of the Cartesian product (a `D`-tuple of indices)
	using element = boost::multi::detail::tuple_prepend_t<index_extension::value_type, typename extents_t<D - 1>::element>;

	extents_t() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) intentionally trivial; default-init by design

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(multi::ssize_t size1)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : allow terse syntax
	: extents_t{index_extension{size1}} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: impl_{ext1} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 2, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2)
	: impl_{ext1, ext2} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 3, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3)
	: impl_{ext1, ext2, ext3} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 4, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4) noexcept
	: impl_{ext1, ext2, ext3, ext4} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 5, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4, index_extension ext5)
	: impl_{ext1, ext2, ext3, ext4, ext5} {}

	template<class T = void, std::enable_if_t<sizeof(T*) && D == 6, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr extents_t(index_extension ext1, index_extension ext2, index_extension ext3, index_extension ext4, index_extension ext5, index_extension ext6)
	: impl_{ext1, ext2, ext3, ext4, ext5, ext6} {}

	template<class... Ts, std::enable_if_t<sizeof...(Ts) == static_cast<std::size_t>(D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; allow terse syntax // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(detail::tuple<Ts...> const& exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: extents_t(exts, std::make_index_sequence<static_cast<std::size_t>(D)>()) {}

	template<class... Ts, std::enable_if_t<sizeof...(Ts) == static_cast<std::size_t>(D), int> = 0, class = decltype(base_{std::declval<::std::tuple<Ts...> >()})>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress noExplicitConstructor ; allow terse syntax // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(::std::tuple<Ts...> exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: impl_{std::move(exts)} {}

	template<
		class... Exts,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) TODO(correaa)
			(sizeof...(Exts) >= 2) && (sizeof...(Exts) == static_cast<std::size_t>(D))
			&& std::conjunction_v<std::is_convertible<Exts, index_extension>...>  // NOLINT(modernize-type-traits) not a fold-expr: MSVC 19.21 (VS2019 16.1) miscompiles `(... && ...)` here with C2059
			&& !std::conjunction_v<std::is_integral<Exts>...>,  // NOLINT(modernize-type-traits) for C++20
			int> = 0
	>
	BOOST_MULTI_HD constexpr extents_t(Exts... exts)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) allow terse syntax
	: impl_{index_extension(exts)...} {}

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
	: impl_{index_extension(static_cast<index_extension::size_type>(exts))...} {}

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
	: impl_{index_extension(exts)...} {}

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

	BOOST_MULTI_HD constexpr extents_t(index_extension const& ext, typename detail::layout_t<D - 1>::extents_type const& other)
	: extents_t(multi::detail::ht_tuple(ext, other.base())) {}

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return impl_; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return impl_; }

	friend constexpr auto operator*(index_extension const& ext, extents_t const& self) -> extents_t<D + 1> {
		// return extents_t<D + 1>(tuple(extension, self.base()));
		return extents_t<D + 1>(ext, self);
	}

	friend BOOST_MULTI_HD auto operator==(extents_t const& self, extents_t const& other) { return self.base() == other.base(); }
	friend BOOST_MULTI_HD auto operator!=(extents_t const& self, extents_t const& other) { return self.base() != other.base(); }

	/// Type to store an index in the leading dimension
	using index        = multi::index;
	/// A tuple type that allows storing `D` indices to locate an element in the structured Cartesian product
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
		return extents_t<D - 1>{impl_.tail()};
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

	constexpr explicit operator bool() const { return !detail::layout_t<D>{*this}.empty(); }  // TODO(correaa) simplifty algorithm

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto to_linear(index const& idx, Indices const&... rest) const {
		auto const sub_extensions = extents_t<D - 1>{this->base().tail()};
		return (idx * sub_extensions.num_elements()) + sub_extensions.to_linear(rest...);
	}

	template<class... Indices>
	BOOST_MULTI_HD constexpr auto operator()(index idx, Indices... rest) const { return to_linear(idx, rest...); }

 private:
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
				return cursor_t<typename multi::detail::layout_t<std::tuple_size_v<Before> + 1>::indexes, DD - 1> (
					apply([n] (auto... idxs) -> auto {return detail::mk_tuple(idxs..., n);}, bef_)
				);
			} else {
				return apply([n] (auto... idxs) -> auto {return detail::mk_tuple(idxs..., n);}, bef_);
			}
		}
	};

 public:
	/// Indexable cursor, pointer‐like objects that multidimensionally indexable (type usually returned by `.home()`)
	using cursor = cursor_t<tuple<>, D>;

	/// Returns a cursor to the home (e.g. top-left) element
	static auto home() -> cursor { return {}; }

	/// Random-access iterator in the leading dimension of the extents object, in general they dereference an extents of lower dimension, for D == 1 to an index value
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
		return impl_[idx];
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

	/// Range with the elements of an extents (Cartesian product of extent integer ranges)
	class elements_t {
		extents_t xs_;
		explicit constexpr elements_t(extents_t const& exts) : xs_{exts} {}

		friend class extents_t;

	 public:
		using difference_type = extents_t::difference_type;

		/// Random-access iterator in the leading dimension, resulting from `begin()/end()`
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

		 private:
			template<class CUT>
			class mk_tup {
				CUT cu_;

			 public:
				constexpr explicit mk_tup(CUT current) : cu_{current} {}
				template<class... Ts>
				constexpr auto operator()(Ts... idxs) const { return detail::mk_tuple(cu_, idxs...); }
			};

		 public:
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
				xs_.base().head().begin(),
				extents_t<D - 1>{xs_.base().tail()}.elements().begin(),
				extents_t<D - 1>{xs_.base().tail()}.elements().begin(),
				extents_t<D - 1>{xs_.base().tail()}.elements().end(),
			};
		}

		constexpr auto end() const {
			return iterator{
				xs_.base().head().end(),
				extents_t<D - 1>{xs_.base().tail()}.elements().begin(),
				extents_t<D - 1>{xs_.base().tail()}.elements().begin(),
				extents_t<D - 1>{xs_.base().tail()}.elements().end(),
			};
		}

		BOOST_MULTI_HD constexpr auto operator[](index idx) const { return begin()[idx]; }

		BOOST_MULTI_HD constexpr auto  size() const noexcept { return xs_.num_elements(); }
		BOOST_MULTI_HD constexpr auto ssize() const noexcept { return this->size(); }
		BOOST_MULTI_HD constexpr auto usize() const noexcept { return static_cast<std::size_t>(xs_.num_elements()); }
	};

	/// yields a random‐access range with all the elements of the extents cartesian product object
	constexpr auto elements() const { return elements_t{*this}; }

	template<class Func>
	BOOST_MULTI_HD constexpr auto element_transformed(Func fun) const { return [fun](auto const&... idxs) -> decltype(auto) { return fun(detail::mk_tuple(idxs...)); } ^(*this); }

	BOOST_MULTI_HD constexpr auto               extension() const { return this->get<0>(); }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0
	[[nodiscard]] BOOST_MULTI_HD constexpr auto extent() const { return this->get<0>(); }     // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0

	BOOST_MULTI_HD constexpr auto size() const noexcept { return this->get<0>().size(); }
	BOOST_MULTI_HD constexpr auto sizes() const {
		return this->apply([](auto const&... exts) -> auto { return multi::detail::mk_tuple(exts.size()...); });
	}

	// this implementation is incorrect, use restriction
	// /// creates an extents object with the indices rotated to the left
	// constexpr auto rotate() const {
	// 	this->apply([](auto const& head, auto const&... rest) -> extents_t { return extents_t(rest..., head); });
	// }

	// this implementation is incorrect, use restriction
	// /// creates an extents object with the indices rotated to the right
	// constexpr auto unrotate() const {
	// 	this->apply([](auto const&... rest, auto const& tail) -> extents_t { return extents_t(tail, rest...); });
	// }

	// this implementation is incorrect, use restriction
	// /// creates an extents object with the first two dimensions transposed
	// constexpr auto transpose() const {
	// 	return this->apply([](auto const& head1, auto const& head2, auto const&... rest) -> extents_t { return extents_t(head2, head1, rest...); });
	// }

	// /// creates an extents object with the first two dimensions transposed
	// constexpr auto transposed() const { return transpose(); }

	[[deprecated]] BOOST_MULTI_HD constexpr auto extensions() const {
		using std::apply;
		return apply([](auto... sizes) -> extents_t { return extents_t(sizes...); }, sizes());
	}
	BOOST_MULTI_HD constexpr auto extents() const {
		using std::apply;
		return apply([](auto... sizes) -> extents_t { return extents_t(sizes...); }, sizes());
	}
	/// A `D`-tuple to hold the sizes of the structured Cartesian product
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
	: impl_{boost::multi::detail::get<I>(tup)...} {}

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
	// friend constexpr auto get(extents_t const& self) -> typename std::tuple_element_t<Index, base_> {
	friend constexpr auto get(extents_t const& self) -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// friend constexpr auto get(extents_t const& self) -> typename std::tuple_element_t<Index, base_> {
	friend constexpr auto get(extents_t& self) -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// friend constexpr auto get(extents_t const& self) -> typename std::tuple_element_t<Index, base_> {
	friend constexpr auto get(extents_t&& self) -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(std::move(self).base());
	}


	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// constexpr auto get() const -> std::tuple_element_t<Index, base_> {
	constexpr auto get() const& -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// constexpr auto get() const -> std::tuple_element_t<Index, base_> {
	constexpr auto get() & -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < D), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// constexpr auto get() const -> std::tuple_element_t<Index, base_> {
	constexpr auto get() && -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(std::move(this->base()));
	}

	template<class F>
	constexpr auto apply(F&& fun) const -> decltype(auto) {
		return std::apply(std::forward<F>(fun), this->base());
	}
};

template<> class extents_t<0> {
	tuple<> impl_;  // NOLINT(misc-non-private-member-variables-in-classes) make private
	using base_ = tuple<>;

 public:
	template<::boost::multi::dimensionality_type DD>
	using projection_type = std::tuple_element_t<static_cast<std::size_t>(DD), base_>;

	template<dimensionality_type> friend class extents_t;

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
	: impl_{tup} {}

	extents_t() = default;

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return impl_; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return impl_; }

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

template<> class extents_t<1> {
	tuple<multi::index_extension> impl_;  // NOLINT(misc-non-private-member-variables-in-classes) make private
	using base_ = tuple<multi::index_extension>;

 public:
	template<::boost::multi::dimensionality_type DD>
	using projection_type = std::tuple_element_t<static_cast<std::size_t>(DD), base_>;

	template<dimensionality_type> friend class extents_t;

	static constexpr auto dimensionality = 1;  // TODO(correaa): consider deprecation

	constexpr static dimensionality_type rank_v = 1;

	using size_type = multi::index_extension::size_type;
	using difference_type = multi::index_extension::difference_type;
	using element = tuple<multi::index_extension::value_type>;
	using index = multi::index;
	using sizes_type = tuple<size_type>;

	[[deprecated("use .extent()")]] constexpr auto extension() const { using std::get; return get<0>(impl_); }
	[[nodiscard]] constexpr auto extent() const { using std::get; return get<0>(impl_); }

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
		return elements_t{get<0>(impl_)};
	}

	using nelems_type = index;

	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax (compatible with std::vector(int) constructor
	BOOST_MULTI_HD constexpr extents_t(multi::ssize_t size)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: impl_(multi::index_extension{0, size}) {}

	template<class T1>
	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int>  // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(tuple<T1> extensions)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: impl_{static_cast<multi::index_extension>(extensions.head())} {}

	// cppcheck-suppress noExplicitConstructor ; to allow passing tuple<int, int> // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(multi::index_extension const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: impl_{other} {}

	BOOST_MULTI_HD constexpr explicit extents_t(base_ tup)
	: impl_{tup} {}

	template<class OtherExtents,
		decltype( multi::detail::implicit_cast<multi::index_extension>(OtherExtents{}.extent()) )* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	BOOST_MULTI_HD constexpr extents_t(OtherExtents const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: impl_{other.extent()} {}

	extents_t() = default;

	BOOST_MULTI_HD constexpr auto base() const& -> base_ const& { return impl_; }
	BOOST_MULTI_HD constexpr auto base() & -> base_& { return impl_; }

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
	// constexpr auto get() const -> std::tuple_element_t<Index, base_> {  // by value, to match the other extents_t<D> specializations (structured bindings)
	constexpr auto get() const -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(this->base());
	}

	template<std::size_t Index, std::enable_if_t<(Index < 1), int> = 0>                       // NOLINT(modernize-use-constraints) TODO(correaa)
	// friend constexpr auto get(extents_t const& self) -> std::tuple_element_t<Index, base_> {  // by value, to match the other extents_t<D> specializations (structured bindings)
	friend constexpr auto get(extents_t const& self) -> multi::index_extension {  // by value (== tuple_element_t<Index, base_>), spelled without a library-trait alias so the signature has no built-in trait (g++-15 / nvcc)
		using boost::multi::detail::get;
		return get<Index>(self.base());
	}
};

template<dimensionality_type D> using iextensions = extents_t<D>;

template<dimensionality_type D> using extensions_t [[deprecated("use extents_t")]] = extents_t<D>;

// template<boost::multi::dimensionality_type D>
// constexpr auto array_size_impl(boost::multi::extents_t<D> const&)
// 	-> std::integral_constant<std::size_t, static_cast<std::size_t>(D)>;

extents_t(multi::ssize_t) -> extents_t<1>;
/// deduces an extents dimension from the number of size arguments
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
	using type = void;  // TODO(correaa) should be undefined!
};

template<std::size_t Index, boost::multi::dimensionality_type D>
struct std::tuple_element<Index, boost::multi::extents_t<D>> {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) to implement structured binding
	// using type = typename std::tuple_element_t<Index, typename boost::multi::extents_t<D>::base_>;
	using type = typename ::boost::multi::extents_t<D>::template projection_type<Index>;
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

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<typename Fn, boost::multi::dimensionality_type D>
constexpr auto
apply(Fn&& fun, boost::multi::extents_t<D> const& exts) noexcept -> decltype(auto) {  // NOLINT(cert-dcl58-cpp,bugprone-std-namespace-modification) workaround
	return exts.apply(std::forward<Fn>(fun));
}

}  // end namespace std

#if defined(__NVCC__) || (defined(__NVCOMPILER) && (__NVCOMPILER_MAJOR__ <= 23))
namespace boost::multi {

// namespace-scope `get` for extents_t: the ADL target that the `detail::tuple` base used to
// provide before extents_t switched from inheritance to composition.
//  - nvcc/EDG does not find the in-class hidden-friend `get<Index>` template through ADL when
//    called with explicit template arguments (`get<0>(exts)`);
//  - being callable qualified (`multi::get<N>(exts)`) it also works in trailing-return SFINAE
//    and in C++17, unlike the `using std::get; get<N>(exts)` two-step;
//  - `std::get` cannot serve this role: adding `get` overloads to namespace std is UB and
//    collides with the libc++ tuple machinery.
template<std::size_t N, ::boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D> const& tp)
	-> decltype(tp.template get<N>()) {
	return tp.template get<N>();
}

template<std::size_t N, ::boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D>& tp)
	-> decltype(tp.template get<N>()) {
	return tp.template get<N>();
}

template<std::size_t N, boost::multi::dimensionality_type D>
constexpr auto get(::boost::multi::extents_t<D>&& tp)
	-> decltype(std::move(tp).template get<N>()) {
	return std::move(tp).template get<N>();
}

}  // end namespace boost::multi
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
namespace std::ranges {  // NOLINT(cert-dcl58-cpp) to enable borrowed, nvcc needs namespace
template<>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::extents_t<1>::elements_t> = true;  // NOLINT(misc-definitions-in-headers)
}  // end namespace std::ranges
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_EXTENTS_HPP
