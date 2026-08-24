// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP
// #pragma once

#include "boost/multi/detail/tuple_zip.hpp"
#include "boost/multi/utility.hpp"  // IWYU pragma: export

#include <cmath>
#include <type_traits>

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#include <ranges>  // IWYU pragma: keep
#include <vector>  // for .to conversion
#endif

#if (__cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L))
#if __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4623)  // assignment operator was implicitly defined as deleted
#pragma warning(disable : 4626)  // assignment operator was implicitly defined as deleted
#pragma warning(disable : 4625)  // copy constructor was implicitly defined as deleted
#endif

namespace boost::multi {

/// When specialized to `true` for an element type, the library
/// treats that type as if it were trivially default-constructible *and* trivially
/// destructible. This is a performance opt-in for types that are
/// technically non-trivial under the language rules (such as std::complex).
template<class Element>
inline constexpr bool force_element_trivial = false;

/// When specialized to `true`, skip the per-element destructor loop on array teardown.
template<class Element>
inline constexpr bool force_element_trivial_destruction = force_element_trivial<Element>;

/// When specialized to `true`, skip the per-element default-construction loop on array allocation.
template<class Element>
inline constexpr bool force_element_trivial_default_construction = force_element_trivial<Element>;

/// Opt-in specializations for `std::complex<T>`, enabled by defining
/// `_BOOST_MULTI_FORCE_TRIVIAL_STD_COMPLEX` at compile time, which treats std::complex
/// types as trivially constructible/destructible.
#ifdef _BOOST_MULTI_FORCE_TRIVIAL_STD_COMPLEX
template<class T>
inline constexpr bool force_element_trivial<std::complex<T>> = std::is_trivially_copyable_v<T> && std::is_trivially_default_constructible_v<T>;  // std::is_trivial_v deprecated in C++26

template<class T>
inline constexpr bool force_element_trivial_destruction<std::complex<T>> = std::is_trivially_default_constructible_v<T>;

template<class T>
inline constexpr bool force_element_trivial_default_construction<std::complex<T>> = std::is_trivially_destructible_v<T>;

template<> inline constexpr bool force_element_trivial<std::complex<double>>                      = true;
template<> inline constexpr bool force_element_trivial_default_construction<std::complex<double>> = true;
template<> inline constexpr bool force_element_trivial_destruction<std::complex<double>>          = true;

template<> inline constexpr bool force_element_trivial<std::complex<float>>                      = true;
template<> inline constexpr bool force_element_trivial_default_construction<std::complex<float>> = true;
template<> inline constexpr bool force_element_trivial_destruction<std::complex<float>>          = true;
#endif

}  // end namespace boost::multi

#include "boost/multi/detail/adl.hpp"  // TODO(correaa) remove instantiation of force_element_trivial in this header
#include "boost/multi/detail/config/ASSERT.hpp"
#include "boost/multi/detail/layout.hpp"  // IWYU pragma: export
// #include "boost/multi/detail/memory.hpp"          // for pointer_traits
#include "boost/multi/detail/operators.hpp"       // for random_iterable
#include "boost/multi/detail/pointer_traits.hpp"  // IWYU pragma: export
#include "boost/multi/detail/serialization.hpp"
#include "boost/multi/detail/types.hpp"  // for dimensionality_type  // IWYU pragma: export

#include <algorithm>  // for copy_n
#include <array>      // for std::array
#include <cstring>    // for std::memset in reinterpret_cast
// #include <functional>  // for std::invoke
#include <iterator>  // for std::next
#include <memory>    // for std::pointer_traits
#include <new>       // for std::launder

#if __has_include(<span>)
#if !defined(_MSVC_LANG) || (_MSVC_LANG > 202002L)
#include <span>
#endif
#if defined(__cpp_lib_span) && __cpp_lib_span >= 202002L && !defined(_MSVC_LANG)
#define BOOST_MULTI_HAS_SPAN
#endif
#endif

#include <utility>  // for forward

#ifdef __NVCC__
#define BOOST_MULTI_FRIEND_CONSTEXPR template<class = void> friend constexpr  // workaround nvcc
#else
#define BOOST_MULTI_FRIEND_CONSTEXPR friend constexpr
#endif

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

#if defined(__has_cpp_attribute) && __has_cpp_attribute(gnu::no_dangling)
#define BOOST_MULTI_NO_DANGLING [[gnu::no_dangling]]
#else
#define BOOST_MULTI_NO_DANGLING
#endif

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#define BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_PUSH() \
	_Pragma("clang diagnostic push")                   \
		_Pragma("clang diagnostic ignored \"-Wunsafe-buffer-usage\"")
#else
#define BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_PUSH()
#endif

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#define BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_POP() _Pragma("clang diagnostic pop")
#else
#define BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_POP()
#endif

namespace boost::multi {

/// Read-only `D`-dimensional subarray-reference (view into a subarray)
///
/// Provides the same interface as `subarray` but prevents modification of the referenced elements.
/// Has reference semantics: cannot be rebound, assignments are deep, and size is immutable.
///
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam ElementPtr Pointer-like type to const elements (default `T const*`)
/// @tparam Layout Layout type describing strides and extensions
template<typename T, dimensionality_type D, typename ElementPtr = T const*, class Layout = layout_t<D>>
class const_subarray;

namespace detail {
template<class T, dimensionality_type D, class... Ts>
auto is_const_subarray_aux(multi::const_subarray<T, D, Ts...> const&) -> std::true_type;
auto is_const_subarray_aux(...) -> std::false_type;

template<class T> struct is_const_subarray : decltype(detail::is_const_subarray_aux(std::declval<T>())){};
template<class T>
constexpr bool is_const_subarray_v = is_const_subarray<T>::value;
}  // end namespace detail

/// Mutable `D`-dimensional view into part or all of an array
///
/// Represents a subregion of a larger array without owning the elements.
/// Has reference semantics: cannot be rebound, assignments are deep, and size is immutable.
/// Invalidated if the originating array is destroyed or resized.
///
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam ElementPtr Pointer-like type to the elements (default `T*`)
/// @tparam Layout type describing strides and extensions
template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>>
class subarray;

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>>
class move_subarray;

namespace detail {
template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
constexpr auto is_subarray_aux(const_subarray<T, D, ElementPtr, Layout> const&) -> std::true_type;
constexpr auto is_subarray_aux(...) -> std::false_type;

template<class A> struct is_subarray : decltype(is_subarray_aux(std::declval<A>())){};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<dimensionality_type D>
struct of_dim {
	template<typename T, class ElementPtr, class Layout>
	static constexpr auto is_subarray_of_dim_aux(subarray<T, D, ElementPtr, Layout> const&) -> std::true_type;
	static constexpr auto is_subarray_of_dim_aux(...) -> std::false_type;

	template<class A> struct is_subarray_of_dim : decltype(is_subarray_of_dim_aux(std::declval<A>())){};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
};
}  // namespace detail

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace detail {
template<class Element, dimensionality_type D, typename ElementPtr, bool IsConst = false, bool IsMove = false, typename Stride = typename std::iterator_traits<ElementPtr>::difference_type, class SubLayout = layout_t<D - 1>>
struct array_iterator;
}  // end namespace detail

template<class Element, dimensionality_type D, typename ElementPtr, bool IsConst = false, bool IsMove = false, typename Stride = typename std::iterator_traits<ElementPtr>::difference_type, class SubLayout = layout_t<D - 1>>
using array_iterator [[deprecated]] = typename detail::array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride, SubLayout>;

namespace detail {
template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D, std::make_signed_t<typename std::pointer_traits<ElementPtr>::size_type>>>
struct array_types : private Layout {  // cppcheck-suppress syntaxError ; false positive in cppcheck
	using element                                      = T;
	using element_type [[deprecated("use ::element")]] = element;  // this follows more closely https://en.cppreference.com/w/cpp/memory/pointer_traits

	using element_ptr       = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element const>;

	/// Pointer-like type that produces an moved element (r-value)
	using element_move_ptr = multi::move_ptr<element, element_ptr>;

	using element_ref = typename std::iterator_traits<element_ptr>::reference;

	using layout_type = Layout;

	using layout_type::rank_v;

	/// Integer type to store dimensionality information (e.g. 1D, 2D, 3D)
	using dimensionality_type = typename layout_type::dimensionality_type;
	using layout_type::dimensionality;

	using layout_type::stride;
	using typename layout_type::stride_type;

	using layout_type::num_elements;
	using layout_type::offset;

	/// Type to store an index in the leading dimension
	using index = typename layout_type::index;

	using typename layout_type::index_extension;

	/// Type that represents a range of indices
	using index_range = typename layout_type::index_range;  // re-export publicly: array_types inherits Layout privately, so MSVC otherwise sees index_range as inaccessible (C2247) in derived subarray classes

	using typename layout_type::strides_type;

	/// returns the layout internal strides of an array as a tuple
	BOOST_MULTI_HD constexpr auto strides() const {  // cppcheck-suppress [functionStatic,duplInheritedMember];
		return layout_type::strides();
	}

	using typename layout_type::difference_type;

	/// Type to hold the size of the array in the leading dimension
	using size_type = typename layout_type::size_type;
	// using layout_type::size;
	/// returns the size of the array in the leading dimension
	BOOST_MULTI_HD constexpr auto size() const noexcept -> size_type { return layout_type::size(); }

	using layout_type::nelems;

	using layout_type::extension;
	using layout_type::extent;

	/// A type to store the extent of an array (the range of valid indices in the leading dimension), returned from `.extent()`.
	using extent_type = typename layout_type::extent_type;

	/// (deprecated) use `extent_type`
	using extension_type [[deprecated("use extent_type")]] = extent_type;  // NOLINT  ; old spelling kept for compatibility

	using layout_type::extensions;
	using layout_type::extents;

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
	using typename layout_type::extensions_type;  // TODO(correaa) remove
#ifdef _MSC_VER
#pragma warning(pop)
#endif

	// using typename layout_type::extents_type;
	/// A type that stores the extents of the array or subarray (returned from `.extents()` and used for `array` constructors)
	using extents_type = typename layout_type::extents_type;

	/// A tuple type that allows storing `D` indices to locate an element in the array
	using indices_type = typename layout_type::extents_type::indices_type;

	/// Returns the index extensions (structured cartesian product of half-open ranges) for all dimensions as an `extents_type`
	/// (`extents_t<D>`), a tuple of `D` `index_extension` values each encoding `[first, last)`.
	/// The result can be passed directly to array constructors or compared for shape equality.
	/// Prefer this over the deprecated `extensions()`.
	[[nodiscard]] BOOST_MULTI_HD constexpr auto extents() const -> extents_type { return static_cast<layout_type const&>(*this).extents(); }  // cppcheck-suppress duplInheritedMember;

	/// Deprecated, prefer `extents()`.
	[[deprecated("use .extents()")]] BOOST_MULTI_HD constexpr auto extensions() const -> extents_type { return static_cast<layout_type const&>(*this).extents(); }  // cppcheck-suppress duplInheritedMember;

	using layout_type::empty;
	using layout_type::is_empty;

	using layout_type::sub;

	using layout_type::sizes;
	using typename layout_type::sizes_type;

	using typename layout_type::indexes;

	// [[deprecated("This is for compatiblity with Boost.MultiArray, you can use `rank` member type or `dimensionality` static member variable")]]
	// static constexpr auto num_dimensions() { return dimensionality; }

 private:
	[[deprecated("This is for compatiblity with Boost.MultiArray, you can use `offsets` member function")]]
	// NOLINTNEXTLINE(readability-identifier-naming) TODO(correaa) to remove
	auto index_bases() const -> std::ptrdiff_t const*;  // = delete;  this function is not implemented, it can give a linker error

 public:
	[[deprecated("This is for compatiblity with Boost.MultiArray, you can use `offsets` member function")]]
	constexpr auto shape() const { return detail::convertible_tuple<decltype(this->sizes())>(this->sizes()); }

	[[deprecated]] auto is_compact() const { return this->layout().is_compact(); }

 private:
	constexpr auto layout_mutable() -> layout_type& { return static_cast<layout_type&>(*this); }  // NOLINT(readability-identifier-naming)

	template<typename, ::boost::multi::dimensionality_type, class> friend struct ::boost::multi::array;
	template<typename, ::boost::multi::dimensionality_type, class> friend struct ::boost::multi::dynamic_array;

 public:
	/// Array value after evaluation through the first index, an object of lower dimension, `multi::array<T, D ‐ 1, P>` or, for `D == 1`, `std::pointer_traits<P>::element_type` (usually `T`)
	using value_type = typename std::conditional_t<
		(D > 1),
		array<element, D - 1, typename multi::pointer_traits<element_ptr>::default_allocator_type>,
		element>;

	using reference = typename std::conditional_t<
		(D > 1),
		subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_ptr>::reference>;

	using const_reference = typename std::conditional_t<
		(D > 1),
		const_subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_const_ptr>::reference>;

	/// returns the base pointer of the array (arithmetic base of the layout, generally the first element)
	BOOST_MULTI_HD constexpr auto base() const -> element_const_ptr { return base_; }  // cppcheck-suppress duplInheritedMember ; to overwrite

	/// returns the base const-pointer of the array (arithmetic base of the layout, generally the first element)
	BOOST_MULTI_HD constexpr auto cbase() const -> element_const_ptr { return base_; }

	/// returns the layout of the array
	BOOST_MULTI_HD constexpr auto layout() const -> layout_type const& { return *this; }

	BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_PUSH()
	// cppcheck-suppress duplInheritedMember ; to overwrite
	[[deprecated("for compatibility with BMA, use .base()")]]
	constexpr auto origin() const& -> decltype(auto) { return base_ + Layout::origin(); }  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	BOOST_MULTI_IGNORED_UNSAFE_BUFFER_USAGE_POP()

 protected:
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // warning C4820:  '7' bytes padding added after data member 'boost::multi::array_types<T,2,ElementPtr,Layout>::base_' [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\test\array_fancyref.cpp.x.vcxproj]
#endif
	BOOST_MULTI_NO_UNIQUE_ADDRESS
	element_ptr base_;  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes,misc-non-private-member-variables-in-classes) : TODO(correaa) try to make it private, [static_]array needs mutation
#ifdef _MSC_VER
#pragma warning(pop)
#endif

	template<class, ::boost::multi::dimensionality_type, typename, bool, bool, typename, class> friend struct detail::array_iterator;

	// using derived = subarray<T, D, ElementPtr, Layout>;
	BOOST_MULTI_HD constexpr explicit array_types(std::nullptr_t) : Layout{}, base_(nullptr) {}

 public:
	array_types() = default;  // cppcheck-suppress uninitMemberVar ; base_ not initialized

	BOOST_MULTI_HD constexpr array_types(layout_type const& lyt, element_ptr data)
	: Layout{lyt}, base_{std::move(data)} {}

 protected:
	template<
		class ArrayTypes,
		typename                                                                                      = std::enable_if_t<!std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>,
		decltype(multi::detail::explicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr>
	// underlying pointers are explicitly convertible
	BOOST_MULTI_HD constexpr explicit array_types(ArrayTypes const& other)
	: Layout{other.layout()}, base_{other.base_} {}

	template<
		class ArrayTypes,
		typename                                                                                      = std::enable_if_t<!std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>,
		decltype(multi::detail::implicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointers are implicitly convertible
	BOOST_MULTI_HD constexpr /*implt*/ array_types(ArrayTypes const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : inherit behavior of underlying pointer
	: Layout{other.layout()}, base_{other.base_} {}

	template<
		typename ElementPtr2,
		typename = decltype(Layout{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().layout()}),
		typename = decltype(element_ptr{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().base_})>
	BOOST_MULTI_HD constexpr explicit array_types(array_types<T, D, ElementPtr2, Layout> const& other)
	: Layout{other.layout()}, base_{other.base_} {}

	template<class T2, ::boost::multi::dimensionality_type D2, class E2, class L2> friend struct array_types;
};
}  // end namespace detail

#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace detail {
/// Pointer-like type to a `subarray`, providing an extra level of indirection
///
/// Corresponds to `subarray` the way a raw pointer corresponds to a reference.
/// Can be rebound (unlike `subarray`) and supports pointer arithmetic.
///
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam ElementPtr Pointer-like type to the elements
/// @tparam Layout Layout type describing strides and extensions
/// @tparam IsConst Whether the pointed-to subarray is const
template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout, bool IsConst>
struct subarray_ptr;
}  // end namespace detail

// template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = multi::layout_t<D>>
// using const_subarray_ptr = detail::subarray_ptr<T, D, ElementPtr, Layout, true>;

namespace detail {
template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = multi::layout_t<D>, bool IsConst = false>
struct subarray_ptr  // : to allow mixin CRTP
: boost::multi::iterator_facade<
	  subarray_ptr<T, D, ElementPtr, Layout, IsConst>, void, std::random_access_iterator_tag,
	  subarray<T, D, ElementPtr, Layout> const&, typename Layout::difference_type> {

 private:
	Layout layout_;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  //'boost::multi::subarray_ptr<double,1,fancy::ptr<double>,boost::multi::layout_t<1,boost::multi::size_type>,true>': '7' bytes padding added after data member 'boost::multi::subarray_ptr<double,1,fancy::ptr<double>,boost::multi::layout_t<1,boost::multi::size_type>,true>::base_'
#endif

	ElementPtr base_;
	// typename std::iterator_traits<ElementPtr>::difference_type offset_;  // TODO(correaa) revive under a macro, unused for now

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

 public:
	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;
	template<typename, multi::dimensionality_type, typename, bool, bool, typename, class> friend struct detail::array_iterator;

	// ~subarray_ptr() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	using pointer         = subarray<T, D, ElementPtr, Layout>*;
	using element         = typename subarray<T, D, ElementPtr, Layout>::decay_type;
	using difference_type = typename Layout::difference_type;

	using value_type = element;

	using reference = std::conditional_t<
		IsConst,
		const_subarray<T, D, ElementPtr, Layout>,
		subarray<T, D, ElementPtr, Layout>>;

	using iterator_category = std::random_access_iterator_tag;

	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr subarray_ptr(std::nullptr_t nil) : layout_{}, base_{nil} /*, offset_{0}*/ {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) terse syntax and functionality by default

	subarray_ptr() = default;  // cppcheck-suppress uninitMemberVar ; base_ is not initialized

	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;

	BOOST_MULTI_HD constexpr subarray_ptr(typename reference::element_ptr base, layout_t<reference::dimensionality - 1> lyt) : layout_{lyt}, base_{base} /*, offset_{0}*/ {}

	template<bool OtherIsConst, std::enable_if_t<!OtherIsConst, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// cppcheck-suppress noExplicitConstructor ; see below
	BOOST_MULTI_HD constexpr /*mplct*/ subarray_ptr(subarray_ptr<T, D, ElementPtr, Layout, OtherIsConst> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : propagate implicitness of pointer
	: layout_{other.layout_}, base_{other.base_} /*, offset_{other.offset_}*/ {}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherLayout, bool OtherIsConst,
		decltype(multi::detail::implicit_cast<typename reference::element_ptr>(std::declval<OtherEPtr>()))* = nullptr  // propagate implicitness of pointer
		>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr /*mplct*/ subarray_ptr(subarray_ptr<OtherT, OtherD, OtherEPtr, OtherLayout, OtherIsConst> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: layout_{other.layout_}, base_{other.base_} /*, offset_{other.offset_}*/ {}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherLayout, bool OtherIsConst,
		decltype(multi::detail::explicit_cast<typename reference::element_ptr>(std::declval<OtherEPtr>()))* = nullptr  // propagate implicitness of pointer
		>
	BOOST_MULTI_HD constexpr explicit subarray_ptr(subarray_ptr<OtherT, OtherD, OtherEPtr, OtherLayout, OtherIsConst> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: layout_{other.layout_}, base_{other.base_} /*, offset_{other.offset_}*/ {}

	template<
		class ElementPtr2,
		std::enable_if_t<std::is_same_v<ElementPtr2, ElementPtr> && (D == 0), int> = 0  // NOLINT(modernize-use-constraints) for C++20
		>
	BOOST_MULTI_HD constexpr explicit subarray_ptr(ElementPtr2 const& other) : layout_{}, base_{other} /*, offset_{0}*/ {}

	BOOST_MULTI_HD constexpr explicit operator bool() const { return static_cast<bool>(base()); }

	// cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto operator*() const -> reference { return reference(layout_, base_); }

	BOOST_MULTI_HD constexpr auto operator->() const {
		class proxy {
			reference ref_;

		 public:
			BOOST_MULTI_HD constexpr explicit proxy(reference&& ref) : ref_{std::move(ref)} {}
			BOOST_MULTI_HD constexpr auto operator->() && -> reference* { return std::addressof(this->ref_); }
		};
		return proxy{operator*()};
	}

	using raw_pointer = reference*;  ///< Type aliases for Thrust introspection

	BOOST_MULTI_HD constexpr auto get() const { return reinterpret_cast<reference* const&>(*this); }  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)

	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> reference { return *(*this + n); }

	BOOST_MULTI_HD constexpr auto operator<(subarray_ptr const& other) const -> bool { return distance_to(other) > 0; }

	BOOST_MULTI_HD constexpr subarray_ptr(typename reference::element_ptr base, Layout const& lyt) : layout_{lyt}, base_{std::move(base)} /*, offset_{0}*/ {}

	template<typename, multi::dimensionality_type, typename, class> friend class const_subarray;

	BOOST_MULTI_HD constexpr auto base() const -> typename reference::element_ptr { return base_; }  // cppcheck-suppress returnByReference;

	template<class OtherSubarrayPtr, std::enable_if_t<!std::is_base_of_v<subarray_ptr, OtherSubarrayPtr>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator==(OtherSubarrayPtr const& other) const
		-> decltype(base_ == other.base_ && layout_ == other.layout_) {
		return base_ == other.base_ && layout_ == other.layout_;
	}

	template<class OtherSubarrayPtr, std::enable_if_t<!std::is_base_of_v<subarray_ptr, OtherSubarrayPtr>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator!=(OtherSubarrayPtr const& other) const
		-> decltype(base_ != other.base_ || layout_ != other.layout_) {
		return base_ != other.base_ || layout_ != other.layout_;
	}

	constexpr auto operator==(subarray_ptr const& other) const -> bool {
		return base_ == other.base_ && layout_ == other.layout_;
	}

	constexpr auto operator!=(subarray_ptr const& other) const -> bool {
		return base_ != other.base_ || layout_ != other.layout_;
	}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherL, bool OtherIsConst,
		std::enable_if_t<!std::is_base_of_v<subarray_ptr, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst>>, int> = 0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
		>
	friend BOOST_MULTI_HD constexpr auto operator==(subarray_ptr const& self, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst> const& other) -> bool {
		BOOST_MULTI_ASSERT((!self || !other) || (self->layout() == other->layout()));  // comparing array ptrs of different provenance is undefined
		return self->base() == other->base();
	}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherL, bool OtherIsConst,
		std::enable_if_t<!std::is_base_of_v<subarray_ptr, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst>>, int> = 0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
		>
	friend BOOST_MULTI_HD constexpr auto operator!=(subarray_ptr const& self, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst> const& other) -> bool {
		BOOST_MULTI_ASSERT((!self || !other) || (self->layout() == other->layout()));  // comparing array ptrs of different provenance is undefined
		return self->base() != other->base();
	}

 protected:
	BOOST_MULTI_HD constexpr void increment() { base_ += layout_.nelems(); }
	BOOST_MULTI_HD constexpr void decrement() { base_ -= layout_.nelems(); }

	BOOST_MULTI_HD constexpr void advance(difference_type n) { base_ += layout_.nelems() * n; }
	BOOST_MULTI_HD constexpr auto distance_to(subarray_ptr const& other) const -> difference_type {
		BOOST_MULTI_ASSERT(layout_.nelems() == other.layout_.nelems());
		// assert( Ref::nelems() == other.Ref::nelems() && Ref::nelems() != 0 );
		// assert( (other.base() - base())%Ref::nelems() == 0);
		BOOST_MULTI_ASSERT(layout_ == other.layout_);
		return (other.base_ - base_) / layout_.nelems();
	}

 public:
	BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> subarray_ptr& {
		advance(n);
		return *this;
	}
};
}  // end namespace detail

namespace detail {
template<class Element, ::boost::multi::dimensionality_type D, typename ElementPtr, bool IsConst, bool IsMove, typename Stride, class SubLayout>
struct array_iterator  // NOLINT(misc-multiple-inheritance) for facades
: boost::multi::iterator_facade<
	  array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride>, void, std::random_access_iterator_tag,
	  subarray<Element, D - 1, ElementPtr> const&, typename layout_t<D - 1>::difference_type>
, multi::decrementable<array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride>>
, multi::incrementable<array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride>>
, multi::affine<array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride>, multi::difference_type>
, multi::totally_ordered2<array_iterator<Element, D, ElementPtr, IsConst, IsMove, Stride>, void> {
	~array_iterator() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	constexpr auto operator=(array_iterator&&)  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		noexcept                                // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
		-> array_iterator& = default;

	array_iterator(array_iterator&&) noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
		= default;                             // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	using difference_type   = typename layout_t<D>::difference_type;
	using element           = Element;
	// using element_type      = Element;  // this creates a problem with std::ranges
	using element_ptr       = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<Element const>;
	using value_type        = typename subarray<Element, D - 1, element_ptr>::decay_type;

	using pointer   = void;  // subarray<element, D - 1, element_ptr>*;
	using reference = std::conditional_t<
		IsConst,
		const_subarray<element, D - 1, element_ptr>,
		subarray<element, D - 1, element_ptr>>;
	using const_reference = const_subarray<element, D - 1, element_ptr>;

	template<class Element2>
	using rebind = array_iterator<std::decay_t<Element2>, D, typename std::pointer_traits<ElementPtr>::template rebind<Element2>, IsConst, IsMove, Stride>;

	using iterator_category = std::random_access_iterator_tag;

	constexpr static dimensionality_type rank_v = D;

	using rank = std::integral_constant<dimensionality_type, D>;  // TODO(correaa) make rank a function for compat with mdspan?

	using ptr_type = subarray_ptr<element, D - 1, element_ptr, layout_t<D - 1>, true>;

	using stride_type = Stride;
	using layout_type = typename reference::layout_type;  // layout_t<D - 1>

	BOOST_MULTI_HD constexpr array_iterator() : ptr_{}, stride_{} {}  // = default;  // TODO(correaa) make = default, now it is not compiling

	template<class, dimensionality_type, class, bool, bool, typename, class> friend struct array_iterator;

	template<
		class EElement, typename PPtr, bool B, typename S, class L,
		decltype(multi::detail::explicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr, false, B, S, L>>().base()))* = nullptr>
	BOOST_MULTI_HD constexpr explicit array_iterator(array_iterator<EElement, D, PPtr, false, B, S, L> const& other)
	: ptr_{element_ptr{other.base()}, other.ptr_->layout()}, stride_{other.stride_} {}

	template<
		class EElement, typename PPtr, bool B, typename S, class L,
		decltype(multi::detail::implicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr, false, B, S, L>>().base()))* = nullptr>  // propagate implicitness of pointer
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr /*mplct*/ array_iterator(array_iterator<EElement, D, PPtr, false, B, S, L> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: ptr_(other.ptr_), stride_{other.stride_} {}

	array_iterator(array_iterator const&)                    = default;
	auto operator=(array_iterator const&) -> array_iterator& = default;

	BOOST_MULTI_HD constexpr explicit operator bool() const { return ptr_->base(); }  // TODO(correaa) implement bool conversion for subarray_ptr
	BOOST_MULTI_HD constexpr auto     operator*() const -> reference { return *ptr_; }

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wlarge-by-value-copy"  // TODO(correaa) can/should it be returned by reference?
#endif

	BOOST_MULTI_HD constexpr auto operator->() const -> decltype(auto) { return ptr_; }

#ifdef __clang__
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_HD constexpr auto operator+(difference_type n) const -> array_iterator {
		array_iterator ret{*this};
		ret += n;
		return ret;
	}
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> reference { return *((*this) + n); }

	template<bool OtherIsConst, std::enable_if_t<(IsConst != OtherIsConst), int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	BOOST_MULTI_HD constexpr auto operator==(array_iterator<Element, D, ElementPtr, OtherIsConst> const& other) const -> bool {
		// BOOST_MULTI_ASSERT( this->stride_ == other.stride_ );
		// BOOST_MULTI_ASSERT( this->ptr_->layout() == other.ptr_->layout() );
		// return (this->ptr_ == other.ptr_) && (this->stride_ == other.stride_) && (*(this->ptr_)).layout() == (*(other.ptr_)).layout();
		// compare the raw element pointers (`base()`) rather than the whole `subarray_ptr`s: comparing
		// `subarray_ptr == subarray_ptr` forms an overload set that, via ADL on the Thrust element pointer,
		// pulls in Thrust's `is_pointer_convertible_v`-constrained `operator==` and hard-errors evaluating
		// `thrust::detail::pointer_element<subarray_ptr>` (incomplete) under CCCL 3 (CUDA 13+).
		return this->base() == other.base() && this->stride_ == other.stride_ && this->ptr_->layout() == other.ptr_->layout();
	}

	template<bool OtherIsConst, std::enable_if_t<(IsConst != OtherIsConst), int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	BOOST_MULTI_HD constexpr auto operator!=(array_iterator<Element, D, ElementPtr, OtherIsConst> const& other) const -> bool {
		return !operator==(other);
	}

	BOOST_MULTI_HD constexpr auto operator==(array_iterator const& other) const -> bool {
		BOOST_MULTI_ASSERT(this->stride_ == other.stride_);
		BOOST_MULTI_ASSERT(this->ptr_->layout() == other.ptr_->layout());
		// compare raw element pointers (`base()`), not the whole `subarray_ptr`s, to avoid pulling in
		// Thrust's `pointer_element<subarray_ptr>` (incomplete) under CCCL 3 — see operator== above.
		return this->base() == other.base() && this->ptr_->layout() == other.ptr_->layout();
	}

	BOOST_MULTI_HD constexpr auto operator!=(array_iterator const& other) const -> bool {
		return !operator==(other);
	}

	BOOST_MULTI_HD constexpr auto operator<(array_iterator const& other) const -> bool {
		return 0 < other - *this;
	}

	BOOST_MULTI_HD constexpr explicit array_iterator(typename subarray<element, D - 1, element_ptr>::element_ptr base, layout_t<D - 1> const& lyt, stride_type stride)
	: ptr_(base, lyt), stride_{stride} {}

	template<class, dimensionality_type, class, class> friend class const_subarray;

	template<class... As>
	BOOST_MULTI_HD constexpr auto operator()(index idx, As... args) const -> decltype(auto) { return this->operator[](idx)(args...); }
	BOOST_MULTI_HD constexpr auto operator()(index idx) const -> decltype(auto) { return this->operator[](idx); }

 private:
	template<class Self, typename Tuple, std::size_t... I>
	static BOOST_MULTI_HD constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...> /*012*/) -> decltype(auto) {
		using std::get;  // for C++17 compatibility
		return std::forward<Self>(self)(get<I>(tuple)...);
	}

 public:
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl) const -> decltype(auto) { return apply_impl_(*this, tpl, std::make_index_sequence<std::tuple_size_v<Tuple>>()); }

 private:
	ptr_type    ptr_;
	stride_type stride_;

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	BOOST_MULTI_HD constexpr void decrement_() { ptr_.base_ -= stride_; }
	BOOST_MULTI_HD constexpr void advance_(difference_type n) { ptr_.base_ += stride_ * n; }  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	BOOST_MULTI_HD constexpr auto base() const -> element_ptr { return ptr_.base_; }
	BOOST_MULTI_HD constexpr auto stride() const -> stride_type { return stride_; }

	friend /*constexpr*/ auto base(array_iterator const& self) -> element_ptr { return self.base(); }
	friend constexpr auto     stride(array_iterator const& self) -> stride_type { return self.stride_; }  // TODO(correaa) remove

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	constexpr auto operator++() -> array_iterator& {
		ptr_.base_ += stride_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}

	constexpr auto operator--() -> array_iterator& {
		ptr_.base_ -= stride_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	friend constexpr auto operator-(array_iterator const& self, array_iterator const& other) -> difference_type {
		BOOST_MULTI_ASSERT(self.stride_ == other.stride_);  // LCOV_EXCL_LINE
		BOOST_MULTI_ASSERT(self.stride_ != 0);              // LCOV_EXCL_LINE
		return (self.ptr_.base() - other.ptr_.base()) / self.stride_;
	}

	constexpr auto operator+=(difference_type n) -> array_iterator& {
		advance_(+n);
		return *this;
	}

	constexpr auto operator-=(difference_type n) -> array_iterator& {
		advance_(-n);
		return *this;
	}
};
}  // end namespace detail

namespace detail {
/// A cursor is a lightweigh type for multidimensional indexing, it is similar to an array object without the size (extents) information
template<typename ElementPtr, dimensionality_type D, class StridesType>
struct cursor_t {
	/// Signed integer for index arithmetic in the leading dimension
	using difference_type = typename std::iterator_traits<ElementPtr>::difference_type;
	/// Tuple type to describe the strides of the array defined by the cursor
	using strides_type    = StridesType;

	/// Element pointer type (e.g. `T*`)
	using element_ptr                                = ElementPtr;
	/// Element reference type (e.g. `T&`)
	using element_ref                                = typename std::iterator_traits<element_ptr>::reference;
	/// Element reference type (e.g. `T`)
	using element                                    = typename std::iterator_traits<element_ptr>::value_type;
	using element_type [[deprecated("use element")]] = typename std::iterator_traits<element_ptr>::value_type;

	/// Pointer type (`void`, for compatibility)
	using pointer   = element_ptr;
	/// Reference type (`void`, for compatibility)
	using reference = element_ref;

	/// Tuple type of indices
	using indices_type = typename extents_t<D>::indices_type;

	cursor_t() = default;

 private:
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // '7' bytes padding added after data member 'boost::multi::array_types<T,2,ElementPtr,Layout>::base_' [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\test\array_fancyref.cpp.x.vcxproj]
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

	strides_type strides_;
	element_ptr  base_;

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif

	template<class, dimensionality_type, class, class> friend class multi::const_subarray;
	template<class, dimensionality_type, class> friend struct detail::cursor_t;

	BOOST_MULTI_HD constexpr cursor_t(element_ptr base, strides_type const& strides) : strides_{strides}, base_{base} {}

	template<class OtherCursor, class = decltype(multi::detail::implicit_cast<element_ptr>(std::declval<OtherCursor>().base()))>
	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr cursor_t(OtherCursor const& other) : strides_{other.strides()}, base_{other.base()} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	template<class OtherCursor>
	BOOST_MULTI_HD constexpr explicit cursor_t(OtherCursor const& other) : strides_{other.strides()}, base_{other.base()} {}

 public:
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	/// Indexing operator, it returns a cursors of lower dimensionality; recursively obtaining an element at the corresponding relative index
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> decltype(auto) {
		using std::get;  // for C++17 compatibility
		if constexpr(D != 1) {
			return cursor_t<ElementPtr, D - 1, std::decay_t<decltype(strides_.tail())>>(
				base_ + get<0>(strides_) * n,  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				strides_.tail()
			);
		} else {
			return base_[get<0>(strides_) * n];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		}
	}
#ifdef __clang__
#pragma clang diagnostic pop
#endif
	/// Function call operators, to obtain indexing, on one or more multiple indexing arguments
	BOOST_MULTI_HD constexpr auto operator()(difference_type n) const -> decltype(auto) {
		return operator[](n);
	}
	template<class... Ns>
	BOOST_MULTI_HD constexpr auto operator()(difference_type n, Ns... rest) const -> decltype(auto) {
		return operator[](n)(rest...);
	}

 private:
	template<class Tuple, std::size_t... I>
	BOOST_MULTI_HD constexpr auto apply_impl_(Tuple const& tup, std::index_sequence<I...> /*012*/) const -> decltype(auto) {
		using std::get;                                   // for C++17 compatibility
		return ((get<I>(tup) * get<I>(strides_)) + ...);  // NOLINT(readability-redundant-parentheses) for fold expression
	}

 public:
	template<class Tuple = indices_type>
	BOOST_MULTI_HD constexpr auto operator+=(Tuple const& tup) -> cursor_t& {
		base_ += apply_impl_(tup, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator*() const -> reference { return *base_; }
	BOOST_MULTI_HD constexpr auto operator->() const -> pointer { return base_; }

	BOOST_MULTI_HD constexpr auto base() const -> pointer { return base_; }
	BOOST_MULTI_HD constexpr auto strides() const -> strides_type { return strides_; }
	template<multi::dimensionality_type DD = 0>
	BOOST_MULTI_HD constexpr auto stride() const {
		using std::get;
		return get<DD>(strides_);
	}
};
}  // end namespace detail

namespace detail {
template<typename Pointer, class LayoutType>
struct elements_range_t;
}  // end namespace detail

namespace detail {
template<typename Pointer, class LayoutType>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
struct elements_iterator_t
// : boost::multi::random_accessable<elements_iterator_t<Pointer, LayoutType>, typename std::iterator_traits<Pointer>::difference_type, typename std::iterator_traits<Pointer>::reference>
{
	/// Signed integer type for iterator arithmetic
	using difference_type   = typename std::iterator_traits<Pointer>::difference_type;
	/// Value type obtained from iterator (e.g. `T`)
	using value_type        = typename std::iterator_traits<Pointer>::value_type;
	/// Pointer type of the element of the sequence (e.g. `T*`)
	using pointer           = Pointer;
	/// Reference type (e.g. `T&`)
	using reference         = std::remove_const_t<typename std::iterator_traits<Pointer>::reference>;  // TODO(correaa) investigate why top-level const reaches here
	/// Category of iterator (generally random access)
	using iterator_category = std::random_access_iterator_tag;
	/// Const pointer type of the element of the sequence (e.g. `T const*`)
	using const_pointer     = typename std::pointer_traits<pointer>::template rebind<value_type const>;

 private:
	/// Layout type of the iterator
	using layout_type = LayoutType;

	pointer                                base_;
	layout_type                            l_;
	difference_type                        n_ = 0;
	extents_t<layout_type::dimensionality> xs_;

	using indices_type = typename extents_t<layout_type::dimensionality>::indices_type;
	indices_type ns_   = {};

	template<typename, class> friend struct elements_iterator_t;
	template<typename, class> friend struct multi::detail::elements_range_t;

	BOOST_MULTI_HD constexpr elements_iterator_t(pointer base, layout_type const& lyt, difference_type n)
	: base_{std::move(base)}, l_{lyt}, n_{n}, xs_{l_.extents()}, ns_{lyt.is_empty() ? indices_type{} : xs_.from_linear(n)} {}

 public:
	elements_iterator_t() = default;

 private:
	/// Arithmetic base pointer of the iterator, typically the base of the original array
	BOOST_MULTI_HD constexpr auto base() -> pointer { return base_; }              // NOLINT(readability-identifier-naming) TODO(correaa) rename
	BOOST_MULTI_HD constexpr auto base() const -> const_pointer { return base_; }  // NOLINT(readability-identifier-naming) TODO(correaa) rename

	/// Layout used to calculate the position of an iterator, typically the layout of the original array
	BOOST_MULTI_HD constexpr auto layout() const -> layout_type { return l_; }  // NOLINT(readability-identifier-naming) TODO(correaa) rename

 public:
	template<class Other, decltype(multi::detail::implicit_cast<pointer>(std::declval<Other>().base_))* = nullptr>
	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr /*impl*/ elements_iterator_t(Other const& other) : elements_iterator_t{other.base_, other.l_, other.n_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	template<class Other>
	BOOST_MULTI_HD constexpr explicit elements_iterator_t(Other const& other) : elements_iterator_t{other.base_, other.l_, other.n_} {}

	elements_iterator_t(elements_iterator_t const&) = default;

	auto operator=(elements_iterator_t const&) -> elements_iterator_t& = default;

	BOOST_MULTI_HD constexpr auto operator++() -> elements_iterator_t& {
		apply([&exts = this->xs_](auto&... idxs) -> auto { return exts.next_canonical(idxs...); }, ns_);
		// std::apply([&xs = this->xs_](auto&... idxs) { return xs.next_canonical(idxs...); }, ns_);
		++n_;
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator--() -> elements_iterator_t& {
		std::apply([&exts = this->xs_](auto&... idxs) -> auto { return exts.prev_canonical(idxs...); }, ns_);
		--n_;
		return *this;
	}

	template<class = void>  // TODO(correaa) lazy instantiation to workaround MSVC linking limitations iterator copy for restrictions
	BOOST_MULTI_HD constexpr auto operator++(int) -> elements_iterator_t {
		elements_iterator_t ret{*this};
		++(*this);
		return ret;
	}

	template<class = void>  // TODO(correaa) lazy instantion to workaround MSVC linking limitations iterator copy for restrictions
	BOOST_MULTI_HD constexpr auto operator--(int) -> elements_iterator_t {
		elements_iterator_t ret{*this};
		--(*this);
		return ret;
	}

	BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> elements_iterator_t& {
		auto linear_n = apply(xs_, ns_);
		ns_           = xs_.from_linear(linear_n + n);
		n_ += n;
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator-=(difference_type n) -> elements_iterator_t& {
		// auto const nn = std::apply(xs_, ns_);
		// ns_ = xs_.from_linear(nn - n);
		n_ -= n;
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator-(elements_iterator_t const& other) const noexcept -> difference_type {
		BOOST_MULTI_ASSERT(base_ == other.base_ && l_ == other.l_);
		return n_ - other.n_;
	}

	BOOST_MULTI_HD constexpr auto operator<(elements_iterator_t const& other) const noexcept -> bool {
		BOOST_MULTI_ASSERT(base_ == other.base_ && l_ == other.l_);
		return n_ < other.n_;
	}

	BOOST_MULTI_HD constexpr auto operator<=(elements_iterator_t const& other) const noexcept -> bool { return ((*this) < other) || ((*this) == other); }

	BOOST_MULTI_HD constexpr auto operator>(elements_iterator_t const& other) const noexcept -> bool { return other < (*this); }
	BOOST_MULTI_HD constexpr auto operator>=(elements_iterator_t const& other) const noexcept -> bool { return !((*this) < other); }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	BOOST_MULTI_HD constexpr auto current() const -> pointer { return base_ + std::apply(l_, ns_); }

	/// Dereference operator, gets a reference to the element pointed by this iterator
	BOOST_MULTI_HD constexpr auto operator*() const -> reference {  // cppcheck-suppress duplInheritedMember ; to overwrite
		return base_[apply(l_, ns_)];                               // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}

	/// Subscript operator, refers to the `n`th element after or before this iterator position of positive or negative `n` respectively (same as `*(it + n)`)
	BOOST_MULTI_HD constexpr auto operator[](difference_type const& n) const -> reference {
		auto const linear_n = apply(xs_, ns_);
		return base_[apply(l_, xs_.from_linear(linear_n + n))];
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_HD constexpr auto operator+(difference_type n) const -> elements_iterator_t {
		auto ret{*this};
		ret += n;
		return ret;
	}

	friend BOOST_MULTI_HD constexpr auto operator+(difference_type n, elements_iterator_t const& self) -> elements_iterator_t {  // `n + it` form, required by std::random_access_iterator
		return self + n;
	}

	BOOST_MULTI_HD constexpr auto operator-(difference_type n) const -> elements_iterator_t {
		auto ret{*this};
		ret -= n;
		return ret;
	}

	BOOST_MULTI_HD constexpr auto operator==(elements_iterator_t const& other) const -> bool {
		BOOST_MULTI_ASSERT(base_ == other.base_ && l_ == other.l_);  // TODO(correaa) calling host function from host device
		return n_ == other.n_;                                       // and base_ == other.base_ and l_ == other.l_;
	}
	BOOST_MULTI_HD constexpr auto operator!=(elements_iterator_t const& other) const -> bool {
		BOOST_MULTI_ASSERT(base_ == other.base_ && l_ == other.l_);  // TODO(correaa) calling host function from host device
		return n_ != other.n_;
	}

	~elements_iterator_t() = default;
};
}  // end namespace detail

namespace detail {
template<typename Pointer, class LayoutType>
struct elements_range_t {
	/// Pointer type to element (e.g. `T*`)
	using pointer = Pointer;

	/// Element type (e.g. `T`)
	using value_type    = typename std::iterator_traits<pointer>::value_type;
	/// Const-pointer type to element (e.g. `T const*`)
	using const_pointer = typename std::pointer_traits<pointer>::template rebind<value_type const>;

	/// Reference type to element (e.g. `T&`)
	using reference       = typename std::iterator_traits<pointer>::reference;
	/// Const-eference type to element (e.g. `T const&`)
	using const_reference = typename std::iterator_traits<const_pointer>::reference;

	/// Integer type that holds the size of the range
	using size_type       = typename std::iterator_traits<pointer>::difference_type;
	/// Signed integer type for index arithmetic
	using difference_type = typename std::iterator_traits<pointer>::difference_type;

	/// Random access iterator type to access each element
	using iterator       = detail::elements_iterator_t<pointer, LayoutType>;
	/// Random access const-iterator type to access each element
	using const_iterator = detail::elements_iterator_t<const_pointer, LayoutType>;

 private:
	/// Element type (e.g. `T`)
	using element                                      = value_type;
	using element_type [[deprecated("use ::element")]] = value_type;

	/// Layout of the range (e.g. `T*`)
	using layout_type = LayoutType;

	pointer     base_;
	layout_type l_;

 public:
	template<class OtherRange, decltype(multi::detail::implicit_cast<pointer>(std::declval<OtherRange>().base_))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*impl*/ elements_range_t(OtherRange const& other) : base_{other.base}, l_{other.l_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) to reproduce the implicitness of the argument
	template<class OtherRange, decltype(multi::detail::explicit_cast<pointer>(std::declval<OtherRange>().base_))* = nullptr>
	constexpr explicit elements_range_t(OtherRange const& other) : elements_range_t{other} {}

	constexpr elements_range_t(pointer base, layout_type const& lyt) : base_{std::move(base)}, l_{lyt} {}

	constexpr auto base() -> pointer { return base_; }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0
	constexpr auto base() const -> const_pointer { return base_; }

	constexpr auto layout() const -> layout_type { return l_; }

 private:
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	constexpr auto at_aux_(difference_type n) const -> reference {
		BOOST_MULTI_ASSERT(!is_empty());
		return base_[std::apply(l_, l_.extents().from_linear(n))];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const& -> const_reference { return at_aux_(n); }
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) && -> reference { return at_aux_(n); }
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) & -> reference { return at_aux_(n); }

	constexpr auto size() const noexcept -> size_type { return l_.num_elements(); }

	using extent_type                                      = multi::extent_t<index>;
	using extension_type [[deprecated("use extent_type")]] = extent_type;  // NOLINT

	[[deprecated("use extent")]] BOOST_MULTI_HD constexpr auto extension() const { return extent_type{0, size()}; }
	[[nodiscard]] BOOST_MULTI_HD constexpr auto                extent() const { return extent_type{0, size()}; }

	[[nodiscard]]
	constexpr auto empty() const -> bool { return l_.empty(); }
	/// Checks if the container has no elements
	constexpr auto is_empty() const -> bool { return l_.is_empty(); }

	elements_range_t(elements_range_t const&) = delete;
	elements_range_t(elements_range_t&&)      = default;

	template<class Range> auto operator==(Range const& other) const -> bool {
		return size() == other.size() && adl_equal(other.begin(), other.end(), begin());
	}

	template<class Range> auto operator!=(Range const& other) const -> bool {
		return size() != other.size() || !adl_equal(other.begin(), other.end(), begin());
	}

	/// swaps with `other` range, element by element (it doesn't rebind, O(N) operation)
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>& other) & noexcept {
		BOOST_MULTI_ASSERT(size() == other.size());
		adl_swap_ranges(begin(), end(), other.begin());
	}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>& other) && noexcept {
		BOOST_MULTI_ASSERT(size() == other.size());
		adl_swap_ranges(begin(), end(), other.begin());
	}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&& other) & noexcept {
		BOOST_MULTI_ASSERT(size() == other.size());
		adl_swap_ranges(begin(), end(), std::move(other).begin());
	}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&& other) && noexcept {
		BOOST_MULTI_ASSERT(size() == other.size());
		adl_swap_ranges(begin(), end(), std::move(other).begin());
	}

	~elements_range_t() = default;

 private:
	BOOST_MULTI_HD constexpr auto begin_aux_() const { return iterator{base_, l_, 0}; }
	BOOST_MULTI_HD constexpr auto end_aux_() const { return iterator{base_, l_, l_.num_elements()}; }

 public:
	/// returns an iterator to the beginning of the range
	BOOST_MULTI_HD constexpr auto begin() const& -> const_iterator { return begin_aux_(); }
	/// returns an iterator to the end of the range
	BOOST_MULTI_HD constexpr auto end() const& -> const_iterator { return end_aux_(); }

	BOOST_MULTI_HD constexpr auto begin() && -> iterator { return begin_aux_(); }
	BOOST_MULTI_HD constexpr auto end() && -> iterator { return end_aux_(); }

	BOOST_MULTI_HD constexpr auto begin() & -> iterator { return begin_aux_(); }
	BOOST_MULTI_HD constexpr auto end() & -> iterator { return end_aux_(); }

	/// yields the front element of a non-empty range.
	/// \pre `!is_empty()`
	BOOST_MULTI_HD constexpr auto front() const& -> const_reference { return *begin(); }

	/// yields the back element of a non-empty range.
	/// \pre `!is_empty()`
	BOOST_MULTI_HD constexpr auto back() const& -> const_reference { return *std::prev(end(), 1); }

	BOOST_MULTI_HD constexpr auto front() && -> reference { return *begin(); }
	BOOST_MULTI_HD constexpr auto back() && -> reference { return *std::prev(end(), 1); }

	BOOST_MULTI_HD constexpr auto front() & -> reference { return *begin(); }
	BOOST_MULTI_HD constexpr auto back() & -> reference { return *std::prev(end(), 1); }

	/// Assignment operator, copies each element, sizes must match
	auto operator=(elements_range_t const& other) -> elements_range_t& {
		if(this == std::addressof(other)) {
			return *this;
		}
		assert(other.size() == this->size());
		if(!is_empty()) {
			adl_copy(other.begin(), other.end(), this->begin());
		}
		return *this;
	}

	/// Assignment operator, copies each element, sizes must match (O(N) operation)
	constexpr auto operator=(elements_range_t&& other) noexcept(false) -> elements_range_t& {  // cannot be =delete in NVCC?  // NOLINT(bugprone-unsafe-to-allow-exceptions)
		if(!is_empty()) {
			adl_copy(other.begin(), other.end(), this->begin());
		}
		(void)std::move(other);
		return *this;
	}

	template<class OtherElementRange, class = decltype(adl_copy(std::begin(std::declval<OtherElementRange&&>()), std::end(std::declval<OtherElementRange&&>()), std::declval<iterator>()))>
	auto operator=(OtherElementRange&& other) & -> elements_range_t& {  // NOLINT(cppcoreguidelines-missing-std-forward) std::forward<OtherElementRange>(other) creates a problem with move-only elements
		BOOST_MULTI_ASSERT(size() == other.size());
		if(!is_empty()) {
			adl_copy(std::begin(other), std::end(other), begin());
		}
		return *this;
	}

	template<class OtherElementRange, class = decltype(adl_copy(std::begin(std::declval<OtherElementRange&&>()), std::end(std::declval<OtherElementRange&&>()), std::declval<iterator>()))>
	constexpr auto operator=(OtherElementRange&& other) && -> elements_range_t& {  // NOLINT(cppcoreguidelines-missing-std-forward) std::forward<OtherElementRange>(other) creates a problem with move-only elements
		operator=(std::forward<OtherElementRange>(other));
		// BOOST_MULTI_ASSERT(size() == other.size());
		// if(!is_empty()) {
		// 	adl_copy(std::begin(other), std::end(other), begin());
		// }
		return *this;
	}

	auto operator=(std::initializer_list<value_type> values) && -> elements_range_t& {
		operator=(values);
		return *this;
	}
	auto operator=(std::initializer_list<value_type> values) & -> elements_range_t& {
		BOOST_MULTI_ASSERT(static_cast<size_type>(values.size()) == size());
		adl_copy_n(values.begin(), values.size(), begin());
		return *this;
	}
};
}  // end namespace detail

template<class It>
[[deprecated("remove")]] BOOST_MULTI_HD constexpr auto ref(It begin, It end)
	-> multi::subarray<typename It::element, It::rank_v, typename It::element_ptr> {
	return multi::subarray<typename It::element, It::rank_v, typename It::element_ptr>{begin, end};
}

/// A `D`-dimensional array whose size is bound at construction and never changes.
///
/// Provides contiguous, allocator-managed storage with no reallocation after construction.
/// Pointers, references, and iterators to elements remain valid for the lifetime of the object.
///
/// @note Assignments require matching extensions; use `array` when resizing or full value semantics is needed
///
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam Alloc Allocator type
template<typename T, ::boost::multi::dimensionality_type D, class Alloc> struct dynamic_array;  // this might be needed by MSVC 14.3 in c++17 mode

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<typename T, ::boost::multi::dimensionality_type D, typename ElementPtr, class Layout>
class const_subarray : public detail::array_types<T, D, ElementPtr, Layout> {
	using types = detail::array_types<T, D, ElementPtr, Layout>;  // TODO(correaa) eliminate

 public:
	using detail::array_types<T, D, ElementPtr, Layout>::rank_v;

	friend class const_subarray<typename types::element, D + 1, typename types::element_ptr>;

	// using typename types::element_type;
	using types::layout;

	/// type that holds the layout of the subarray or array
	using layout_type = Layout;

	/// A tuple type to hold all the indices necessary to get to the element of the array
	// using indices_type = typename extents_t<D>::indices_type;

	/// A type to hold the size of an array or subarray (size in the leading dimension) (usually signed)
	using size_type = typename detail::array_types<T, D, ElementPtr, Layout>::size_type;

	/// returns the internal layout information of the array
	BOOST_MULTI_HD constexpr auto layout() const -> layout_type { return detail::array_types<T, D, ElementPtr, Layout>::layout(); }  // cppcheck-suppress duplInheritedMember

	/// Assingment operators (disabled because the const-subarray is immutable)
	auto operator=(const_subarray const&) -> const_subarray& = delete;
	auto operator=(const_subarray&&) -> const_subarray&      = delete;

	const_subarray() = default;
	BOOST_MULTI_HD constexpr const_subarray(layout_type const& layout, ElementPtr const& base)
	: detail::array_types<T, D, ElementPtr, Layout>{layout, base} {}

 protected:
	// using types::types;
	BOOST_MULTI_HD constexpr explicit const_subarray(std::nullptr_t nil) : types{nil} {}

	template<typename, ::boost::multi::dimensionality_type, class Alloc> friend struct dynamic_array;

	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;

	// TODO(correaa) vvv consider making it explicit (seems that in C++23 it can prevent auto s = a[0];)
	// const_subarray(const_subarray const&) = default;  // NOTE: reference type cannot be copied. perhaps you want to return by std::move or std::forward if you got the object from a universal reference argument

 public:
	const_subarray(const_subarray const&) = delete;

	/// Element type (`T`)
	using element                                      = T;
	/// (deprecated, use `element`)
	using element_type [[deprecated("use ::element")]] = typename types::element;
	/// Pointer to element type (usually `T*`)
	using element_ptr                                  = typename types::element_ptr;
	/// Pointer to immutable element type (usually `T const*`)
	using element_const_ptr                            = typename types::element_const_ptr;
	/// Reference type to an element (usually `T&`)
	using element_ref                                  = typename types::element_ref;
	/// `const`-qualified reference type to an element (usually `T const&`)
	using element_cref                                 = typename std::iterator_traits<element_const_ptr>::reference;

 private:
	using index_gen [[deprecated("here to fulfill backward-compatible MultiArray concept")]]    = char*;
	using extent_gen [[deprecated("here to fulfill backward-compatible  MultiArray concept")]]  = void;
	using extent_range [[deprecated("here to fulfill backward-compatible MultiArray concept")]] = void;

	// private:
	using elements_iterator  = detail::elements_iterator_t<element_ptr, layout_type>;
	using celements_iterator = detail::elements_iterator_t<element_const_ptr, layout_type>;

	using const_elements_range = detail::elements_range_t<element_const_ptr, layout_type>;
	using elements_range       = detail::elements_range_t<element_ptr, layout_type>;

	constexpr auto elements_aux_() const { return elements_range(this->base_, this->layout()); }

	template<class TT> using il_ = std::initializer_list<TT>;

 public:
	const_subarray(const_subarray&&) noexcept = default;  // lints(readability-redundant-access-specifiers)

	// NOLINTBEGIN(google-explicit-constructor,hicpp-explicit-conversions,modernize-use-constraints,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) for C++20
	template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	explicit const_subarray(il_<TT> const& il_1d) : const_subarray(multi::layout(il_1d), multi::base(il_1d)) {}
	// template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	// explicit const_subarray(il_<il_<TT>> const& il_2d) : const_subarray(multi::layout(il_2d), multi::base(il_2d)) {}
	// template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	// explicit const_subarray(il_<il_<il_<TT>>> const& il_3d) : const_subarray(multi::layout(il_3d), multi::base(il_3d)) {}
	// NOLINTEND(google-explicit-constructor,hicpp-explicit-conversions,modernize-use-constraints,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) for C++20

#if defined(__cpp_deleted_function) && (__cpp_deleted_function >= 202403L) && (__cplusplus > 202302L)
#define BOOST_MULTI_DELETE(ReasoN) delete (ReasoN)
#else
#define BOOST_MULTI_DELETE(ReasoN) delete
#endif

	// NOLINTBEGIN(google-explicit-constructor,modernize-use-constraints)
	template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	explicit const_subarray(il_<TT>&& /*il_1d*/) = BOOST_MULTI_DELETE("temporary init-list dangles");
	template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	explicit const_subarray(il_<il_<TT>>&& /*il_2d*/) = BOOST_MULTI_DELETE("temporary init-list dangles");
	template<class TT = T, std::enable_if_t<std::is_convertible_v<decltype(multi::base(std::declval<il_<TT> const&>())), ElementPtr>, int> = 0>
	explicit const_subarray(il_<il_<il_<TT>>>&& /*il_3d*/) = BOOST_MULTI_DELETE("temporary init-list dangles");
	// NOLINTEND(google-explicit-constructor,modernize-use-constraints)

#undef BM_DELETE

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
 private:
	auto to_mdspan_aux_() const {
		using std::apply;
		auto shape = apply(
			[](auto... sizes) -> auto { return std::dextents<std::size_t, D>{static_cast<std::size_t>(sizes)...}; },
			this->sizes()
		);

		auto strides = apply(
			[](auto... strds) -> auto { return std::array<std::size_t, D>{static_cast<std::size_t>(strds)...}; },
			this->strides()
		);

		return std::mdspan<T, std::dextents<std::size_t, D>, std::layout_stride>{
			this->base_, std::layout_stride::mapping(shape, strides)
		};
	}

 public:
	operator std::mdspan<T const, std::dextents<std::size_t, D>, std::layout_stride>() const& { return to_mdspan_aux_(); }
#endif

	/// possibly moves the contents
	friend BOOST_MULTI_HD constexpr auto move(const_subarray const& self) -> const_subarray { return const_subarray(self.layout(), self.base_); }

	/// returns a random-access range with all the elements of the array
	constexpr auto elements() const& { return const_elements_range(this->base(), this->layout()); }  // cppcheck-suppress duplInheritedMember ; to overwrite

 private:
	constexpr auto const_elements_() const -> const_elements_range { return elements_aux_(); }

 public:
	~const_subarray() = default;  // this lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	BOOST_MULTI_FRIEND_CONSTEXPR auto sizes(const_subarray const& self) noexcept -> typename const_subarray::sizes_type { return self.sizes(); }  // needed by nvcc
#if (defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ < 12 || (__CUDACC_VER_MAJOR__ == 12 && __CUDACC_VER_MINOR__ < 4))) || \
	(defined(__NVCOMPILER) /*&& (__NVCOMPILER_MAJOR__ < 23 || (__NVCOMPILER_MAJOR__ == 23 && __NVCOMPILER_MINOR__ < 11))*/)
	template<class = void> friend constexpr auto size(const_subarray const& self) noexcept -> typename const_subarray::size_type { return self.size(); }
#endif

	friend constexpr auto dimensionality(const_subarray const& /*self*/) { return D; }

 private:
	using default_allocator_type = typename multi::pointer_traits<const_subarray::element_ptr>::default_allocator_type;

 public:
	/// returns the associated allocator
	BOOST_MULTI_HD constexpr auto get_allocator() const -> default_allocator_type {
		using multi::get_allocator;
		return get_allocator(this->base());
	}

	// BOOST_MULTI_FRIEND_CONSTEXPR auto get_allocator(const_subarray const& self) -> default_allocator_type { return self.get_allocator(); }

	/// Associated type that this reference decays to (when true copies are needed)
	using decay_type = array<std::decay_t<typename types::element>, D, typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type>;

	/// materializes an independent, owning `array` copy of this view with the associated array-value type (use unary prefix `+` as a shortcut)
	constexpr auto decay() const& -> decay_type {  // cppcheck-suppress duplInheritedMember ; to overwrite
		decay_type ret{*this};
		return ret;
	}

	/// materializes an independent, owning `array` copy of this view with the associated array-value type (use unary prefix `+` as a shortcut)
	[[deprecated]] constexpr auto copy() const& -> decay_type {  // cppcheck-suppress duplInheritedMember ; to overwrite
		decay_type ret{*this};
		return ret;
	}

	/// materializes an independent, owning `array` copy of this view with the associated array-value type
	constexpr auto operator+() const -> decay_type { return decay(); }

	/// reference to a subarray of lower dimension (or to an element for `D == 1`)
	using reference = typename std::conditional_t<
		(D > 1),
		const_subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_ptr>::reference>;

	/// const-reference to a subarray of lower dimension (or to an element for `D == 1`)
	using const_reference = typename std::conditional_t<
		(D > 1),
		const_subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_const_ptr>::reference>;

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	BOOST_MULTI_HD constexpr auto at_aux_(index idx) const {
		BOOST_MULTI_ASSERT((/*this->stride() == 0 ||*/ (this->extent().contains(idx))) && ("out of bounds"));

		// clang-format off
	#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif
		return reference(
			this->layout().sub(),
			this->base_ + (/*idx **/ this->layout().stride() * idx - this->layout().offset())  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		);  // cppcheck-suppress syntaxError ; bug in cppcheck 2.5
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif
	}

 public:
 
	// clang-format off
	#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif
	// clang-format on

	BOOST_MULTI_HD constexpr auto operator[](index idx) const& -> const_reference {  // cppcheck-suppress duplInheritedMember ; to overwrite
		BOOST_MULTI_ASSERT((this->stride() == 0 || (this->extent().contains(idx))) && ("out of bounds"));
		return const_reference(
			this->layout().sub(),
			this->base_ + (idx * this->layout().stride() - this->layout().offset())  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		);                                                                           // cppcheck-suppress syntaxError ; bug in cppcheck 2.5
	}

	// clang-format off
	#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
	#pragma clang diagnostic pop
	#endif
	// clang-format on

	/// yields the front subarray of leading dimension (or element for `D == 1`) of a non-empty array.
	/// \pre `!is_empty()`
	constexpr auto front() const& -> const_reference {
		assert(!this->is_empty());
		return *begin();
	}
	/// yields the back subarray of leading dimension (or element for `D == 1`) of a non-empty array.
	/// \pre `!is_empty()`
	constexpr auto back() const& -> const_reference {
		assert(!this->is_empty());
		return *(end() - 1);
	}

	using typename types::index;

	/// yields an equivalent subarray with a specific starting index
	constexpr auto reindexed(index first) const& {  // NOLINT(readability-identifier-naming) TODO(correaa) decide if making it public
		return const_subarray(this->layout().reindex(first), types::base_);
	}
	constexpr auto reindexed(index first) & {  // NOLINT(readability-identifier-naming) TODO(correaa) decide if making it public
		return const_subarray(this->layout().reindex(first), types::base_);
	}
	constexpr auto reindexed(index first) && { return const_subarray(this->layout().reindex(first), types::base_); }  // NOLINT(readability-identifier-naming) TODO(correaa) decide if making it public

	// TODO(correaa) : implement reindexed_aux
	template<class... Indexes>
	constexpr auto reindexed(index first, Indexes... idxs) const& -> const_subarray {  // NOLINT(readability-identifier-naming) TODO(correaa) decide if making it public
		return reindexed(first).rotated().reindexed(idxs...).unrotated();
	}

 private:
	constexpr auto taked_aux_(difference_type n) const {
		BOOST_MULTI_ASSERT(n <= this->size());
		return const_subarray(this->layout().take(n), this->base_);
	}

 public:
	/// yields an array-view of the same dimensionality taking the first count subarrays in the leading dimension
	constexpr auto taked(difference_type count) const& { return taked_aux_(count).as_const(); }  // cppcheck-suppress duplInheritedMember;

 private:
	BOOST_MULTI_HD constexpr auto halved_aux_() const {
		auto new_layout = this->layout().halve();
		return subarray<T, D + 1, element_ptr>(new_layout, this->base_);
	}

	// public:
	BOOST_MULTI_HD constexpr auto halved() const& -> const_subarray<T, D + 1, element_ptr> { return halved_aux_(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove

// private:
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	constexpr auto dropped_aux_(difference_type n) const {
		BOOST_MULTI_ASSERT(n <= this->size());
		typename types::layout_type const new_layout(
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride() * (this->size() - n)
		);

		return subarray<T, D, ElementPtr, typename types::layout_type>(new_layout, this->base_ + n * this->layout().stride());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	/// yields a subarray of the same dimensionally in which the first `count` elements are dropped.
	constexpr auto dropped(difference_type count) const& { return dropped_aux_(count).as_const(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

 private:
	BOOST_MULTI_HD constexpr auto sliced_aux_(index first, index last) const {
		// TODO(correaa) remove first == last condition
		BOOST_MULTI_ASSERT(((first == last) || this->extent().contains(first)) && ("sliced first out of bounds"));
		BOOST_MULTI_ASSERT(((first == last) || this->extent().contains(last - 1)) && ("sliced last  out of bounds"));
		typename types::layout_type new_layout = this->layout();
		new_layout.nelems()                    = this->stride() * (last - first);                               // TODO(correaa) : reconstruct layout instead of mutating it
		BOOST_MULTI_ASSERT(this->base_ || ((first * this->layout().stride() - this->layout().offset()) == 0));  // it is UB to offset a nullptr

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

		return const_subarray(new_layout, this->base_ + (first * this->layout().stride() - this->layout().offset()));  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

 public:
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) const& -> const_subarray { return sliced_aux_(first, last); }

 private:
	constexpr auto blocked(index first, index last) const& { return sliced(first, last).reindexed(first).as_const(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or make private
	constexpr auto blocked(index first, index last) & { return sliced(first, last).reindexed(first).as_const(); }       // NOLINT(readability-identifier-naming) TODO(correaa) remove or make private

	using iextension = typename const_subarray::index_extension;

	// public:
	/// yields a stencil of the original array with the same dimensionality
	constexpr auto stenciled(iextension iex) & { return blocked(iex.first(), iex.last()).as_const(); }                                                                                            // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1) & { return stenciled(iex).rotated().stenciled(iex1).unrotated().as_const(); }                                                       // TODO(correaa) fix const  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2) & -> const_subarray { return stenciled(iex).rotated().stenciled(iex1, iex2).unrotated(); }                         // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3) & -> const_subarray { return stenciled(iex).rotated().stenciled(iex1, iex2, iex3).unrotated(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs) & -> const_subarray { return stenciled(iex).rotated().stenciled(iex1, iex2, iex3, iexs...).unrotated(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename

	constexpr auto stenciled(iextension iex) && { return blocked(iex.first(), iex.last()).as_const(); }                                                                                            // TODO(correaa) fix const  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1) && { return stenciled(iex).rotated().stenciled(iex1).unrotated().as_const(); }                                                       // TODO(correaa) fix const  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2) && { return stenciled(iex).rotated().stenciled(iex1, iex2).unrotated().as_const(); }                                // TODO(correaa) fix const  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3) && -> const_subarray { return stenciled(iex).rotated().stenciled(iex1, iex2, iex3).unrotated(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs) && -> const_subarray { return stenciled(iex).rotated().stenciled(iex1, iex2, iex3, iexs...).unrotated(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename

	constexpr auto stenciled(iextension iex) const& { return blocked(iex.first(), iex.last()).as_const(); }                                                                                     // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1) const& { return stenciled(iex).rotated().stenciled(iex1).unrotated().as_const(); }                                                // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2) const& { return stenciled(iex).rotated().stenciled(iex1, iex2).unrotated().as_const(); }                         // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3) const& { return stenciled(iex).rotated().stenciled(iex1, iex2, iex3).unrotated().as_const(); }  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename

	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs) const& {  // NOLINT(readability-identifier-naming) TODO(correaa) remove or rename
		return stenciled(iex).rotated().stenciled(iex1, iex2, iex3, iexs...).unrotated().as_const();
	}

	constexpr auto strided_aux_(difference_type step) const {
		assert(this->size() % step == 0);
		typename types::layout_type const new_layout{this->layout().sub(), this->layout().stride() * step, this->layout().offset(), this->layout().nelems()};
		return subarray<T, D, ElementPtr, typename types::layout_type>(new_layout, types::base_);
	}

 public:
	/// A subarray-view of the array with skipping `step` in the leading dimension
	constexpr auto strided(difference_type step) const& { return strided_aux_(step).as_const(); }

	/// A subarray-view from index `first` to index `last` (not inclusive) skipping `step` in the leading dimension
	constexpr auto sliced(
		typename types::index first, typename types::index last, typename types::index step
	) const& -> const_subarray {
		return sliced(first, last).strided(step);
	}

 private:
	// NOLINTNEXTLINE(readability-identifier-naming) TODO(correaa) to remove
	BOOST_MULTI_HD constexpr auto range(index_range irng) const& -> decltype(auto) { return sliced(irng.front(), irng.front() + irng.size()); }  // cppcheck-suppress duplInheritedMember;

 public:
	[[deprecated("is_flattable will be a property of the layout soon")]]
	constexpr auto is_flattable() const -> bool {
		return layout().is_flattable();
	}

 private:
	auto flattened_aux_() const {
		auto new_layout = this->layout().flatten(this->base_);
		return multi::subarray<T, D - 1, ElementPtr, decltype(new_layout)>(new_layout, this->base_);
	}

 public:
	auto flattened() const& { return flattened_aux_().as_const(); }

	constexpr auto flatted() const& {
		assert(this->layout().is_flattable());
		multi::layout_t<D - 1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return const_subarray<T, D - 1, ElementPtr>{new_layout, this->base_};
	}

 private:
	constexpr auto broadcasted() const& {  // NOLINT(readability-identifier-naming)  TODO(correaa) remove?
		// TODO(correaa) introduce a broadcasted_layout?
		multi::layout_t<D + 1> const new_layout(layout(), 0, 0);  //, (std::numeric_limits<size_type>::max)());  // paren for MSVC macros
		return const_subarray<T, D + 1, typename const_subarray::element_const_ptr>{new_layout, types::base_};
	}

	constexpr auto diagonal_aux_() const {
		using boost::multi::detail::get;
		auto                   square_size = (std::min)(get<0>(this->sizes()), get<1>(this->sizes()));  // paren for MSVC macros
		multi::layout_t<D - 1> new_layout{(*this)({0, square_size}, {0, square_size}).layout().sub()};
		new_layout.nelems() += (*this)({0, square_size}, {0, square_size}).layout().nelems();  // TODO(correaa) : don't use mutation
		new_layout.stride() += (*this)({0, square_size}, {0, square_size}).layout().stride();  // TODO(correaa) : don't use mutation
		return subarray<T, D - 1, typename const_subarray::element_ptr>(new_layout, types::base_);
	}

 public:
	/// Subarray of lower dimension that represents the main diagonal of a square array
	template<class Dummy = void, std::enable_if_t<(D > 1) && sizeof(Dummy*), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	constexpr auto diagonal() const& { return this->diagonal_aux_().as_const(); }

 private:
	BOOST_MULTI_HD constexpr auto partitioned_aux_(multi::ssize_t n) const {
		BOOST_MULTI_ASSERT(n != 0);
		// vvv TODO(correaa) should be size() here?
		BOOST_MULTI_ASSERT((this->layout().nelems() % n) == 0);  // if you get an assertion here it means that you are partitioning an array with an incommunsurate partition
		multi::layout_t<D + 1> new_layout{this->layout(), this->layout().nelems() / n, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= n;
		return subarray<T, D + 1, element_ptr>(new_layout, types::base_);
	}

 public:
	/// yields a subarray of higher dimension by splitting the leading dimension into `n` equal-sized partitions (`n` must divide `size()`)
	BOOST_MULTI_HD constexpr auto partitioned(size_type n) const& -> const_subarray<T, D + 1, element_ptr> { return partitioned_aux_(n); }

 private:
	BOOST_MULTI_HD constexpr auto chunked_aux_(size_type count) const {
		BOOST_MULTI_ASSERT(this->size() % count == 0);
		return partitioned_aux_(this->size() / count);
	}

 public:  // in Mathematica this is called Partition https://reference.wolfram.com/language/ref/Partition.html in RangesV3 it is called chunk
	/// produces a subarray of higher dimension by chunking in the leading dimension (if `count` doesn't divide `size`, so elements are left out at the end)
	BOOST_MULTI_HD constexpr auto chunked(size_type count) const& -> const_subarray<T, D + 1, element_ptr> { return chunked_aux_(count); }

	// constexpr auto tiled(size_type count) const& {
	// 	BOOST_MULTI_ASSERT(count != 0);
	// 	struct divided_type {
	// 		const_subarray<T, D + 1, element_ptr> quotient;
	// 		const_subarray<T, D, element_ptr>     remainder;
	// 	};
	// 	return divided_type(
	// 		this->taked(this->size() - (this->size() % count)).chunked(count),
	// 		this->dropped(this->size() - (this->size() % count))
	// 	);
	// }

 private:
	constexpr auto reversed_aux_() const { return const_subarray(layout().reverse(), types::base_); }

 public:
	/// yields a view of the array with the leading dimension in the reverse order (e.g. `a.reversed()[i][j] == a[a.size() - 1 - i][j]`)
	constexpr auto reversed() const& { return reversed_aux_().as_const(); }
	constexpr auto reversed() & -> const_subarray { return reversed_aux_(); }
	constexpr auto reversed() && -> const_subarray { return reversed_aux_(); }

 private:
	BOOST_MULTI_HD constexpr auto transposed_aux_() const {
		return const_subarray(layout().transpose(), types::base_);  // TODO(correaa) should be subarray(...)?
	}

 public:
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"  // TODO(correaa) for latex documentation in MrDocs
#endif
	/// A transpose view $A^\mathrm{T}$, that exchanges the first two indices
	BOOST_MULTI_HD constexpr auto transposed() const& -> const_subarray { return transposed_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD auto operator~(const_subarray const& self) -> const_subarray { return self.transposed(); }

 private:
	BOOST_MULTI_HD constexpr auto rotated_aux_() const {
		return const_subarray(layout().rotate(), types::base_);
	}
	BOOST_MULTI_HD constexpr auto unrotated_aux_() const {
		return const_subarray(layout().unrotate(), types::base_);
	}

 public:
	/// yields a view of the array where the indices are rotated (to the left) (e.g. `a.rotated()[i][j][k] == a[j][k][i]`)
	BOOST_MULTI_HD constexpr auto rotated() const& -> const_subarray { return rotated_aux_(); }
	/// yields a view of the array where the indices are unrotated (to the right, opposite to `.rotated()`) (e.g. `a.unrotated()[i][j][k]== a[k][i][j]`, `a.rotated().unrotated()` is the same as `a`)
	BOOST_MULTI_HD constexpr auto unrotated() const& -> const_subarray { return unrotated_aux_(); }

 private:
	BOOST_MULTI_HD constexpr auto unordered_aux_() const {
		return const_subarray(layout().sort(), types::base_);
	}

 public:
	/// yields a view in which index access is unordered (an arbitrary transposition of indices generally to optimize access)
	BOOST_MULTI_HD constexpr auto unordered() const {
		return unordered_aux_().as_const();
	}

 private:
	template<typename, ::boost::multi::dimensionality_type, typename, class> friend class const_subarray;

	BOOST_MULTI_HD constexpr auto paren_aux_() const& { return const_subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_); }

 public:
	/// Subarray returning operator (takes multiple parameters, the number of parameters is equal or lower than the number of dimensions, individual arguments can be single indices or ranges)
	BOOST_MULTI_HD constexpr auto operator()() const& -> const_subarray { return paren_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	BOOST_MULTI_HD constexpr auto operator[]() const& -> const_subarray { return paren_aux_(); }
#endif

 private:
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index_range rng, As... args) const& { return range(rng).rotated().paren_aux_(args...).unrotated(); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(intersecting_range<index> inr, As... args) const& -> decltype(auto) { return paren_aux_(intersection(this->extent(), inr), args...); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args) const& -> decltype(auto) { return operator[](idx).paren_aux_(args...); }

	template<class... As> BOOST_MULTI_HD constexpr auto brckt_aux_(index_range rng, As... args) const& { return range(rng).rotated().paren_aux_(args...).unrotated(); }
	template<class... As> BOOST_MULTI_HD constexpr auto brckt_aux_(intersecting_range<index> inr, As... args) const& -> decltype(auto) { return paren_aux_(intersection(this->extent(), inr), args...); }
	template<class... As> BOOST_MULTI_HD constexpr auto brckt_aux_(index idx, As... args) const& -> decltype(auto) { return operator[](idx).paren_aux_(args...); }

 public:
	// vvv DO NOT remove default parameter `= irange` : the default template parameters below help interpret the expression `{first, last}` syntax as index ranges
	template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                       operator()(A1 arg1) const& -> decltype(auto) { return paren_aux_(arg1); }                                // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator()(A1 arg1, A2 arg2) const& -> decltype(auto) { return paren_aux_(arg1, arg2); }                 // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator()(A1 arg1, A2 arg2, A3 arg3) const& -> decltype(auto) { return paren_aux_(arg1, arg2, arg3); }  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) const& -> decltype(auto) { return paren_aux_(arg1, arg2, arg3, arg4, args...); }

	// clang-format off
	#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	// vvv DO NOT remove default parameter `= irange` : the default template parameters below help interpret the expression `{first, last}` syntax as index ranges
	// template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                    operator[](A1 arg1) const& -> decltype(auto) { return paren_aux_(arg1); }                             // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator[](A1 arg1, A2 arg2) const& -> decltype(auto) { return brckt_aux_(arg1, arg2); }                 // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator[](A1 arg1, A2 arg2, A3 arg3) const& -> decltype(auto) { return brckt_aux_(arg1, arg2, arg3); }  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator[](A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) const& -> decltype(auto) { return brckt_aux_(arg1, arg2, arg3, arg4, args...); }
	#endif
	// clang-format on

 private:
	template<typename Tuple, std::size_t... I> BOOST_MULTI_HD constexpr auto apply_impl_(Tuple const& tuple, std::index_sequence<I...> /*012*/) const& -> decltype(auto) {
		using std::get;
		return this->operator()(get<I>(tuple)...);
	}

 public:
	/// When evaluated on a tuple object this is equivalent to `.operator()(get<0>(tup), get<1>(tup), ...)`. (The argument type typical has `tuple_size<Tuple> == D`)
	template<typename Tuple = typename const_subarray::indices_type> BOOST_MULTI_HD constexpr auto apply(Tuple const& tuple) const& -> decltype(auto) { return apply_impl_(tuple, std::make_index_sequence<std::tuple_size_v<Tuple>>{}); }

	/// Random-access iterator in the leading dimension (return type of `.begin()` or `.end()` from a constant subarray)
	using iterator = detail::array_iterator<element, D, element_ptr, false, false, typename layout_type::stride_type, typename layout_type::sub_type>;  ///< Random access iterator across the leading dimension (e.g. returned by `begin`/`end`)

	/// Random-access const-iterator in the leading dimension (return type of `.begin()` or `.end()` from a constant subarray)
	using const_iterator = detail::array_iterator<element, D, element_ptr, true, false, typename layout_type::stride_type, typename layout_type::sub_type>;  ///< Random access const-iterator across the leading dimension

	const_subarray(const_iterator first, const_iterator last)
	: const_subarray(layout_type(first->layout(), first.stride(), 0, (last - first) * first->size()), first.base()) {
		BOOST_MULTI_ASSERT(first->layout() == last->layout());
	}

 private:
	template<class It>
	friend BOOST_MULTI_HD constexpr auto ref(It begin, It end) -> multi::subarray<typename It::element, It::rank_v, typename It::element_ptr>;

	using const_ptr = detail::subarray_ptr<T, D, ElementPtr, Layout, true>;
	using ptr       = detail::subarray_ptr<T, D, ElementPtr, Layout, false>;

 public:
	/// For `D == 1` it is the pointer type to the elements, otherwise it is void
	using pointer =
		typename std::conditional_t<
			(D > 1),
			void,
			typename std::iterator_traits<ElementPtr>::pointer>;

	/// For `D == 1` it is the pointer type to the elements (as immutables), otherwise it is void
	using const_pointer =
		typename std::conditional_t<
			(D > 1),
			void,
			typename std::iterator_traits<ElementPtr>::pointer>;

 private:
	constexpr auto addressof_aux_() const { return ptr(this->base_, this->layout()); }

	// NOLINTBEGIN(readability-identifier-naming)  // TODO(correaa) rename as private or remove
	constexpr auto addressof_() && -> ptr { return addressof_aux_(); }             // cppcheck-suppress duplInheritedMember;
	constexpr auto addressof_() const&& -> const_ptr { return addressof_aux_(); }  // cppcheck-suppress duplInheritedMember;

	constexpr auto addressof_() & -> ptr { return addressof_aux_(); }             // cppcheck-suppress duplInheritedMember;
	constexpr auto addressof_() const& -> const_ptr { return addressof_aux_(); }  // cppcheck-suppress duplInheritedMember;
	// NOLINTEND(readability-identifier-naming)

 public:
	// NOLINTBEGIN(google-runtime-operator) //NOSONAR
	// operator& is not defined for r-values anyway
	constexpr auto operator&() && { return addressof_(); }       // cppcheck-suppress duplInheritedMember;  //NOSONAR
	constexpr auto operator&() const&& { return addressof_(); }  // cppcheck-suppress duplInheritedMember;  //NOSONAR

	constexpr auto operator&() & { return addressof_(); }       // cppcheck-suppress duplInheritedMember;  //NOSONAR
	constexpr auto operator&() const& { return addressof_(); }  // cppcheck-suppress duplInheritedMember;  //NOSONAR
	// NOLINTEND(google-runtime-operator)

 private:
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	BOOST_MULTI_HD constexpr auto begin_aux_() const { return iterator(types::base_, this->sub(), this->stride()); }
	BOOST_MULTI_HD constexpr auto end_aux_() const { return iterator(types::base_ + this->nelems(), this->sub(), this->stride()); }  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	/// returns an iterator to the beginning (in the leading dimension)
	BOOST_MULTI_HD constexpr auto begin() const& -> const_iterator { return begin_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite  ///< returns a iterator to the beginning

	/// returns an iterator to the end (in the leading dimension)
	BOOST_MULTI_HD constexpr auto end() const& -> const_iterator { return end_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

	/// returns an const-iterator to the beginning
	BOOST_MULTI_HD constexpr auto cbegin() const& { return begin(); }

	/// returns an const-iterator to the end
	BOOST_MULTI_HD constexpr auto cend() const& { return end(); }

	/// Indexable cursor, pointer-like objects that multidimensionally indexable (type usually returned by `.home()`)
	using cursor = detail::cursor_t<typename const_subarray::element_ptr, D, typename const_subarray::strides_type>;

	/// Indexable const-cursor, pointer-like objects that multidimensionally indexable (type usually returned by `.home() const`)
	using const_cursor = detail::cursor_t<typename const_subarray::element_const_ptr, D, typename const_subarray::strides_type>;

 private:
	BOOST_MULTI_HD constexpr auto home_aux_() const { return cursor(this->base_, this->strides()); }

 public:
	/// returns a cursor pointing to the top corner element of the array
	BOOST_MULTI_HD constexpr auto home() const& -> const_cursor { return home_aux_(); }

	template<
		class Range,
		std::enable_if_t<!has_extents<std::decay_t<Range>>::value, int> = 0,
		//  std::enable_if_t<not multi::is_implicitly_convertible_v<subarray, Range>, int> =0,
		class                                                           = decltype(Range(std::declval<typename const_subarray::const_iterator>(), std::declval<typename const_subarray::const_iterator>()))>
	constexpr explicit operator Range() const { return Range(begin(), end()); }  // for example std::vector(it, ti, alloc = {})

	template<class TT, class... As>
	friend constexpr auto operator==(const_subarray const& self, const_subarray<TT, D, As...> const& other) -> bool {
		return (self.extent() == other.extent()) && (self.elements() == other.elements());
	}
	template<class TT, class... As>
	friend constexpr auto operator!=(const_subarray const& self, const_subarray<TT, D, As...> const& other) -> bool {
		return (self.extent() != other.extent()) || (self.elements() != other.elements());
	}

	constexpr auto operator==(const_subarray const& other) const -> bool {
		return (this->extent() == other.extent()) && (this->elements() == other.elements());
	}
	constexpr auto operator!=(const_subarray const& other) const -> bool {
		return (this->extent() != other.extent()) || (this->elements() != other.elements());
	}

	friend constexpr auto lexicographical_compare(const_subarray const& self, const_subarray const& other) -> bool {
		if(self.extent().first() > other.extent().first()) {
			return true;
		}
		if(self.extent().first() < other.extent().first()) {
			return false;
		}
		return adl_lexicographical_compare(
			self.begin(), self.end(),
			other.begin(), other.end()
		);
	}

	constexpr auto operator<(const_subarray const& other) const& -> bool { return lexicographical_compare(*this, other); }
	constexpr auto operator<=(const_subarray const& other) const& -> bool { return *this == other || lexicographical_compare(*this, other); }
	constexpr auto operator>(const_subarray const& other) const& -> bool { return other < *this; }

	/// yields a view of the array in which the internal representation is static_cast to another type (and/or pointer)
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>, std::enable_if_t<std::is_const_v<typename std::pointer_traits<P2>::element_type>, int> = 0  // NOLINT(modernize-use-constraints) TODO(correaa)
			 >
	constexpr auto static_array_cast() const& {                                    // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), static_cast<P2>(this->base_));  // TODO(correaa) might violate constness
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>, std::enable_if_t<!std::is_const_v<typename std::pointer_traits<P2>::element_type>, int> = 0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
			 >
	[[deprecated("violates constness")]]
	constexpr auto static_array_cast() const& {                                    // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), static_cast<P2>(this->base_));  // TODO(correaa) might violate constness
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() && {  // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), static_cast<P2>(this->base_));
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() & {  // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), static_cast<P2>(this->base_));
	}

 private:
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>, class... Args>
	constexpr auto static_array_cast_(Args&&... args) const& {  // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), P2{this->base_, std::forward<Args>(args)...});
	}

 public:
	/// a view with the elements transformed
	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) const& {
		return static_array_cast_<
			//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
			std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
			transform_ptr<
				//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
				std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
				UF, element_const_ptr, std::invoke_result_t<UF const&, element_cref>>>(std::forward<UF>(fun));
	}

	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) & {
		return static_array_cast_<
			std::decay_t<std::invoke_result_t<UF const&, element_ref>>,
			transform_ptr<
				std::decay_t<std::invoke_result_t<UF const&, element_ref>>,
				UF, element_ptr, std::invoke_result_t<UF const&, element_ref>>>(std::forward<UF>(fun));
	}
	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) && { return element_transformed(std::forward<UF>(fun)); }

	/// yields a view of the array containing a specific member of original element type
	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2 const>,
		class Element = typename const_subarray::element,
		class PM      = T2 Element::*>
	constexpr auto member_cast(PM member) const& -> subarray<T2, D, P2> {
		static_assert(sizeof(T) % sizeof(T2) == 0, "array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
												   "Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements.");

		return subarray<T2, D, P2>{this->layout().scale(sizeof(T), sizeof(T2)), static_cast<P2>(&(this->base_->*member))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM      = T2 Element::*>
	constexpr auto member_cast(PM member) & -> subarray<T2, D, P2> {
		static_assert(sizeof(T) % sizeof(T2) == 0, "array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
												   "Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");

		return subarray<T2, D, P2>{this->layout().scale(sizeof(T), sizeof(T2)), static_cast<P2>(&(this->base_->*member))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM      = T2 Element::*>
	constexpr auto member_cast(PM member) && -> subarray<T2, D, P2> {
		return this->member_cast<T2, P2, Element, PM>(member);
	}

 private:
	template<class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>>
	using rebind = subarray<std::decay_t<T2>, D, P2>;

 public:
	/// creates a view of the array with element references with const-removed
	template<
		class T2 = std::remove_const_t<T>,
		class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		std::enable_if_t<    // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
			std::is_same_v<  // check that pointer family is not changed
				typename std::pointer_traits<P2>::template rebind<T2>,
				typename std::pointer_traits<element_ptr>::template rebind<T2>> &&
				std::is_same_v<  // check that only constness is changed
					std::remove_const_t<typename std::pointer_traits<P2>::element_type>, std::remove_const_t<typename const_subarray::element>>,
			int> = 0>
	constexpr auto const_array_cast() const {
		if constexpr(std::is_pointer_v<P2>) {
			return rebind<T2, P2>(this->layout(), const_cast<P2>(this->base_));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		} else {
			return rebind<T2, P2>(this->layout(), reinterpret_cast<P2 const&>(this->base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)  //NOSONAR
		}
	}

	/// yields a const-element view of the subarray, preventing modification of elements.
	constexpr auto as_const() const { return const_subarray(this->layout(), this->base_); }
	// return rebind<element, element_const_ptr>{this->layout(), this->base()};

 private:
	template<class T2, class P2>
	constexpr auto reinterpret_array_cast_aux_() const -> rebind<T2, P2> {
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return rebind<T2, P2>(
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)     // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		);
	}

 public:
	/// yields a view of the subarray where elements are reinterpreted as a different type (elements must have compatible size)
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2 const>>
	constexpr auto reinterpret_array_cast() const& { return reinterpret_array_cast_aux_<T2, P2>().as_const(); }

	template<
		class T2,
		class P2 =
			std::conditional_t<
				std::is_const_v<typename std::pointer_traits<typename const_subarray::element_ptr>::element_type>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2 const>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>>>
	constexpr auto reinterpret_array_cast(size_type count) const& {
		static_assert(sizeof(T) % sizeof(T2) == 0, "error: reinterpret_array_cast is limited to integral stride values");

		BOOST_MULTI_ASSERT(sizeof(T) == sizeof(T2) * static_cast<std::size_t>(count));

		if constexpr(std::is_pointer_v<ElementPtr>) {
			using void_ptr_like = std::conditional_t<
				std::is_const_v<typename std::pointer_traits<decltype(this->base_)>::element_type>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<void const>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<void>>;
			return const_subarray<T2, D + 1, P2>(
				layout_t<D + 1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				static_cast<P2>(static_cast<void_ptr_like>(this->base_))  // NOLINT(bugprone-casting-through-void) direct reinterepret_cast doesn't work here for some exotic pointers (e.g. thrust::pointer)
			);
		} else {  // TODO(correaa) try to unify both if-branches
			return const_subarray<T2, D + 1, P2>(
				layout_t<D + 1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				reinterpret_cast<P2 const&>(this->base_)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
			);
		}
	}

 private:
	/// serializes the subarray elements to an Boost.Serialization-compatible archive
	template<class Archive>                                     // NOLINTNEXTLINE(readability-identifier-naming)
	auto serialize(Archive& arxiv, unsigned int /*version*/) {  // cppcheck-suppress duplInheritedMember ; to overwrite
		using AT = multi::archive_traits<Archive>;
		// if(version == 0) {
		//  std::for_each(this->begin(), this->end(), [&](reference&& item) {arxiv & AT    ::make_nvp("item", std::move(item));});
		// } else {
		std::for_each(this->elements().begin(), this->elements().end(), [&](element const& elem) -> void { arxiv& AT ::make_nvp("elem", elem); });
		// }
		//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv & cereal::make_nvp("item", item);});
		//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv &                          item ;});
	}

	friend class boost::serialization::access;
	friend class cereal::access;
};

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

template<class T>
BOOST_MULTI_HD constexpr auto move(T&& val) noexcept -> decltype(auto) {
	if constexpr(has_member_move<T>::value) {
		return std::forward<T>(val).move();
	} else {
		return std::move(std::forward<T>(val));
	}
}

/// Movable D‐dimensional view into part or all of an array (elements can be moved when dereferenced and assigned)
template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout>
class move_subarray : public subarray<T, D, ElementPtr, Layout> {
	// cppcheck-suppress noExplicitConstructor ; see below
	BOOST_MULTI_HD constexpr move_subarray(subarray<T, D, ElementPtr, Layout>& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  TODO(correa) check if this is necessary
	: subarray<T, D, ElementPtr, Layout>(other.layout(), other.base()) {}

	friend class subarray<T, D, ElementPtr, Layout>;

 public:
	using subarray<T, D, ElementPtr, Layout>::operator[];
	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> decltype(auto) {  // cppcheck-suppress duplInheritedMember ; to overwrite
		return multi::move(subarray<T, D, ElementPtr, Layout>::operator[](idx));
	}

	using subarray<T, D, ElementPtr, Layout>::begin;
	using subarray<T, D, ElementPtr, Layout>::end;

	BOOST_MULTI_HD constexpr auto begin() && { return this->mbegin(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto end() && { return this->mend(); }      // cppcheck-suppress duplInheritedMember ; to overwrite
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout>
class subarray : public const_subarray<T, D, ElementPtr, Layout> {
	// cppcheck-suppress noExplicitConstructor ; see below
	BOOST_MULTI_HD constexpr subarray(const_subarray<T, D, ElementPtr, Layout> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  TODO(correa) check if this is necessary
	: subarray(other.layout(), other.base_) {}

	template<typename, multi::dimensionality_type, typename, class> friend class subarray;
	template<typename, multi::dimensionality_type, typename, class, bool> friend struct detail::subarray_ptr;

	template<class, multi::dimensionality_type, class, bool, bool, typename, class> friend struct detail::array_iterator;

 public:
	using typename const_subarray<T, D, ElementPtr, Layout>::size_type;

	subarray(subarray const&) = delete;

	/// yields a subarray whose elements are marked for move
	BOOST_MULTI_HD constexpr auto move() { return move_subarray<T, D, ElementPtr, Layout>(*this); }

	friend BOOST_MULTI_HD constexpr auto move(subarray& self) { return self.move(); }
	friend BOOST_MULTI_HD constexpr auto move(subarray&& self) { return std::move(self).move(); }

	/// Iterator in the leading dimension that mark elements as movable
	using move_iterator = detail::array_iterator<T, D, ElementPtr, false, true>;

	using typename const_subarray<T, D, ElementPtr, Layout>::element;
	using typename const_subarray<T, D, ElementPtr, Layout>::element_ptr;
	using typename const_subarray<T, D, ElementPtr, Layout>::element_const_ptr;

	/// Subarray reference after binding first index, `multi::subarray<T, D - 1, P>` or, for `D == 1`, `std::pointer_traits<P>::reference` (usually `T&`)
	using reference = typename std::conditional_t<
		(D > 1),
		subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_ptr>::reference>;

	/// Subarray immutable reference after binding first index, `multi::const_subarray<T, D - 1, P>` or, for `D == 1`, `std::pointer_traits<P>::reference` (usually `T const&`)
	using const_reference = typename std::conditional_t<
		(D > 1),
		const_subarray<element, D - 1, element_ptr>,
		typename std::iterator_traits<element_const_ptr>::reference>;

	subarray(subarray&&) noexcept = default;
	~subarray()                   = default;

	template<
		class Other,
		class                                                                                                  = std::enable_if_t<!std::is_base_of_v<subarray, Other> && !std::is_base_of_v<Other, subarray>>,  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-type-traits)  TODO(correaa) in C++20
		decltype(multi::detail::implicit_cast<typename subarray::element_ptr>(typename Other::element_ptr{}))* = nullptr,
		decltype(std::declval<Other const&>().base())*                                                         = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr /*mplct*/ subarray(Other&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor,cppcoreguidelines-missing-std-forward) : to reproduce the implicitness of the argument
	: subarray(other.layout(), other.base()) {}

 private:
	using ptr = void;  // detail::subarray_ptr<T, D, ElementPtr, Layout, false>;

 public:
	// clang-format off
	#ifdef __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wlarge-by-value-copy"  // good to know
	#endif

	// cppcheck-suppress duplInheritedMember ; to overwrite  // NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() && { return detail::subarray_ptr<T, D, ElementPtr, Layout, false>(this->base_, this->layout()); }  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR
	// cppcheck-suppress duplInheritedMember ; to overwrite  // NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() & { return detail::subarray_ptr<T, D, ElementPtr, Layout, false>(this->base_, this->layout()); }  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR

	#ifdef __clang__
	#pragma clang diagnostic pop
	#endif
	// clang-format on

	using const_subarray<T, D, ElementPtr, Layout>::operator&;

	using const_subarray<T, D, ElementPtr, Layout>::const_subarray;

	using const_subarray<T, D, ElementPtr, Layout>::elements;
	constexpr auto                         elements() & { return this->elements_aux_(); }   // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_NO_DANGLING constexpr auto elements() && { return this->elements_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::flattened;
	constexpr auto flattened() & { return this->flattened_aux_(); }   // cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto flattened() && { return this->flattened_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::begin;
	// cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto begin() && noexcept { return this->begin_aux_(); }
	BOOST_MULTI_HD constexpr auto begin() & noexcept { return this->begin_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::end;
	BOOST_MULTI_HD constexpr auto end() && noexcept { return this->end_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	// cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto end() & noexcept { return this->end_aux_(); }

	/// returns an move-iterator (moves on dereference) to the beginning in the leading dimension
	BOOST_MULTI_HD constexpr auto mbegin() { return move_iterator{this->begin()}; }
	/// returns an move-iterator (moves on dereference) to the ending in the leading dimension
	BOOST_MULTI_HD constexpr auto mend() { return move_iterator{this->end()}; }

	using const_subarray<T, D, ElementPtr, Layout>::home;
	BOOST_MULTI_HD constexpr auto home() && { return this->home_aux_(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto home() & { return this->home_aux_(); }   // cppcheck-suppress duplInheritedMember ; to overwrite

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
	operator std::mdspan<T const, std::dextents<std::size_t, D>, std::layout_stride>() const& { return this->to_mdspan_aux_(); }
	operator std::mdspan<T, std::dextents<std::size_t, D>, std::layout_stride>() & { return this->to_mdspan_aux_(); }
	operator std::mdspan<T, std::dextents<std::size_t, D>, std::layout_stride>() && { return this->to_mdspan_aux_(); }
#endif

	// /// assigns `.size()` values after an iterator into the array.
	// template<class It> constexpr auto assign(It first) & -> It {  // cppcheck-suppress duplInheritedMember ; to overwrite
	// 	adl_copy_n(first, this->size(), begin());
	// 	std::advance(first, this->size());
	// 	return first;
	// }
	// template<class It> BOOST_MULTI_HD constexpr auto assign(It first) && -> It { return assign(first); }  // cppcheck-suppress duplInheritedMember ; to overwrite

	// template<class TT = typename subarray::element>
	// constexpr auto fill(TT const& value) & -> decltype(auto) {
	// 	return adl_fill_n(this->begin(), this->size(), value), *this;
	// }
	// constexpr auto fill() & -> decltype(auto) { return fill(typename subarray::element{}); }

	// template<class TT = typename subarray::element>
	// [[deprecated]] constexpr auto fill(TT const& value) && -> decltype(auto) { return std::move(this->fill(value)); }
	// [[deprecated]] constexpr auto fill() && -> decltype(auto) {
	// 	return std::move(*this).fill(typename subarray::element{});
	// }

	using const_subarray<T, D, ElementPtr, Layout>::strided;
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	constexpr auto strided(difference_type step) && { return this->strided_aux_(step); }
	constexpr auto strided(difference_type step) & { return this->strided_aux_(step); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::taked;
	constexpr auto taked(difference_type count) && -> subarray { return this->taked_aux_(count); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto taked(difference_type count) & -> subarray { return this->taked_aux_(count); }   // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::dropped;
	// cppcheck-suppress-begin duplInheritedMember ; to ovewrite
	constexpr auto dropped(difference_type count) && -> subarray { return this->dropped_aux_(count); }
	constexpr auto dropped(difference_type count) & -> subarray { return this->dropped_aux_(count); }
	// cppcheck-suppress-end duplInheritedMember ; to ovewrite

	using const_subarray<T, D, ElementPtr, Layout>::rotated;
	// cppcheck-suppress-begin duplInheritedMember ; to ovewrite
	BOOST_MULTI_HD constexpr auto rotated() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::rotated(); }
	BOOST_MULTI_HD constexpr auto rotated() & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::rotated(); }
	// cppcheck-suppress-end duplInheritedMember ; to ovewrite

	using const_subarray<T, D, ElementPtr, Layout>::unrotated;
	BOOST_MULTI_HD constexpr auto unrotated() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unrotated(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto unrotated() & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unrotated(); }   // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::transposed;
	BOOST_MULTI_HD constexpr auto transposed() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::transposed(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto transposed() & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::transposed(); }   // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::unordered;
	BOOST_MULTI_HD constexpr auto unordered() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unordered(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto unordered() & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unordered(); }   // cppcheck-suppress duplInheritedMember ; to overwrite

	// BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	// auto operator~ (subarray const& self) { return self.transposed(); }
	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD auto operator~(subarray& self) { return self.transposed(); }
	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD auto operator~(subarray&& self) { return std::move(self).transposed(); }

	using const_subarray<T, D, ElementPtr, Layout>::reindexed;

	template<class... Indexes>
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto reindexed(index first, Indexes... idxs) & -> subarray {
		return const_subarray<T, D, ElementPtr, Layout>::reindexed(first, idxs...);
		// return ((this->reindexed(first).rotated()).reindexed(idxs...)).unrotated();
	}
	template<class... Indexes>
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto reindexed(index first, Indexes... idxs) && -> subarray {
		return const_subarray<T, D, ElementPtr, Layout>::reindexed(first, idxs...);
		// return ((std::move(*this).reindexed(first).rotated()).reindexed(idxs...)).unrotated();
	}

	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto base() const& -> typename subarray::element_const_ptr { return this->base_; }
	BOOST_MULTI_HD constexpr auto base() & -> ElementPtr { return this->base_; }
	BOOST_MULTI_HD constexpr auto base() && -> ElementPtr { return this->base_; }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto operator=(const_subarray<T, D, ElementPtr, Layout> const& other) & -> subarray& {
		if(this == std::addressof(other)) {
			return *this;
		}
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		this->elements() = other.elements();
		return *this;
	}

	/// swaps every corresponding element of the array references, extents must match. O(N) operation
	constexpr void swap(subarray& other) & noexcept {
		BOOST_MULTI_ASSERT(this->extents() == other.extents());
		adl_swap_ranges(this->elements().begin(), this->elements().end(), other.elements().begin());
	}

	constexpr void swap(subarray&& other) & noexcept {
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		adl_swap_ranges(this->elements().begin(), this->elements().end(), std::move(other).elements().begin());
	}

	constexpr void swap(subarray&& other) && noexcept { return swap(std::move(other)); }
	constexpr void swap(subarray& other) && noexcept { return swap(other); }

	friend constexpr void swap(subarray&& self, subarray&& other) noexcept { std::move(self).swap(std::move(other)); }
	friend constexpr void swap(subarray&& self, subarray& other) noexcept { std::move(self).swap(other); }
	friend constexpr void swap(subarray& self, subarray&& other) noexcept { self.swap(std::move(other)); }
	friend constexpr void swap(subarray& self, subarray& other) noexcept { self.swap(other); }

	// fix mutation
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...> const& other) && -> decltype(auto) {operator=(          other ); return *this;}
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...> const& other) & -> subarray& {
		BOOST_MULTI_ASSERT(other.extents() == this->extents());
		this->elements() = other.elements();
		return *this;
	}

	// fix mutation
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...>&& other) && -> subarray& {
		operator=(std::move(other));
		return *this;
	}
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...>&& other) & -> subarray& {
		BOOST_MULTI_ASSERT(this->extents() == other.extents());
		this->elements() = std::move(other).elements();
		return *this;
	}

	template<
		class Range,
		class                                              = std::enable_if_t<!std::is_base_of_v<subarray, Range>>,  // NOLINT(modernize-type-traits)  TODO(correaa) in C++20
		class                                              = std::enable_if_t<!detail::is_subarray<Range>::value>,   // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
		std::enable_if_t<!has_elements<Range>::value, int> = 0>                                                      // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator=(Range const& rng) & -> subarray& {                                                      // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		BOOST_MULTI_ASSERT(this->size() == static_cast<size_type>(adl_size(rng)));                                   // TODO(correaa) or use std::cmp_equal?
		if(adl_size(rng)) {
			adl_copy_n(adl_begin(rng), adl_size(rng), this->begin());
		}
		return *this;
	}

	template<
		class MultiRange,
		class                                                  = std::enable_if_t<!std::is_base_of_v<subarray, MultiRange>>,  // NOLINT(modernize-type-traits)  TODO(correaa) in C++20
		class                                                  = std::enable_if_t<!detail::is_subarray<MultiRange>::value>,   // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
		std::enable_if_t<has_elements<MultiRange>::value, int> = 0>                                                           // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator=(MultiRange const& mrng) & -> subarray& {                                                         // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		BOOST_MULTI_ASSERT(this->extents() == mrng.extents());                                                                // TODO(correaa) or use std::cmp_equal?
		adl_copy_n(mrng.elements().begin(), this->num_elements(), this->elements().begin());
		return *this;
	}

	template<
		class Range,
		class = std::enable_if_t<!std::is_base_of_v<subarray, Range>>,  // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
		class = std::enable_if_t<!detail::is_subarray<Range>::value>    // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
		>
	constexpr auto operator=(Range const& rng) && -> subarray& {
		operator=(rng);
		return *this;
	}

	template<class TT, class... As>
	constexpr auto operator=(const_subarray<TT, D, As...> const& other) && -> subarray& {
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		this->elements() = other.elements();
		return *this;
	}

	template<class TT, class... As>
	constexpr auto operator=(subarray<TT, D, As...>&& other) & -> subarray& {
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		this->elements() = std::move(other).elements();
		return *this;
	}

	/// Assignment operators (right-hand side must be of the same dimensionality)
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto operator=(const_subarray<T, D, ElementPtr, Layout> const& other) const&& -> subarray&;  // for std::indirectly_writable

	constexpr auto operator=(subarray const& other) & -> subarray& {
		if(this == std::addressof(other)) {
			return *this;
		}
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		this->elements() = other.elements();
		return *this;
	}

	constexpr auto operator=(subarray&& other) & noexcept(false) -> subarray& {  // NOLINT(bugprone-unsafe-to-allow-exceptions) TODO(correaa) make conditionally noexcept
		// if(this == std::addressof(other)) { return *this; }
		BOOST_MULTI_ASSERT(this->extent() == other.extent());
		this->elements() = std::move(other).elements();
		return *this;
	}

	auto operator=(std::initializer_list<typename subarray::value_type> values) && -> subarray& {
		operator=(values);
		return *this;
	}
	auto operator=(std::initializer_list<typename subarray::value_type> values) & -> subarray& {
		BOOST_MULTI_ASSERT(static_cast<size_type>(values.size()) == this->size());
		if(values.size() != 0) {
			adl_copy_n(values.begin(), values.size(), this->begin());
		}
		return *this;
	}

	// BOOST_MULTI_HD constexpr auto operator[](index idx) const&    { return static_cast<typename subarray::const_reference>(this->at_aux_(idx)); }  // TODO(correaa) use return type to cast
	using const_subarray<T, D, ElementPtr, Layout>::operator[];
	// BOOST_MULTI_HD constexpr auto operator[](index idx) const& { return const_subarray<T, D, ElementPtr, Layout>::operator[](idx); }

	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> typename subarray::reference { return this->at_aux_(idx); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto operator[](index idx) & -> typename subarray::reference { return this->at_aux_(idx); }   // cppcheck-suppress duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::diagonal;

	// template<class Dummy = void, std::enable_if_t<(D > 1) && sizeof(Dummy*), int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	// cppcheck-suppress-begin duplInheritedMember ; to override
	constexpr auto diagonal() & { return this->diagonal_aux_(); }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0
	constexpr auto diagonal() && { return this->diagonal_aux_(); }
	// cppcheck-suppress-end duplInheritedMember ; to override

	using const_subarray<T, D, ElementPtr, Layout>::sliced;
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::sliced(first, last); }
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::sliced(first, last); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

 private:
	using const_subarray<T, D, ElementPtr, Layout>::range;
	// NOLINTNEXTLINE(readability-identifier-naming)
	BOOST_MULTI_HD constexpr auto range(index_range irng) && -> decltype(auto) { return std::move(*this).sliced(irng.front(), irng.front() + irng.size()); }  // cppcheck-suppress duplInheritedMember;
	// NOLINTNEXTLINE(readability-identifier-naming)
	BOOST_MULTI_HD constexpr auto range(index_range irng) & -> decltype(auto) { return sliced(irng.front(), irng.front() + irng.size()); }  // cppcheck-suppress duplInheritedMember;

	using const_subarray<T, D, ElementPtr, Layout>::paren_aux_;

	BOOST_MULTI_HD constexpr auto paren_aux_() & { return subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_); }  // cppcheck-suppress duplInheritedMember;
	BOOST_MULTI_HD constexpr auto paren_aux_() && { return subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_); }

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx) & -> decltype(auto) { return operator[](idx); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx) && -> decltype(auto) { return operator[](idx); }

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args) & -> decltype(auto) { return operator[](idx).paren_aux_(args...); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args) && -> decltype(auto) { return operator[](idx).paren_aux_(args...); }

	template<class... As>
	BOOST_MULTI_HD constexpr auto paren_aux_(index_range irng, As... args) & {
		return this->range(irng).rotated().paren_aux_(args...).unrotated();
	}
	template<class... As>
	BOOST_MULTI_HD constexpr auto paren_aux_(index_range irng, As... args) && {
		return std::move(*this).range(irng).rotated().paren_aux_(args...).unrotated();
	}

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(intersecting_range<index> inr, As... args) & -> decltype(auto) { return paren_aux_(intersection(this->extent(), inr), args...); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(intersecting_range<index> inr, As... args) && -> decltype(auto) { return paren_aux_(intersection(this->extent(), inr), args...); }

 public:
	using const_subarray<T, D, ElementPtr, Layout>::operator();
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto operator()() & -> subarray { return this->paren_aux_(); }
	BOOST_MULTI_HD constexpr auto operator()() && -> subarray { return this->paren_aux_(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                       operator()(A1 arg1) & -> decltype(auto) { return this->paren_aux_(arg1); }
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator()(A1 arg1, A2 arg2) & -> decltype(auto) { return this->paren_aux_(arg1, arg2); }
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator()(A1 arg1, A2 arg2, A3 arg3) & -> decltype(auto) { return this->paren_aux_(arg1, arg2, arg3); }
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) & -> decltype(auto) { return this->paren_aux_(arg1, arg2, arg3, arg4, args...); }

	template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                       operator()(A1 arg1) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1); }                                                                    // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator()(A1 arg1, A2 arg2) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2); }                                                     // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator()(A1 arg1, A2 arg2, A3 arg3) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2, arg3); }                                      // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2, arg3, arg4, args...); }  // NOLINT(whitespace/line_length) pattern line

	// clang-format off
	#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	// template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                       operator[](A1 arg1) & -> decltype(auto) { return this->paren_aux_(arg1); }
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator[](A1 arg1, A2 arg2) & -> decltype(auto) { return this->paren_aux_(arg1, arg2); }
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator[](A1 arg1, A2 arg2, A3 arg3) & -> decltype(auto) { return this->paren_aux_(arg1, arg2, arg3); }
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator[](A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) & -> decltype(auto) { return this->paren_aux_(arg1, arg2, arg3, arg4, args...); }

	// template<class A1 = irange> BOOST_MULTI_HD constexpr auto                                                                       operator[](A1 arg1) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1); }                                                                    // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange> BOOST_MULTI_HD constexpr auto                                                    operator[](A1 arg1, A2 arg2) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2); }                                                     // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange> BOOST_MULTI_HD constexpr auto                                 operator[](A1 arg1, A2 arg2, A3 arg3) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2, arg3); }                                      // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator[](A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) && -> decltype(auto) { return std::move(*this).paren_aux_(arg1, arg2, arg3, arg4, args...); }  // NOLINT(whitespace/line_length) pattern line
	#endif
	// clang-format on
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

 private:
	template<class Self, typename Tuple, std::size_t... I>
	static BOOST_MULTI_HD constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...> /*012*/) -> decltype(auto) {
		using std::get;  // for C++17 compatibility
		return std::forward<Self>(self)(get<I>(tuple)...);
	}

 public:
	using const_subarray<T, D, ElementPtr, Layout>::apply;
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	template<typename Tuple = typename subarray::indices_type> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl) && -> decltype(auto) { return apply_impl_(std::move(*this), tpl, std::make_index_sequence<std::tuple_size_v<Tuple>>()); }
	template<typename Tuple = typename subarray::indices_type> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl) & -> decltype(auto) { return apply_impl_(*this, tpl, std::make_index_sequence<std::tuple_size_v<Tuple>>()); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	using const_subarray<T, D, ElementPtr, Layout>::partitioned;
	// cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto partitioned(size_type size) & -> subarray<T, D + 1, typename subarray::element_ptr> { return this->partitioned_aux_(size); }

	// cppcheck-suppress duplInheritedMember ; to overwrite
	BOOST_MULTI_HD constexpr auto partitioned(size_type size) && -> subarray<T, D + 1, typename subarray::element_ptr> { return this->partitioned_aux_(size); }

	// using const_subarray<T, D, ElementPtr, Layout>::flatted;
	constexpr auto flatted() const& {
		assert(this->layout().is_flattable());
		multi::layout_t<D - 1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return const_subarray<T, D - 1, ElementPtr>{new_layout, this->base_};
	}

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto flatted() & {
		multi::layout_t<D - 1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return subarray<T, D - 1, ElementPtr>(new_layout, this->base_);
	}

	constexpr auto flatted() && { return this->flatted(); }  // cppcheck-suppress duplInheritedMember ; to override

	using const_subarray<T, D, ElementPtr, Layout>::reinterpret_array_cast;

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast() & {  // cppcheck-suppress duplInheritedMember ; to override
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return subarray<T2, D, P2>(
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)     // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		);
	}

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto reinterpret_array_cast() && {
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return subarray<T2, D, P2>(
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)     // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		);
	}

 private:
	template<typename P2>
	constexpr static auto reinterpret_pointer_cast_(ElementPtr const& ptr) -> decltype(auto) {
		if constexpr(std::is_pointer_v<ElementPtr>) {
			return static_cast<P2>(static_cast<void*>(ptr));  // NOLINT(bugprone-casting-through-void) direct reinterepret_cast doesn't work here
		} else {
			return reinterpret_cast<P2 const&>(ptr);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
		}
	}

 public:
	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto reinterpret_array_cast(size_type count) & {
		static_assert(sizeof(T) % sizeof(T2) == 0, "error: reinterpret_array_cast is limited to integral stride values");

		BOOST_MULTI_ASSERT(sizeof(T) == sizeof(T2) * static_cast<std::size_t>(count));

		layout_t<D + 1> const lyt1{this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count};
		auto const            lyt2 = lyt1.rotate();
		return subarray<T2, D + 1, P2>(
			lyt2,
			reinterpret_pointer_cast_<P2>(this->base_)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
		);
	}

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto reinterpret_array_cast(size_type count) && {
		static_assert(sizeof(T) % sizeof(T2) == 0, "error: reinterpret_array_cast is limited to integral stride values");

		BOOST_MULTI_ASSERT(sizeof(T) == sizeof(T2) * static_cast<std::size_t>(count));

		return subarray<T2, D + 1, P2>(
			layout_t<D + 1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
			reinterpret_pointer_cast_<P2>(this->base_)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
		);
	}

 private:
	using element_move_ptr = multi::move_ptr<typename subarray::element, typename subarray::element_ptr>;

 public:
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	/// An array view in which element references are r-values
	constexpr auto element_moved() & { return subarray<T, D, typename subarray::element_move_ptr, Layout>(this->layout(), element_move_ptr{this->base_}); }
	constexpr auto element_moved() && { return element_moved(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	/// serializes data to generic archive (e.g. Boost or Cereal archive)
	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int /*version*/) {  // cppcheck-suppress duplInheritedMember ; to override
		using AT = multi::archive_traits<Archive>;
		// if(version == 0) {
		//  std::for_each(this->begin(), this->end(), [&](typename subarray::reference item) {arxiv & AT    ::make_nvp("item", item);});
		// } else {
		std::for_each(this->elements().begin(), this->elements().end(), [&](typename subarray::element& elem) -> void { arxiv& AT ::make_nvp("elem", elem); });
		//}
		//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv & cereal::make_nvp("item", item);});
		//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv &                          item ;});
	}
};

// template<class Subarray> auto diagonal(Subarray&& sarr)
// 	-> decltype(std::forward<Subarray>(sarr).diagonal()) {
// 	return std::forward<Subarray>(sarr).diagonal();
// }

namespace detail {

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class Element, typename Ptr> struct array_iterator<Element, 0, Ptr> {};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<class Element, typename Ptr, bool IsConst, bool IsMove, typename Stride>
struct array_iterator<Element, 1, Ptr, IsConst, IsMove, Stride>  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) stride_ is not initialized in some constructors
: boost::multi::iterator_facade<
	  array_iterator<Element, 1, Ptr, IsConst, IsMove, Stride>,
	  Element, std::random_access_iterator_tag,
	  std::conditional_t<
		  IsConst,
		  typename std::iterator_traits<typename std::pointer_traits<Ptr>::template rebind<Element const>>::reference,
		  std::conditional_t<
			  IsMove,
			  std::add_rvalue_reference_t<std::decay_t<typename std::iterator_traits<Ptr>::reference>>,
			  typename std::iterator_traits<Ptr>::reference>>,
	  multi::difference_type>  //
{
	using affine = multi::affine<detail::array_iterator<Element, 1, Ptr, IsConst, IsMove, Stride>, multi::difference_type>;

	using pointer = std::conditional_t<
		IsConst,
		typename std::pointer_traits<Ptr>::template rebind<Element const>,
		Ptr>;

	auto segment() const {
		return subarray<Element, 1, Ptr>(
			layout_t<1>(
				layout_t<0>(extents_t<0>{}),
				stride().stride2(),
				0,
				stride().nelems2()
			),  // scalar sub-layout (num_elements()==1); a default `{}` would have num_elements()==0 and break `elements()`/`operator==`
			this->stride().segment_base(this->ptr_)
		);
	}

	// using segment_iterator = ::boost::multi::detail::array_iterator<Element, 2, Ptr, IsConst, IsMove, typename Stride::stride2_type>;

	// using segment_iterator = array_iterator<Element, 2, Ptr, IsConst, IsMove, typename Stride::stride2_type>;

	auto outer() const {
		return array_iterator<Element, 2, Ptr, IsConst, IsMove, typename Stride::stride2_type>(
			this->stride().segment_base(this->ptr_),
			layout_t<1>(
				layout_t<0>(extents_t<0>{}),
				stride().stride2(),
				0,
				stride().nelems2()
			),
			stride().stride1()
		);
	}

	auto local() {
		return array_iterator<Element, 1, Ptr, IsConst, IsMove, typename Stride::stride2_type>(
			this->base(), layout_t<0>(extents_t<0>{}), stride().stride2()
		);
	}

 private:
	using reference_aux = std::conditional_t<
		IsConst,
		typename std::iterator_traits<typename std::pointer_traits<Ptr>::template rebind<Element const>>::reference,
		typename std::iterator_traits<Ptr>::reference>;

 public:
	using stride_type = Stride;  // multi::index;

	using reference = std::conditional_t<
		IsMove,
		std::add_rvalue_reference_t<std::decay_t<reference_aux>>,
		reference_aux>;

	using difference_type = typename affine::difference_type;

	using iterator_category = typename stride_traits<Stride>::category;
	using iterator_concept  = typename stride_traits<Stride>::category;

	using element      = Element;
	using element_type = typename std::pointer_traits<Ptr>::element_type;  // workaround for clang 15 and libc++ in c++20 mode

	template<class Element2>
	using rebind = array_iterator<std::decay_t<Element2>, 1, typename std::pointer_traits<Ptr>::template rebind<Element2>, IsConst, IsMove, Stride>;

	static constexpr dimensionality_type dimensionality = 1;

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && (!defined(__clang__) || __clang_major__ != 10)
	// template<class T = void,
	//  std::enable_if_t<sizeof(T*) && std::is_base_of_v<std::contiguous_iterator_tag, iterator_category>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	constexpr explicit operator Ptr() const& {
		static_assert(std::is_base_of_v<std::contiguous_iterator_tag, iterator_category>, "iterator must be continuous");
		return ptr_;
	}
#endif

	array_iterator()  = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
	using layout_type = multi::layout_t<0>;

	template<
		bool OtherIsConst, std::enable_if_t<!OtherIsConst, int> = 0  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
		>
	// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	BOOST_MULTI_HD constexpr array_iterator(array_iterator<Element, 1, Ptr, OtherIsConst, IsMove, Stride> const& other)
	: ptr_{other.base()}, stride_{other.stride()} {}

	template<
		class Other,
		decltype(multi::detail::implicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().base())*                          = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr /*mplct*/ array_iterator(Other const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to reproduce the implicitness of the argument
	: ptr_{other.base()}, stride_{other.stride()} {}

	template<
		class Other,
		decltype(multi::detail::explicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().data_)*                           = nullptr>
	constexpr explicit array_iterator(Other const& other)
	: ptr_{other.data_}, stride_{other.stride_} {}

	template<class, dimensionality_type, class, bool, bool, typename, class> friend struct array_iterator;

	template<
		class EElement, typename PPtr,
		typename = decltype(multi::detail::implicit_cast<Ptr>(std::declval<array_iterator<EElement, 1, PPtr>>().data_))>
	BOOST_MULTI_HD constexpr /*impl*/ array_iterator(array_iterator<EElement, 1, PPtr> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to reproduce the implicitness of original pointer
	: ptr_{other.base()}, stride_{other.stride_} {}

	BOOST_MULTI_HD constexpr explicit operator bool() const { return static_cast<bool>(this->ptr_); }

	BOOST_MULTI_HD constexpr auto operator[](typename array_iterator::difference_type n) const -> decltype(auto) {
		return *((*this) + n);
	}

	constexpr auto operator->() const { return static_cast<pointer>(ptr_); }

	using element_ptr = Ptr;

	static constexpr dimensionality_type rank_v = 1;

	using rank = std::integral_constant<dimensionality_type, rank_v>;

	BOOST_MULTI_HD constexpr explicit array_iterator(typename subarray<element, 0, element_ptr>::element_ptr base, layout_t<0> const& /*lyt*/, Stride stride)
	: ptr_(std::move(base) /*, lyt*/), stride_{stride} {}

 private:
	friend class const_subarray<Element, 1, Ptr>;  // TODO(correaa) fix template parameters

	element_ptr ptr_;
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820)  // warning C4820:  '7' bytes padding added after data member 'boost::multi::array_types<T,2,ElementPtr,Layout>::base_' [C:\Gitlab-Runner\builds\t3_1sV2uA\0\correaa\boost-multi\build\test\array_fancyref.cpp.x.vcxproj]
#endif
	BOOST_MULTI_NO_UNIQUE_ADDRESS
	stride_type stride_;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

 public:
	BOOST_MULTI_HD constexpr auto operator+(difference_type n) const { return array_iterator{*this} += n; }
	BOOST_MULTI_HD constexpr auto operator-(difference_type n) const { return array_iterator{*this} -= n; }

	BOOST_MULTI_HD constexpr auto base() const { return static_cast<pointer>(ptr_); }

	[[deprecated("use base() for iterator")]]
	BOOST_MULTI_HD constexpr auto data() const { return base(); }

	BOOST_MULTI_FRIEND_CONSTEXPR auto base(array_iterator const& self) { return self.base(); }

	BOOST_MULTI_HD constexpr auto stride() const -> stride_type { return stride_; }
	friend constexpr auto         stride(array_iterator const& self) -> stride_type { return self.stride_; }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
	BOOST_MULTI_HD constexpr auto operator++() -> array_iterator& {
		ptr_ += stride_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator--() -> array_iterator& {
		ptr_ -= stride_;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> array_iterator& {
		ptr_ = ptr_ + stride_ * n;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator-=(difference_type n) -> array_iterator& {
		ptr_ -= stride_ * n;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_HD constexpr auto operator-(array_iterator const& other) const -> difference_type {
		BOOST_MULTI_ASSERT(stride() != 0);
		BOOST_MULTI_ASSERT(stride() == other.stride());
		BOOST_MULTI_ASSERT((ptr_ - other.ptr_) % stride() == 0);
		return (ptr_ - other.ptr_) / stride();  // with struct-overflow=3 error: assuming signed overflow does not occur when simplifying `X - Y > 0` to `X > Y` [-Werror=strict-overflow]
	}

	BOOST_MULTI_HD constexpr auto operator==(array_iterator const& other) const noexcept {
		BOOST_MULTI_ASSERT(this->stride_ == other.stride_);
		return this->ptr_ == other.ptr_;
	}

	BOOST_MULTI_HD constexpr auto operator!=(array_iterator const& other) const noexcept {
		BOOST_MULTI_ASSERT(this->stride_ == other.stride_);
		return this->ptr_ != other.ptr_;
	}

	template<bool OtherIsConst, std::enable_if_t<OtherIsConst != IsConst, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr auto operator==(array_iterator<Element, 1, Ptr, OtherIsConst> const& other) const -> bool {
		BOOST_MULTI_ASSERT(this->stride_ == other.stride_);
		BOOST_MULTI_ASSERT(stride_ != 0);
		return this->ptr_ == other.ptr_;
	}

	template<bool OtherIsConst, std::enable_if_t<OtherIsConst != IsConst, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD constexpr auto operator!=(array_iterator<Element, 1, Ptr, OtherIsConst> const& other) const -> bool {
		return !operator==(other);
	}

	BOOST_MULTI_HD constexpr auto operator<(array_iterator const& other) const -> bool {
		return 0 < other - *this;
	}

	BOOST_MULTI_HD constexpr auto operator*() const noexcept -> reference {
		return static_cast<reference>(*ptr_);
	}

	// BOOST_MULTI_HD constexpr auto segment() const -> segment_type {
	// }
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

}  // end namespace detail

// template<class Element, dimensionality_type D, typename... Ts>
// using iterator = array_iterator<Element, D, Ts...>;

/// ... a specialization for zero dimensions
template<typename T, typename ElementPtr, class Layout>
class const_subarray<T, 0, ElementPtr, Layout>
: public detail::array_types<T, 0, ElementPtr, Layout> {
 public:
	using types = detail::array_types<T, 0, ElementPtr, Layout>;
	using types::types;

	using element_type [[deprecated("use ::element")]] = typename types::element;
	using element                                      = typename types::element;
	using element_ref                                  = typename std::iterator_traits<typename const_subarray::element_ptr>::reference;
	using element_cref                                 = typename std::iterator_traits<typename const_subarray::element_const_ptr>::reference;
	using iterator                                     = detail::array_iterator<T, 0, ElementPtr>;

	using size_type = typename detail::array_types<T, 0, ElementPtr, Layout>::size_type;

	using layout_type = Layout;

	constexpr auto operator=(element const& elem) & -> const_subarray& {
		//  MULTI_MARK_SCOPE(std::string{"multi::operator= D=0 from "}+typeid(T).name()+" to "+typeid(T).name() );
		adl_copy_n(&elem, 1, this->base_);
		return *this;
	}
	constexpr auto operator=(element const& elem) && -> const_subarray& {
		operator=(elem);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	constexpr auto operator==(element const& elem) const -> bool {
		BOOST_MULTI_ASSERT(this->num_elements() == 1);
		return adl_equal(&elem, std::next(&elem, this->num_elements()), this->base());
	}
	constexpr auto operator!=(element const& elem) const { return !operator==(elem); }

	template<class Range0>
	constexpr auto operator=(Range0 const& rng) & -> const_subarray& {
		adl_copy_n(&rng, 1, this->base_);
		return *this;
	}

	// constexpr auto elements_at(size_type idx [[maybe_unused]]) const& -> element_cref {
	// 	BOOST_MULTI_ASSERT(idx < this->num_elements());
	// 	return *(this->base_);
	// }
	// constexpr auto elements_at(size_type idx [[maybe_unused]]) && -> element_ref {
	// 	BOOST_MULTI_ASSERT(idx < this->num_elements());
	// 	return *(this->base_);
	// }
	// constexpr auto elements_at(size_type idx [[maybe_unused]]) & -> element_ref {
	// 	BOOST_MULTI_ASSERT(idx < this->num_elements());
	// 	return *(this->base_);
	// }

	constexpr auto operator!=(const_subarray const& other) const { return !adl_equal(other.base_, other.base_ + 1, this->base_); }
	constexpr auto operator==(const_subarray const& other) const { return adl_equal(other.base_, other.base_ + 1, this->base_); }

	constexpr auto operator<(const_subarray const& other) const {
		return adl_lexicographical_compare(
			this->base_, this->base_ + this->num_elements(),
			other.base_, other.base_ + other.num_elements()
		);
	}

	using decay_type = typename types::element;

	BOOST_MULTI_HD constexpr auto operator()() const& -> element_ref { return *(this->base_); }  // NOLINT(hicpp-explicit-conversions)

	constexpr operator element_ref() && noexcept { return *(this->base_); }       // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax
	constexpr operator element_ref() & noexcept { return *(this->base_); }        // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax
	constexpr operator element_cref() const& noexcept { return *(this->base_); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax

	constexpr auto elements() const&;

	constexpr auto begin() const& = delete;
	constexpr auto end() const&   = delete;

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	constexpr auto operator[]() const& = delete;
#else
	template<class IndexType>
	constexpr auto operator[](IndexType const&) const& = delete;
#endif

	auto           diagonal() const     = delete;
	constexpr auto sliced() const&      = delete;
	constexpr auto partitioned() const& = delete;
	constexpr auto flattened() const    = delete;

	constexpr auto strided(difference_type) const& = delete;

	constexpr auto taked(difference_type) const&   = delete;
	constexpr auto dropped(difference_type) const& = delete;

	BOOST_MULTI_HD constexpr auto reindexed() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto rotated() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto unrotated() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto unordered() const& { return operator()(); }

	auto transposed() const&              = delete;
	auto flatted() const&                 = delete;
	auto range() const& -> const_subarray = delete;

	// a lightweigh type for multidimensional indexing, it has the indexing interface of an array but without size (extents) information
	using cursor       = detail::cursor_t<typename const_subarray::element_ptr, 0, typename const_subarray::strides_type>;
	// a lightweigh type for multidimensional indexing (const version), it has the indexing interface of an array but without size (extents) information
	using const_cursor = detail::cursor_t<typename const_subarray::element_const_ptr, 0, typename const_subarray::strides_type>;

 private:
	BOOST_MULTI_HD constexpr auto home_aux_() const { return cursor(this->base_, this->strides()); }

 public:
	BOOST_MULTI_HD constexpr auto home() const& -> const_cursor { return home_aux_(); }

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	auto paren_aux_() const& { return operator()(); }

 public:
	template<class Tuple>
	BOOST_MULTI_HD constexpr auto apply(Tuple const& /*unused*/) const {
		static_assert(std::tuple_size_v<Tuple> == 0);
		return operator()();
	}

	BOOST_MULTI_HD constexpr auto operator&() const& {  // NOLINT(google-runtime-operator)
		return detail::subarray_ptr<T, 0, typename std::pointer_traits<ElementPtr>::template rebind<T const>, Layout, false>(this->base_, this->layout());
	}  // NOLINT(google-runtime-operator) extend semantics  //NOSONAR

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast() const& {
		return const_subarray<T2, 0, P2>(
			typename const_subarray::layout_type{this->layout()},
			reinterpret_pointer_cast<P2>(this->base_)
		);
	}

 private:
	constexpr auto broadcasted() const& {  // NOLINT(readability-identifier-naming) TODO(correaa) remove?
		multi::layout_t<1> const new_layout(this->layout(), 0, 0);
		return subarray<T, 1, typename const_subarray::element_const_ptr>(new_layout, types::base_);
	}

 public:
	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int const /*version*/) const {
		using AT        = multi::archive_traits<Archive>;
		auto&  element_ = *(this->base_);
		arxiv& AT::make_nvp("element", element_);
		//  arxiv & cereal::make_nvp("element", element_);
		//  arxiv &                             element_ ;
	}
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

/// ... a specialization for one dimension
template<typename T, typename ElementPtr, class Layout>
class const_subarray<T, 1, ElementPtr, Layout>  // NOLINT(misc-multiple-inheritance) to define operators via CRTP
: public multi::random_iterable<const_subarray<T, 1, ElementPtr, Layout>>
, public detail::array_types<T, 1, ElementPtr, Layout> {
 public:
	~const_subarray() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	template<class TT, std::enable_if_t<std::is_same_v<ElementPtr, TT const*>, int> = 0>      // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-use-constraints) for C++20
	explicit BOOST_MULTI_HD constexpr const_subarray(std::initializer_list<TT> const& il_1d)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) this constructs a reference to the init list
	: detail::array_types<T, 1, ElementPtr, Layout>(
		  layout_type(multi::extents_t<1>({0, static_cast<size_type>(std::size(il_1d))})),
		  std::data(il_1d)
	  ) {
	}

	// boost serialization needs `delete(...)`. void boost::serialization::extended_type_info_typeid<T>::destroy(const void*) const [with T = boost::multi::subarray<double, 1, double*, boost::multi::layout_t<1> >]
	// void operator delete(void* ptr) noexcept = delete;
	// void operator delete(void* ptr, void* place ) noexcept = delete;  // NOLINT(bugprone-easily-swappable-parameters)

	static constexpr dimensionality_type rank_v = 1;

	using types = detail::array_types<T, dimensionality_type{1}, ElementPtr, Layout>;
	using types::types;

	using rank        = std::integral_constant<dimensionality_type, rank_v>;
	using layout_type = Layout;
	using ref_        = const_subarray;

	using element_type [[deprecated("use ::element")]] = T;
	using element                                      = T;
	using element_ptr                                  = typename types::element_ptr;
	using element_const_ptr                            = typename std::pointer_traits<ElementPtr>::template rebind<element const>;
	using element_move_ptr                             = multi::move_ptr<element, element_ptr>;
	using element_ref                                  = typename types::element_ref;
	using element_cref                                 = typename std::iterator_traits<element_const_ptr>::reference;

	/// `std::allocator_traits<Allocator>::const_pointer` for 1D arrays
	using const_pointer   = element_const_ptr;
	/// `std::allocator_traits<Allocator>::pointer` for 1D arrays
	using pointer         = element_ptr;
	using const_reference = typename detail::array_types<T, dimensionality_type{1}, ElementPtr, Layout>::const_reference;
	using reference       = typename detail::array_types<T, dimensionality_type{1}, ElementPtr, Layout>::reference;

	using default_allocator_type = typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type;

	using typename detail::array_types<T, 1, ElementPtr, Layout>::size_type;

	BOOST_MULTI_HD constexpr auto get_allocator() const -> default_allocator_type { return default_allocator_of(const_subarray::base()); }
	BOOST_MULTI_FRIEND_CONSTEXPR
	auto get_allocator(const_subarray const& self) -> default_allocator_type { return self.get_allocator(); }

	using decay_type = array<std::decay_t<typename types::element>, dimensionality_type{1}, typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type>;

	/// Returns a const-element view of the subarray, preventing modification of elements
	constexpr auto as_const() const { return const_subarray(this->layout(), this->base_); }

	constexpr auto decay() const -> decay_type { return decay_type{*this}; }
	constexpr auto copy() const -> decay_type { return decay_type{*this}; }

	constexpr auto operator+() const { return decay(); }

	using basic_const_array = const_subarray<
		T, 1,
		typename std::pointer_traits<ElementPtr>::template rebind<typename const_subarray::element const>,
		Layout>;

 protected:
	template<class A> constexpr void intersection_assign(A&& other) && { intersection_assign(std::forward<A>(other)); }
	template<class A> constexpr void intersection_assign(A&& other) & {  // NOLINT(cppcoreguidelines-rvalue-reference-param-not-moved,cppcoreguidelines-missing-std-forward) false positive clang-tidy 17
		std::for_each(
			intersection(this->extent(), other.extent()).begin(),
			intersection(this->extent(), other.extent()).end(),
			[&](auto const idx) -> void { operator[](idx) = std::forward<A>(other)[idx]; }
		);
	}

	template<typename, ::boost::multi::dimensionality_type, typename EP, class LLayout> friend class const_subarray;
	template<typename, ::boost::multi::dimensionality_type, class Alloc> friend struct dynamic_array;  // TODO(correaa) check if this is necessary

	template<class T2, class P2, class TT, dimensionality_type DD, class PP>
	friend constexpr auto static_array_cast(subarray<TT, DD, PP> const&) -> decltype(auto);

 public:
	const_subarray(const_subarray const&) = delete;

	friend constexpr auto sizes(const_subarray const& self) noexcept -> typename const_subarray::sizes_type { return self.sizes(); }  // needed by nvcc
	friend constexpr auto size(const_subarray const& self) noexcept -> typename const_subarray::size_type { return self.size(); }     // needed by nvcc

	const_subarray(const_subarray&&) noexcept = default;  // in C++ 14 this was necessary to return array references from functions
	// in c++17 things changed and non-moveable non-copyable types can be returned from functions and captured by auto

 private:
 	using const_ptr = detail::subarray_ptr<T, 1, ElementPtr, Layout, true>;
	using ptr       = detail::subarray_ptr<T, 1, ElementPtr, Layout, false>;

	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;
	template<class, dimensionality_type D, class, bool, bool, typename, class> friend struct detail::array_iterator;

 public:
	friend constexpr auto dimensionality(const_subarray const& /*self*/) -> dimensionality_type { return 1; }

	BOOST_MULTI_HD constexpr auto operator&() const& { return const_ptr{this->base_, this->layout()}; }  // NOLINT(google-runtime-operator) extend semantics  //NOSONAR

	BOOST_MULTI_HD constexpr void assign(std::initializer_list<typename const_subarray::value_type> values) const {
		BOOST_MULTI_ASSERT(values.size() == static_cast<std::size_t>(this->size()));
		if(values.size() != 0) {
			assign(values.begin(), values.end());
		}
	}
	template<class It>
	constexpr auto assign(It first) & -> It {
		adl_copy_n(first, this->size(), this->begin());
		std::advance(first, this->size());
		return first;
	}
	template<class It>
	constexpr auto assign(It first) && -> It { return assign(first); }
	template<class It>
	constexpr void assign(It first, It last) & {
		BOOST_MULTI_ASSERT(std::distance(first, last) == this->size());
		(void)last;
		assign(first);
	}
	template<class It>
	constexpr void assign(It first, It last) && { assign(first, last); }

	// constexpr auto operator=(const_subarray     &&) const& noexcept -> const_subarray const&;  // UNIMPLEMENTABLE! TO PASS THE viewable_range CONCEPT!!!, can't be = delete;
	constexpr auto operator=(const_subarray&&) & noexcept -> const_subarray&;  // UNIMPLEMENTABLE! TO PASS THE viewable_range CONCEPT!!!, can't be = delete;
	constexpr auto operator=(const_subarray const&) const -> const_subarray const& = delete;

	template<
		class ECPtr,
		class = std::enable_if_t<std::is_same_v<element_const_ptr, ECPtr> && !std::is_same_v<element_const_ptr, element_ptr>>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
		>
	constexpr auto operator=(const_subarray<T, 1L, ECPtr, Layout> const& other) const&& -> const_subarray& {
		assert(0);
		operator=(other);
		return *this;
	}  // required by https://en.cppreference.com/w/cpp/iterator/indirectly_writable for std::ranges::copy_n

	// a lightweigh type for multidimensional indexing, it has the indexing interface of an array but without size (extents) information	
	using cursor       = detail::cursor_t<typename const_subarray::element_ptr, 1, typename const_subarray::strides_type>;
	// a lightweigh type for multidimensional indexing (const version), it has the indexing interface of an array but without size (extents) information
	using const_cursor = detail::cursor_t<typename const_subarray::element_const_ptr, 1, typename const_subarray::strides_type>;

	auto diagonal() const = delete;

	auto flattened() const& = delete;  // { return flattened_aux_().as_const(); }

 private:
	BOOST_MULTI_HD constexpr auto home_aux_() const { return cursor(this->base_, this->strides()); }

 public:
	BOOST_MULTI_HD constexpr auto home() const& -> const_cursor { return home_aux_(); }

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	BOOST_MULTI_HD constexpr auto at_aux_(index idx) const -> typename const_subarray::reference {  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment
		// NOLINTNEXTLINE(readability-static-accessed-through-instance) can be static
		if constexpr(std::is_integral_v<decltype(this->stride())>) {
			BOOST_MULTI_ASSERT((this->stride() == 0 || (this->extent().contains(idx))) && ("out of bounds"));
		}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
		return *(this->base_ + (this->stride() * idx - this->offset()));  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,readability-static-accessed-through-instance) can be static

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif
	}

	constexpr auto broadcasted() const& {  // NOLINT(readability-identifier-naming) TODO(correaa) remove?
		// multi::layout_t<1> const self_layout{this->layout()};
		// TODO(correaa) introduce a broadcasted_layout?
		multi::layout_t<2> const new_layout(this->layout(), 0, 0, 1);  // , (std::numeric_limits<size_type>::max)()};
		return const_subarray<T, 2, ElementPtr, multi::layout_t<2>>(new_layout, types::base_);
	}

 public:
	constexpr auto repeated(size_type n) && {
		auto exts = this->extents();  // mull-ignore: cxx_init_const

		return [self = std::move(*this)](auto /*idx*/, auto... rest) -> decltype(auto) { return detail::invoke_square(self, rest...); } ^ /*(*/ n* exts /*)*/;
	}

	constexpr auto repeated(size_type n) const& {
		return [this](auto /*idx*/, auto... rest) -> decltype(auto) { return detail::invoke_square(*this, rest...); } ^ /*(*/ n* this->extents() /*)*/;
	}

	// template<template<class...> class Container = std::vector, class... As>
	// constexpr auto to(As const&... args) const& {
	// 	using inner_value_type = typename const_subarray::value_type;
	// 	using container_type   = Container<inner_value_type, As...>;
	// 	return container_type(this->begin(), this->end(), args...);
	// }

	BOOST_MULTI_HD constexpr auto operator[](index idx) const& -> typename const_subarray::const_reference { return at_aux_(idx); }  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment

	BOOST_MULTI_HD constexpr auto front() const& -> const_reference { return *begin(); }
	BOOST_MULTI_HD constexpr auto back() const& -> const_reference { return *std::prev(end(), 1); }

 private:
	template<class Self, typename Tuple, std::size_t... I, const_subarray* = nullptr>
	static constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...> /*012*/) -> decltype(auto) {
		using std::get;  // for C++17 compatibility
		return std::forward<Self>(self)(get<I>(tuple)...);
	}

 public:
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tuple) const& -> decltype(auto) { return apply_impl_(*this, tuple, std::make_index_sequence<std::tuple_size_v<Tuple>>()); }

	// template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value == 0), int> = 0> BOOST_MULTI_HD constexpr auto operator[](Tuple const& /*empty*/) const& -> decltype(auto) { return *this; }  // NOLINT(modernize-use-constraints) for C++20
	// template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value == 1), int> = 0> BOOST_MULTI_HD constexpr auto operator[](Tuple const& indices) const& -> decltype(auto) {                    // NOLINT(modernize-use-constraints) for C++20
	// 	using std::get;
	// 	return operator[](get<0>(indices));
	// }

	// template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value > 1), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// BOOST_MULTI_HD constexpr auto operator[](Tuple const& indices) const& -> decltype(operator[](std::get<0>(indices))[detail::tuple_tail(indices)]) {
	// 	using std::get;  // for C++17 compatibility
	// 	return operator[](get<0>(indices))[detail::tuple_tail(indices)];
	// }

// Warning C4459 comes from boost::multi_array having a namespace indices which collides with the variable name?
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4459)
#endif

	[[deprecated("BMA compat, finish impl")]] BOOST_MULTI_HD constexpr auto operator[](std::tuple<irange> const& indices) const& {
		using std::get;
		return (*this)({get<0>(indices).front(), get<0>(indices).back() + 1});
	}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

	BOOST_MULTI_HD constexpr auto elements_at(size_type idx) const& -> decltype(auto) {
		BOOST_MULTI_ASSERT(idx < this->num_elements());
		return operator[](idx);
	}
	BOOST_MULTI_HD constexpr auto elements_at(size_type idx) && -> decltype(auto) {
		BOOST_MULTI_ASSERT(idx < this->num_elements());
		return operator[](idx);
	}
	BOOST_MULTI_HD constexpr auto elements_at(size_type idx) & -> decltype(auto) {
		BOOST_MULTI_ASSERT(idx < this->num_elements());
		return operator[](idx);
	}

	constexpr auto reindexed(index first) && { return reindexed(first); }
	constexpr auto reindexed(index first) & { return const_subarray{this->layout().reindex(first), types::base_}; }

 private:
	BOOST_MULTI_HD constexpr auto taked_aux_(difference_type count) const {
		BOOST_MULTI_ASSERT(count <= this->size());  // calculating size is expensive that is why
		typename types::layout_type const new_layout(
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride() * count
		);
		return const_subarray{new_layout, this->base_};
	}

 public:
	constexpr auto taked(difference_type count) const& -> const_subarray<T, 1, ElementPtr, Layout> { return taked_aux_(count); }

 private:
	BOOST_MULTI_HD constexpr auto dropped_aux_(difference_type count) const {

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
		return const_subarray(
			this->layout().drop(count), this->base_ + /*(*/ count * this->layout().stride() /*- this->layout().offset()*/ /*)*/  // TODO(correaa) fix need for offset  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		);
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif
	}

 public:
	constexpr auto dropped(difference_type count) const& -> const_subarray { return dropped_aux_(count); }

 private:
	BOOST_MULTI_HD constexpr auto sliced_aux_(index first, index last) const {

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return const_subarray{this->layout().slice(first, last), this->base_ + /*(*/ first * this->layout().stride() /*- this->layout().offset())*/};  // TODO(correaa) fix need for offset

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif
	}

 public:
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) const& -> basic_const_array { return basic_const_array{sliced_aux_(first, last)}; }
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) & -> const_subarray { return sliced_aux_(first, last); }
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) && -> const_subarray { return sliced_aux_(first, last); }

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
 private:
	auto to_mdspan_aux_() const {
		using std::apply;
		auto shape = apply(
			[](auto... sizes) -> auto { return std::dextents<std::size_t, 1>{static_cast<std::size_t>(sizes)...}; },
			this->sizes()
		);

		auto strides = apply(
			[](auto... strds) -> auto { return std::array<std::size_t, 1>{static_cast<std::size_t>(strds)...}; },
			this->strides()
		);

		return std::mdspan<T, std::dextents<std::size_t, 1>, std::layout_stride>{
			this->base_, std::layout_stride::mapping(shape, strides)
		};
	}
	auto to_mdspan_() const& {
		return std::mdspan<T const, std::dextents<std::size_t, 1>, std::layout_stride>(to_mdspan_aux_());
	}

 public:
	/// A (const_)subarray can be converted implicitly to a strided mdspan object
	operator std::mdspan<T const, std::dextents<std::size_t, 1>, std::layout_stride>() const& { return to_mdspan_(); }
#endif

 private:
	using elements_iterator  = detail::elements_iterator_t<element_ptr, layout_type>;
	using celements_iterator = detail::elements_iterator_t<element_const_ptr, layout_type>;

	using elements_range       = detail::elements_range_t<element_ptr, layout_type>;
	using const_elements_range = detail::elements_range_t<element_const_ptr, layout_type>;

	constexpr auto elements_aux_() const { return elements_range{this->base_, this->layout()}; }

 public:
	constexpr auto elements() & -> elements_range { return elements_aux_(); }
	constexpr auto elements() && -> elements_range { return elements_aux_(); }
	constexpr auto elements() const& -> const_elements_range { return const_elements_range{this->base(), this->layout()}; }  // TODO(correaa) simplify

	constexpr auto celements() const -> const_elements_range { return elements_aux_(); }

	constexpr auto blocked(index first, index last) & -> const_subarray {
		return sliced(first, last).reindexed(first);
	}
	constexpr auto stenciled(typename const_subarray::index_extension ext) -> const_subarray {
		return blocked(ext.first(), ext.last());
	}

 private:
	constexpr auto strided_aux_(difference_type step) const {
		assert(this->size() % step == 0);
		auto const new_layout = typename types::layout_type{this->layout().sub(), this->layout().stride() * step, this->layout().offset(), this->layout().nelems()};
		return subarray<T, 1, ElementPtr, Layout>(new_layout, types::base_);
	}

 public:
	constexpr auto strided(difference_type step) const& -> const_subarray { return strided_aux_(step); }

	BOOST_MULTI_HD constexpr auto sliced(index first, index last, difference_type stride) const& -> basic_const_array { return sliced(first, last).strided(stride); }

 private:
	// NOLINTNEXTLINE(readability-identifier-naming)
	BOOST_MULTI_HD constexpr auto range(index_range const& rng) const& { return sliced(rng.front(), rng.last()); }

	BOOST_MULTI_HD constexpr auto paren_aux_() const& { return const_subarray(this->layout(), this->base_); }

	BOOST_MULTI_HD constexpr auto paren_aux_(index idx) const& -> decltype(auto) { return operator[](idx); }

	BOOST_MULTI_HD constexpr auto paren_aux_(index_range const& rng) const& { return range(rng); }

 public:
	BOOST_MULTI_HD constexpr auto operator()() const& { return paren_aux_(); }
#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	BOOST_MULTI_HD constexpr auto operator[]() const& -> const_subarray { return paren_aux_(); }
#endif

	/// yields a subarray that is one dimension lower at index `idx`
	BOOST_MULTI_HD constexpr auto operator()(index idx) const -> decltype(auto) { return operator[](idx); }

	/// Subarray spanning the given index range `rng` along the outermost dimension
	BOOST_MULTI_HD constexpr auto operator()(index_range const& rng) const& { return range(rng); }

 private:
	constexpr auto paren_aux_(intersecting_range<index> const& rng) const& -> decltype(auto) { return paren_aux_(intersection(this->extent(), rng)); }

 public:
	BOOST_MULTI_HD constexpr auto operator()(intersecting_range<index> const& isrange) const& -> decltype(auto) { return paren_aux_(isrange); }

	template<class... Args>
	BOOST_MULTI_HD constexpr auto operator()(Args&&... args) const& -> decltype(paren_(std::declval<const_subarray const&>(), std::forward<Args>(args)...)) {
		return paren_(*this, std::forward<Args>(args)...);
	}

 private:
	BOOST_MULTI_HD constexpr auto halved_aux_() const {
		auto new_layout = this->layout().halve();
		return subarray<T, 2, element_ptr>(new_layout, this->base_);
	}

 public:
	BOOST_MULTI_HD constexpr auto halved() const& -> const_subarray<T, 2, element_ptr> { return halved_aux_(); }

 private:
	BOOST_MULTI_HD constexpr auto partitioned_aux_(size_type size) const {
		BOOST_MULTI_ASSERT(size != 0);
		BOOST_MULTI_ASSERT((this->layout().nelems() % size) == 0);  // TODO(correaa) remove assert? truncate left over? (like mathematica)
		multi::layout_t<2> new_layout{this->layout(), this->layout().nelems() / size, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= size;  // TODO(correaa) : don't use mutation
		return subarray<T, 2, element_ptr>(new_layout, types::base_);
	}

 public:
	/// Yields a view of higher dimension by splitting the leading dimension into `size` equal-sized partitions (`size() % size` must be zero)
	BOOST_MULTI_HD constexpr auto partitioned(size_type size) const& -> const_subarray<T, 2, element_ptr> { return partitioned_aux_(size); }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
 private:
	BOOST_MULTI_HD constexpr auto splitted_aux_() const {
		multi::layout_t<1> const lyt1({}, this->layout().stride(), 0, this->layout().nelems() / this->layout().stride() / 2 * this->layout().stride());
		multi::layout_t<1> const lyt2({}, this->layout().stride(), 0, ((this->layout().nelems() / this->layout().stride()) + 1) / 2 * this->layout().stride());
		return  // std::array<subarray<T, 1, element_ptr>, 2>
			std::pair<subarray<T, 1, element_ptr>, subarray<T, 1, element_ptr>>(
				subarray<T, 1, element_ptr>(lyt1, types::base_),
				subarray<T, 1, element_ptr>(lyt2, types::base_ + lyt1.nelems())  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
			);
	}
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	BOOST_MULTI_HD constexpr auto splitted() const {
		return std::pair<const_subarray<T, 1, element_ptr>, const_subarray<T, 1, element_ptr>>(
			std::get<0>(splitted_aux_()),
			std::get<1>(splitted_aux_())
		);
	}

 private:
	BOOST_MULTI_HD constexpr auto chunked_aux_(size_type size) const {
		BOOST_MULTI_ASSERT(this->size() % size == 0);
		return partitioned_aux_(this->size() / size);
	}

 public:  // in Mathematica this is called Partition https://reference.wolfram.com/language/ref/Partition.html in RangesV3 it is called chunk
	BOOST_MULTI_HD constexpr auto chunked(size_type size) const& -> const_subarray<T, 2, element_ptr> { return chunked_aux_(size); }

	constexpr auto tiled(size_type count) const& {
		BOOST_MULTI_ASSERT(count != 0);

		struct divided_type {
			const_subarray<T, 2, element_ptr> quotient;
			const_subarray<T, 1, element_ptr> remainder;
		};

		return divided_type{
			this->taked(this->size() - (this->size() % count)).chunked(count),
			this->dropped(this->size() - (this->size() % count)),
		};
	}

 private:
	constexpr auto reversed_aux_() const -> const_subarray {
		auto new_layout = this->layout().reverse();
		return const_subarray(new_layout, types::base_);
	}

 public:
	constexpr auto reversed() const& { return reversed_aux_().as_const(); }
	constexpr auto reversed() & { return reversed_aux_(); }
	constexpr auto reversed() && { return reversed_aux_(); }

	BOOST_MULTI_HD constexpr auto rotated() const& { return operator()(); }    // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.9
	BOOST_MULTI_HD constexpr auto unrotated() const& { return operator()(); }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.9
	BOOST_MULTI_HD constexpr auto unordered() const& { return operator()(); }

	auto transposed() const& = delete;
	auto flatted() const&    = delete;

	using iterator       = typename multi::detail::array_iterator<element, 1, typename types::element_ptr, false, false, typename layout_type::stride_type>;
	using const_iterator = typename multi::detail::array_iterator<element, 1, typename types::element_ptr, true, false, typename layout_type::stride_type>;

	// /// Iterator in the leading dimension that mark elements as movable
	// using move_iterator  = typename multi::detail::array_iterator<element, 1, typename types::element_ptr, false, true>;

	using reverse_iterator [[deprecated]]       = std::reverse_iterator<iterator>;
	using const_reverse_iterator [[deprecated]] = std::reverse_iterator<const_iterator>;

	struct [[deprecated("BMA compatibility")]] index_gen {
		auto operator[](irange const& rng) const { return std::make_tuple(rng); }
	};
	using extent_gen [[deprecated("BMA compatibility")]]   = std::array<irange, 1>;
	using extent_range [[deprecated("BMA compatibility")]] = irange;

	template<
		class Range,
		std::enable_if_t<!has_extents<std::decay_t<Range>>::value, int>         = 0,
		std::enable_if_t<!detail::is_subarray<std::decay_t<Range>>::value, int> = 0,
		class                                                                   = decltype((void)std::declval<Range>().begin(), std::declval<Range>().end()),
		class                                                                   = decltype(Range{std::declval<typename const_subarray::const_iterator>(), std::declval<typename const_subarray::const_iterator>()})>
	constexpr explicit operator Range() const {
		// vvv Range{...} needed by Windows GCC?
		return Range{begin(), end()};  // e.g. std::vector(it, it, alloc = {})
	}

 private:
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	BOOST_MULTI_HD constexpr auto begin_aux_() const { return iterator{this->base_, this->layout().sub(), this->stride()}; }
	BOOST_MULTI_HD constexpr auto end_aux_() const { return iterator{this->base_ + types::nelems(), this->layout().sub(), this->stride()}; }  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

 public:
	BOOST_MULTI_HD constexpr auto begin() const& -> const_iterator { return begin_aux_(); }  ///< returns a iterator to the beginning
	BOOST_MULTI_HD constexpr auto end() const& -> const_iterator { return end_aux_(); }      ///< returns a iterator to the end

#if defined(__GNUC__) && !defined(__EDG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	[[deprecated("implement as negative stride")]] constexpr auto rbegin() const& { return const_reverse_iterator(end()); }  // TODO(correaa) implement as negative stride?
	[[deprecated("implement as negative stride")]] constexpr auto rend() const& { return const_reverse_iterator(begin()); }  // TODO(correaa) implement as negative stride?
#if defined(__GNUC__) && !defined(__EDG__)
#pragma GCC diagnostic pop
#endif

	/// returns an (explicitly const-)iterator to the beginning
	BOOST_MULTI_HD constexpr auto cbegin() const& -> const_iterator { return begin(); }
	/// returns an (explicitly const-)iterator to the end
	BOOST_MULTI_HD constexpr auto cend() const& -> const_iterator { return end(); }

	template<class It> constexpr auto assign(It first) && -> decltype(adl_copy_n(first, std::declval<size_type>(), std::declval<iterator>()), void()) {
		return adl_copy_n(first, this->size(), std::move(*this).begin()), void();
	}

	friend constexpr auto operator==(const_subarray const& self, const_subarray const& other) -> bool {
		return self.extent() == other.extent() && self.elements() == other.elements();
	}

	friend constexpr auto operator!=(const_subarray const& self, const_subarray const& other) -> bool {
		return self.extent() != other.extent() || self.elements() != other.elements();
	}

	template<class OtherT, typename OtherEP, class OtherLayout>
	friend constexpr auto operator==(const_subarray const& self, const_subarray<OtherT, 1, OtherEP, OtherLayout> const& other) -> bool {
		return self.extent() == other.extent() && self.elements() == other.elements();
	}

	template<class TT, typename EEPP, class LL>
	friend constexpr auto operator!=(const_subarray const& self, const_subarray<TT, 1, EEPP, LL> const& other) -> bool {
		return self.extent() != other.extent() || self.elements() != other.elements();
	}

	friend constexpr auto operator<(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(self, other); }
	friend constexpr auto operator>(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(other, self); }  // NOLINT(readability-suspicious-call-argument)

	friend constexpr auto operator<=(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(self, other) || self == other; }
	friend constexpr auto operator>=(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(other, self) || self == other; }  // NOLINT(readability-suspicious-call-argument)

	template<class Range, typename = std::enable_if_t<!detail::is_const_subarray_v<Range>>, typename = decltype(void(std::declval<Range const&>().extents()), std::declval<Range const&>().elements())>
	friend constexpr auto operator==(const_subarray const& self, Range const& other) -> bool {
		return self.extents() == other.extents() && self.elements() == other.elements();
	}

	template<class Range, typename = std::enable_if_t<!detail::is_const_subarray_v<Range>>, typename = decltype(void(std::declval<Range const&>().extents()), std::declval<Range const&>().elements())>
	friend constexpr auto operator==(Range const& other, const_subarray const& self) -> bool {
		return self.extents() == other.extents() && self.elements() == other.elements();
	}

	template<class Range, typename = std::enable_if_t<!detail::is_const_subarray_v<Range>>, typename = decltype(void(std::declval<Range const&>().extents()), std::declval<Range const&>().elements())>
	friend constexpr auto operator!=(const_subarray const& self, Range const& other) -> bool {
		return self.extents() != other.extents() || self.elements() != other.elements();
	}

 private:
	template<class A1, class A2>
	static constexpr auto lexicographical_compare_(A1 const& self, A2 const& other) -> bool {  // NOLINT(readability-suspicious-call-argument)
		if(self.extent().first() > other.extent().first()) {
			return true;
		}
		if(self.extent().first() < other.extent().first()) {
			return false;
		}
		return adl_lexicographical_compare(adl_begin(self), adl_end(self), adl_begin(other), adl_end(other));
	}

 public:
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() const -> subarray<T2, 1, P2, Layout> {  // name taken from std::static_pointer_cast
		return {this->layout(), static_cast<P2>(this->base_)};
	}
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>, class... Args>
	constexpr auto static_array_cast(Args&&... args) const -> subarray<T2, 1, P2, Layout> {  // name taken from std::static_pointer_cast
		return subarray<T2, 1, P2, Layout>(
			this->layout(), P2{this->base_, std::forward<Args>(args)...}
		);
	}

	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) const& {
		return static_array_cast<
			//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
			std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
			transform_ptr<
				//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
				std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
				UF, element_const_ptr, std::invoke_result_t<UF const&, element_cref>>>(std::forward<UF>(fun));
	}
	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) & {
		return static_array_cast<
			std::decay_t<std::invoke_result_t<UF const&, element_ref>>,
			transform_ptr<
				std::decay_t<std::invoke_result_t<UF const&, element_ref>>,
				UF, element_ptr, std::invoke_result_t<UF const&, element_ref>>>(std::forward<UF>(fun));
	}
	template<class UF>
	BOOST_MULTI_HD constexpr auto element_transformed(UF&& fun) && { return element_transformed(std::forward<UF>(fun)); }

	template<
		class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM      = T2 std::decay_t<Element>::*>
	constexpr auto member_cast(PM member) const {
		static_assert(sizeof(T) % sizeof(T2) == 0, "array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
												   "Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) reinterpret is what the function does. alternative for GCC/NVCC
		auto&& ref1 = (*(reinterpret_cast<typename const_subarray::element* const&>(const_subarray::base_))).*member;  // ->*pm;
		auto*  ptr1 = &ref1;                                                                                           //-V::537 ptr1 is reinterpreted (not dereferenced) below to support fancy pointer types

		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa) find a better way
		P2 ptr2 = reinterpret_cast<P2&>(ptr1);  // NOSONAR
#else
		auto ptr2 = static_cast<P2>(&(this->base_->*member));  // this crashes nvcc 11.2-11.4 and some? gcc compiler
#endif
		return subarray<T2, 1, P2>(this->layout().scale(sizeof(T), sizeof(T2)), ptr2);
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast() const& {
		BOOST_MULTI_ASSERT(this->layout().stride() * static_cast<size_type>(sizeof(T)) % static_cast<size_type>(sizeof(T2)) == 0);

		return const_subarray<T2, 1, P2>(
			layout_type{this->layout().sub(), this->layout().stride() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2)), this->layout().offset() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2)), this->layout().nelems() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2))},
			reinterpret_pointer_cast<P2>(this->base_)
		);
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2 const>>
	constexpr auto reinterpret_array_cast(size_type n) const& -> subarray<std::decay_t<T2>, 2, P2> {  // TODO(correaa) : use rebind for return type
		static_assert(sizeof(T) % sizeof(T2) == 0, "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

		return subarray<std::decay_t<T2>, 2, P2>(
				   layout_t<2>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, n),
				   reinterpret_pointer_cast<P2>(this->base())
		)
			.rotated();
	}

	template<class Archive>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		using AT = multi::archive_traits<Archive>;
		std::for_each(this->begin(), this->end(), [&](reference& item) -> void { arxiv& AT ::make_nvp("item", item); });
		//  std::for_each(this->begin(), this->end(), [&](auto&&     item) {arxiv & cereal::make_nvp("item", item);});
		//  std::for_each(this->begin(), this->end(), [&](auto&&     item) {arxiv &                          item ;});
	}
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class T2, class P2, class Array, class... Args>
constexpr auto static_array_cast(Array&& self, Args&&... args) -> decltype(auto) {
	return std::forward<Array>(self).template static_array_cast<T2, P2>(std::forward<Args>(args)...);
}

namespace detail {
template<class T, dimensionality_type D, typename Ptr = T*>
struct array_ptr;

template<class T, dimensionality_type D, typename Ptr = T*>
using array_cptr = array_ptr<T, D, typename std::pointer_traits<Ptr>::template rebind<T const>>;

}  // namespace detail

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

/// A `D`-dimensional view of a contiguous, pre-existing memory buffer.
///
/// Does not own or manage the elements it references.
/// Has reference semantics: cannot be rebound after construction, assignments are deep,
/// and the size is fixed for the object's lifetime.
/// Invalidated if the underlying buffer is deallocated or moved.
///
/// @tparam T Element type
/// @tparam D Dimensionality (non-negative)
/// @tparam ElementPtr Pointer-like type to the elements (default `T*`)
/// @tparam Layout Layout type describing strides and extensions
template<
	typename T, dimensionality_type D, typename ElementPtr = T*,
	class Layout =
		std::conditional_t<
			(D == 1),
			// contiguous_layout<>,  // 1, typename std::pointer_traits<ElementPtr>::difference_type>,
			multi::layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>,
			multi::layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>>>
class array_ref : public subarray<T, D, ElementPtr, Layout> {
	using subarray_layout = Layout;

	using subarray_base = subarray<T, D, ElementPtr, Layout>;

 public:
	~array_ref() = default;  // lints(cppcoreguidelines-special-member-functions)

	/// Type to describe the layout of an array, that results from `.layout()` member.
	using layout_type = typename subarray_base::layout_type;

	/// Type for random-access iteration in the leading dimension, that results from `.begin()`/`.end()` members
	using iterator = typename subarray_base::iterator;

	using typename subarray_base::size_type;

	constexpr array_ref() = delete;  // because reference cannot be unbound

	// [[deprecated("references are not copyable, use auto&&")]]
	array_ref(array_ref const&) = delete;   // don't try to use `auto` for references, use `auto&&` or explicit value type
	array_ref(array_ref&&)      = default;  // movable (shallow handle move) so a temporary models std::ranges::view, e.g. for std::views::zip; copy stays deleted, so `auto x = ref;` is still an error

	array_ref(iterator, iterator) = delete;

	// return type removed for MSVC
	friend constexpr auto sizes(array_ref const& self) noexcept /*-> typename array_ref::sizes_type*/ { return self.sizes(); }  // needed by nvcc
	friend constexpr auto size(array_ref const& self) noexcept /*-> typename array_ref::size_type*/ { return self.size(); }     // needed by nvcc

	constexpr auto flatted() const& {
		assert(this->layout().is_flattable());
		multi::layout_t<D - 1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return const_subarray<T, D - 1, ElementPtr>{new_layout, this->base_};
	}

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto flatted() & {
		assert(this->layout().is_flattable());
		multi::layout_t<D - 1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return subarray<T, D - 1, ElementPtr>(new_layout, this->base_);
	}

	constexpr auto flatted() && { return this->flatted(); }  // cppcheck-suppress duplInheritedMember ; to override

#if defined(BOOST_MULTI_HAS_SPAN) && !defined(__NVCC__)
	/// conversion to `std::span` (`D == 1` only)
	template<class U = typename array_ref::element, std::enable_if_t<std::is_convertible_v<typename array_ref::element_const_ptr, U const*> && (D == 1), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr explicit operator std::span<U const>() const { return std::span<U const>(this->data_elements(), this->size()); }

	template<class U = typename array_ref::element, std::enable_if_t<std::is_convertible_v<typename array_ref::element_ptr, U*> && (D == 1), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr explicit operator std::span<U>() { return std::span<U>(this->data_elements(), this->size()); }
#endif

	template<class OtherPtr, class = std::enable_if_t<!std::is_same_v<OtherPtr, ElementPtr>>, decltype(multi::detail::explicit_cast<ElementPtr>(std::declval<OtherPtr>()))* = nullptr>
	constexpr explicit array_ref(array_ref<T, D, OtherPtr>&& other)
	: subarray_base(other.layout(), ElementPtr{std::move(other).base()}) {}  // cppcheck-suppress internalAstError ; bug in cppcheck 2.13.0

	template<class OtherPtr, class = std::enable_if_t<!std::is_same_v<OtherPtr, ElementPtr>>, decltype(multi::detail::implicit_cast<ElementPtr>(std::declval<OtherPtr>()))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	constexpr /*mplct*/ array_ref(array_ref<T, D, OtherPtr>&& other)         // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor,bugprone-use-after-move,hicpp-invalid-access-moved)
	: subarray_base(other.layout(), ElementPtr{std::move(other).base()}) {}  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)

	constexpr array_ref(ElementPtr dat, ::boost::multi::extents_t<D> const& exts) noexcept  // TODO(correa) eliminate this ctor
	: subarray_base(typename subarray_base::layout_type(exts), dat) {}

	/// constructs a D-dimensional view of the contiguous range starting at p and ending at least after the size of the multidimensional array (product of sizes).
	explicit constexpr array_ref(::boost::multi::extents_t<D> exts, ElementPtr dat) noexcept
	: subarray_base{typename array_ref::layout_type(exts), dat} {}

#if defined(BOOST_MULTI_HAS_SPAN) && !defined(__NVCC__)
#ifdef __cpp_lib_span
	template<class Dummy = void*, std::enable_if_t<sizeof(Dummy) && (D == 1), int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	// explicit converts a more primitive type into a more powerful type (also removing explicit generates a problem with nvc++ 26)
	explicit constexpr array_ref(std::span<typename array_ref::element>&& data_ref)
	: array_ref({static_cast<typename array_ref::size_type>(data_ref.size())}, data_ref.data()) { (void)std::move(data_ref); }
#endif
#endif

	// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)  // compatibility with legacy c-arrays
	template<
		class Array,
		std::enable_if_t<  // NOLINT(modernize-use-constraints) for C++20
			!std::is_array_v<Array> && !std::is_base_of_v<array_ref, std::decay_t<Array>> && std::is_convertible_v<decltype(multi::data_elements(std::declval<Array&>())), ElementPtr>, int> = 0>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax and because a reference to c-array can be represented as an array_ref
	constexpr array_ref(Array& array)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax and because a reference to c-array can be represented as an array_ref
	: array_ref(multi::data_elements(array), extensions(array)) {}
	// NOLINTEND(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

	template<class TT = void, std::enable_if_t<sizeof(TT*) && D == 0, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax and because a reference to c-array can be represented as an array_ref
	constexpr array_ref(  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax and because a reference to c-array can be represented as an array_ref
		T& elem           // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	)
	: array_ref(&elem, {}) {}

	template<class TT, std::size_t N>
	// cppcheck-suppress noExplicitConstructor ; see below
	constexpr array_ref(TT (&arr)[N])  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : for backward compatibility // NOSONAR
	: array_ref(
		  ::boost::multi::extents(arr),
		  ::boost::multi::data_elements(arr)
	  ) {}

	template<class TT, std::size_t N>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr array_ref(std::array<TT, N>& arr) : array_ref(::boost::multi::extents(arr), ::boost::multi::data_elements(arr)) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) array_ptr is more general than pointer c-array support legacy c-arrays  // NOSONAR

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) bug in clang-tidy 19?
	template<class TT, std::enable_if_t<std::is_same_v<typename array_ref::value_type, TT>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// cppcheck-suppress noExplicitConstructor
	explicit array_ref(std::initializer_list<TT> il_1d)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: array_ref(
		  (il_1d.size() == 0) ? nullptr
							  : il_1d.begin(),  // TODO(correaa) simplify conditional by still using a il pointer in empty case?
		  typename array_ref::extents_type{static_cast<typename array_ref::size_type>(il_1d.size())}
	  ) {}

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) bug in clang-tidy 19?
	template<class TT, std::enable_if_t<std::is_same_v<typename array_ref::value_type, TT>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// cppcheck-suppress noExplicitConstructor
	explicit array_ref(std::initializer_list<TT>&& il_1d) = delete;

	using subarray_base::operator=;

 private:
	constexpr auto addressof_aux_() const { return detail::array_ptr<T, D, ElementPtr>(this->base_, this->extents()); }

	constexpr auto addressof_() && -> detail::array_ptr<T, D, ElementPtr> { return addressof_aux_(); }       // cppcheck-suppress duplInheritedMember;
	constexpr auto addressof_() & -> detail::array_ptr<T, D, ElementPtr> { return addressof_aux_(); }        // cppcheck-suppress duplInheritedMember;
	constexpr auto addressof_() const& -> detail::array_cptr<T, D, ElementPtr> { return addressof_aux_(); }  // cppcheck-suppress duplInheritedMember;

 public:
	// operator& is not defined for r-values anyway
	// NOLINTNEXTLINE(google-runtime-operator)
	constexpr auto operator&() && { return addressof_(); }  // cppcheck-suppress duplInheritedMember;  // NOLINT(runtime/operator)  //NOSONAR
	// [[deprecated("controversial")]]
	// constexpr auto operator&() & { return addressof_(); }  // NOLINT(runtime/operator) //NOSONAR
	// // [[deprecated("controversial")]]
	// constexpr auto operator&() const& { return addressof_(); }  // NOLINT(runtime/operator) //NOSONAR

 private:
	template<class It> constexpr auto copy_elements_(It first) {
		return adl_copy_n(first, this->num_elements(), this->data_elements());
	}

 public:
	/// pointer to a contiguous range of `.num_elements()` that contains the elements of the array
	BOOST_MULTI_HD constexpr auto data_elements() const& { return static_cast<typename array_ref::element_const_ptr>(array_ref::base_); }

	template<class TT, class... As, std::enable_if_t<!std::is_base_of_v<array_ref, array_ref<TT, D, As...>>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator=(array_ref<TT, D, As...> const& other) && -> array_ref& {                                   // if MSVC complains here, it probably needs /EHsc /permissive- for C++17 mode
		BOOST_MULTI_ASSERT(this->extents() == other.extents());
		array_ref::copy_elements_(other.data_elements());
		return *this;
	}

	constexpr auto operator=(array_ref const& other) & -> array_ref& {
		if(this == std::addressof(other)) {
			return *this;
		}  // lints(cert-oop54-cpp)
		// TODO(correaa) assert on extensions, not on num elements
		BOOST_MULTI_ASSERT(this->num_elements() == other.num_elements());
		array_ref::copy_elements_(other.data_elements());
		return *this;
	}

	constexpr auto operator=(array_ref const& other) && -> array_ref& {
		if(this == std::addressof(other)) {
			return *this;
		}  // lints(cert-oop54-cpp)
		operator=(other);
		return *this;
	}

	constexpr auto operator=(array_ref&& other) & noexcept(std::is_nothrow_copy_assignable_v<T>)  // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor,cppcoreguidelines-noexcept-move-operations)  //NOSONAR(cppS5018)
		-> array_ref& {
		if(this == std::addressof(other)) {
			return *this;
		}  // lints(cert-oop54-cpp)
		operator=(std::as_const(other));
		return *this;
	}

	constexpr auto operator=(array_ref&& other) && noexcept(std::is_nothrow_copy_assignable_v<T>)  // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor,cppcoreguidelines-noexcept-move-operations)
		-> array_ref& {
		if(this == std::addressof(other)) {
			return *this;
		}  // lints(cert-oop54-cpp)
		operator=(std::as_const(other));
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
	auto operator=(array_ref<TT, DD, As...> const& other) & -> array_ref& {
		BOOST_MULTI_ASSERT(this->extents() == other.extents());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr auto operator=(array_ref<TT, DD, As...> const& other) && -> array_ref& {
		this->operator=(other);
		return *this;  // lints (cppcoreguidelines-c-copy-assignment-signature)
	}

 private:
	using elements_type  = array_ref<typename array_ref::element, 1, typename array_ref::element_ptr>;
	using celements_type = array_ref<typename array_ref::element, 1, typename array_ref::element_const_ptr>;

	constexpr auto elements_aux_() const {
		return elements_type(
			this->base_,
			static_cast<typename elements_type::extents_type>(multi::iextension(this->num_elements()))
		);
	}

 public:
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	/// yields a random-access range with all the elements of the array (no copies or allocations are made, O(1) operation)
	constexpr auto elements() const& noexcept -> celements_type { return elements_aux_(); }
	constexpr auto elements() & noexcept -> elements_type { return elements_aux_(); }
	constexpr auto elements() && noexcept -> elements_type { return elements_aux_(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	// friend constexpr auto elements(array_ref& self) -> elements_type { return self.elements(); }
	// friend constexpr auto elements(array_ref&& self) -> elements_type { return std::move(self).elements(); }
	// friend constexpr auto elements(array_ref const& self) -> celements_type { return self.elements(); }

 private:
	constexpr auto celements_() const& { return celements_type{array_ref::data_elements(), array_ref::num_elements()}; }  // cppcheck-suppress duplInheritedMember ; to overwrite

 public:
	// // cppcheck-suppress-begin duplInheritedMember ; to overwrite
	// constexpr auto element_moved() & { return array_ref<T, D, typename array_ref::element_move_ptr, Layout>(this->extents(), typename array_ref::element_move_ptr{this->base_}); }
	// constexpr auto element_moved() && { return element_moved(); }
	// // cppcheck-suppress-end duplInheritedMember ; to overwrite

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	template<typename TT, class... As>
	friend constexpr auto operator==(array_ref const& self, array_ref<TT, D, As...> const& other) -> bool {
		if(self.extents() != other.extents()) {
			return false;
		}
		return adl_equal(
			other.data_elements(), other.data_elements() + other.num_elements(),  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) TODO(correaa) use span?
			self.data_elements()
		);
	}

#ifdef _MSC_VER

	// Workaround for a standard library bug in MSVC 14.3 and greater
	/*
	compile-c-c++ ..\..\..\bin.v2\libs\boost-multi\test\allocator.test\msvc-14.3\debug\cxxstd-20-iso\threading-multi\allocator.obj
	allocator.cpp
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): error C3889: call to object of class type 'std::equal_to<void>': no matching call operator found
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(742): note: could be 'unknown-type std::equal_to<void>::operator ()(_Ty1 &&,_Ty2 &&) noexcept(<expr>) const'
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): note: Failed to specialize function template 'unknown-type std::equal_to<void>::operator ()(_Ty1 &&,_Ty2 &&) noexcept(<expr>) const'
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): note: With the following template arguments:
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): note: '_Ty1=const _Ty &'
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): note: '_Ty2=const _Ty &'
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5646): note: the template instantiation context (the oldest one first) is
	allocator.cpp(141): note: see reference to function template instantiation 'bool std::operator ==<boost::multi::array<int,2,std::allocator<int>>,std::allocator<boost::multi::array<int,2,std::allocator<int>>>>(const std::vector<boost::multi::array<int,2,std::allocator<int>>,std::allocator<boost::multi::array<int,2,std::allocator<int>>>> &,const std::vector<boost::multi::array<int,2,std::allocator<int>>,std::allocator<boost::multi::array<int,2,std::allocator<int>>>> &)' being compiled
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\vector(2275): note: see reference to function template instantiation 'bool std::equal<const _Ty*,const _Ty*>(const _InIt1,const _InIt1,const _InIt2)' being compiled
			with
			[
				_Ty=boost::multi::array<int,2,std::allocator<int>>,
				_InIt1=const boost::multi::array<int,2,std::allocator<int>> *,
				_InIt2=const boost::multi::array<int,2,std::allocator<int>> *
			]
	C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.50.35717\include\xutility(5663): note: see reference to function template instantiation 'bool std::equal<_InIt1,_InIt2,std::equal_to<void>>(const _InIt1,const _InIt1,const _InIt2,_Pr)' being compiled
			with
			[
				_InIt1=const boost::multi::array<int,2,std::allocator<int>> *,
				_InIt2=const boost::multi::array<int,2,std::allocator<int>> *,
				_Pr=std::equal_to<void>
			]

		call "..\..\..\bin.v2\standalone\msvc\msvc-14.3\msvc-setup.bat"  >nul
	 cl /Zm800 -nologo "allocator.cpp" -c -Fo"..\..\..\bin.v2\libs\boost-multi\test\allocator.test\msvc-14.3\debug\cxxstd-20-iso\threading-multi\allocator.obj"     -TP /wd4675 /EHs /std:c++20 /GR /Zc:throwingNew /Z7 /Od /Ob0 /W4 /WX /MDd /Zc:forScope /Zc:wchar_t /Zc:inline -DBOOST_ALL_NO_LIB=1 -DBOOST_COBALT_USE_STD_PMR=1 "-I..\..\.."
	*/

	friend constexpr auto operator==(array_ref const& self, array_ref const& other) -> bool {
		if(self.extents() != other.extents()) {
			return false;
		}
		return adl_equal(
			other.data_elements(), other.data_elements() + other.num_elements(),
			self.data_elements()
		);
	}

#endif

	template<typename TT, class... As>
	friend constexpr auto operator!=(array_ref const& self, array_ref<TT, D, As...> const& other) -> bool {
		if(self.extents() != other.extents()) {
			return true;
		}
		return !adl_equal(
			other.data_elements(), other.data_elements() + other.num_elements(),  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) TODO(correaa) use span?
			self.data_elements()
		);
		// return ! operator==(self, other);  // commented due to bug in nvcc 22.11
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_HD constexpr auto data_elements() & -> typename array_ref::element_ptr { return array_ref::base_; }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0
	BOOST_MULTI_HD constexpr auto data_elements() && -> typename array_ref::element_ptr { return array_ref::base_; }

	friend constexpr auto data_elements(array_ref&& self) -> typename array_ref::element_ptr { return std::move(self).data_elements(); }

	/// for `D == 1` only, gives a pointer to a memory range of `.size()`, for compatibility with `std::vector` and `std::span`.
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data() const& { return data_elements(); }  // NOLINT(modernize-use-constraints) for C++20
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data() && { return data_elements(); }      // NOLINT(modernize-use-constraints) for C++20
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data() & { return data_elements(); }       // NOLINT(modernize-use-constraints) for C++20

	// // TODO(correaa) : find a way to use [[deprecated("use data_elements()")]] for friend functions
	// friend constexpr auto data(array_ref const& self) -> typename array_ref::element_ptr { return self.data_elements(); }
	// friend constexpr auto data(array_ref& self) -> typename array_ref::element_ptr { return self.data_elements(); }
	// friend constexpr auto data(array_ref&& self) -> typename array_ref::element_ptr { return std::move(self).data_elements(); }

	using decay_type = typename array_ref::decay_type;

	/// materializes an independent, owning `array` copy of this view with the associated array-value type (use unary prefix `+` as a shortcut)
	constexpr auto decay() const& -> decay_type const& { return static_cast<decay_type const&>(*this); }  // cppcheck-suppress duplInheritedMember ; to override

 private:
	template<class TTN, std::size_t DD = 0>
	void check_sizes_() const {
		using std::get;  // for C++17 compatibility
		if(size_type{get<DD>(this->sizes())} != size_type{std::extent_v<TTN, unsigned{DD}>}) {
			throw std::bad_cast{};
		}
		if constexpr(DD + 1 != D) {
			check_sizes_<TTN, DD + 1>();
		}
	}

	template<class TT> static auto launder_(TT* pointer) -> TT* {
#if defined(__cpp_lib_launder) && (__cpp_lib_launder >= 201606L)
		return std::launder(pointer);
#else
		return pointer;
#endif
	}

	template<typename, ::boost::multi::dimensionality_type, class> friend struct ::boost::multi::array;

	template<class TTN>
	constexpr auto to_carray_() & -> TTN& {
		check_sizes_<TTN>();
		return *launder_(reinterpret_cast<TTN*>(array_ref::base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	template<typename To>
	constexpr auto to_carray_() const& -> To const& {
		check_sizes_<To>();
		return *launder_(reinterpret_cast<To const*>(array_ref::base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

 public:
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	template<typename To, std::enable_if_t<std::is_array_v<To>, int> = 0>        // NOLINT(modernize-use-constraints) for C++20
	constexpr explicit operator To const&() const& { return to_carray_<To>(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

	template<typename To, std::enable_if_t<std::is_array_v<To>, int> = 0>  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-use-constraints) for C++20
	constexpr explicit operator To&() && { return to_carray_<To>(); }      // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)

	template<typename To, std::enable_if_t<std::is_array_v<To>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr explicit operator To&() & { return to_carray_<To>(); }       // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

 private:
	template<class Ar>
	auto serialize_structured_(Ar& arxiv, unsigned int const version) {
		subarray_base::serialize(arxiv, version);
	}

	template<class Archive>
	auto serialize_flat_(Archive& arxiv, unsigned int const /*version*/) {
		using AT = multi::archive_traits<Archive>;
		arxiv& AT::make_nvp("elements", AT::make_array(this->data_elements(), static_cast<std::size_t>(this->num_elements())));
	}

	//  template<class Ar, class AT = multi::archive_traits<Ar>>
	//  auto serialize_binary_if(std::true_type, Ar& ar) {
	//      ar & AT::make_nvp("binary_data", AT::make_binary_object(this->data_elements(), static_cast<std::size_t>(this->num_elements())*sizeof(typename array_ref::element)));
	//  }
	//  template<class Ar>
	//  auto serialize_binary_if(std::false_type, Ar& ar) {return serialize_flat(ar);}

 public:
	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int const version) {  // cppcheck-suppress duplInheritedMember ;
		serialize_flat_(arxiv, version);
		//      serialize_structured_(ar, version);
		//      switch(version) {
		//          case static_cast<unsigned int>( 0): return serialize_flat_(arxiv);
		//          case static_cast<unsigned int>(-1): return serialize_structured_(arxiv, version);
		//      //  case 2: return serialize_binary_if(std::is_trivially_copy_assignable<typename array_ref::element>{}, arxiv);
		//          default:
		//              if( this->num_elements() <= version ){serialize_structured_(arxiv, version);}
		//              else                                 {serialize_flat_       (arxiv         );}
		//      }
	}
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

/// Convenience alias for a D‐dimensional view of a contiguous, pre‐existing memory constant memory buffer
template<class T, dimensionality_type D, class Ptr = typename std::pointer_traits<T*>::template rebind<T const>>
using array_cref = array_ref<std::decay_t<T>, D, Ptr>;

template<class T, dimensionality_type D, class Ptr = T*>
using array_mref = array_ref<
	std::decay_t<T>, D,
	std::move_iterator<Ptr>>;

template<class TT, std::size_t N>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) interact with legacy  // NOSONAR
constexpr auto ref(TT (&arr)[N]) {
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) interact with legacy
	return array_ref<std::remove_all_extents_t<TT[N]>, std::rank_v<TT[N]>>(arr);
}

namespace detail {
template<class T, dimensionality_type D, typename Ptr /*= T* */>
struct array_ptr
: detail::subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_type, false> {
	using basic_ptr = detail::subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_type, false>;

	constexpr array_ptr(Ptr data, multi::extents_t<D> extensions)
	: basic_ptr{data, typename array_ref<T, D, Ptr>::layout_type(extensions)} {}

	constexpr explicit array_ptr(std::nullptr_t nil) : array_ptr{nil, multi::extents_t<D>{}} {}

	template<typename CArray>
	// cppcheck-suppress constParameterPointer ;  workaround cppcheck 2.11
	constexpr explicit array_ptr(CArray* data) : array_ptr{data_elements(*data), extensions(*data)} {}

	template<
		class TT, std::size_t N,
		std::enable_if_t<std::is_convertible_v<decltype(data_elements(std::declval<TT (&)[N]>())), Ptr>, int> = 0  // NOLINT(modernize-use-constraints,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) support legacy c-arrays TODO(correaa) for C++20
		>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr array_ptr(TT (*array)[N]) : array_ptr{data_elements(*array), extensions(*array)} {}  // NOLINT(modernize-use-constraints,google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) array_ptr is more general than pointer c-array support legacy c-arrays  TODO(correaa) for C++20  // NOSONAR

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto operator*() const -> array_ref<T, D, Ptr> {
		return array_ref<T, D, Ptr>(static_cast<detail::subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_type, false> const&>(*this)->extents(), this->base());
		// return array_ref<T, D, Ptr>((*static_cast<detail::subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_type, false> const&>(*this)).extents(), this->base());  // NOLINT(readability-redundant-parentheses) bug in clang-tidy
	}
};

template<class T, typename Ptr>
class /*[[deprecated("no good uses found")]]*/ array_ptr<T, 0, Ptr> {  // TODO(correaa) make it private mutable member
	mutable multi::array_ref<T, 0, Ptr> ref_;                          // TODO(correaa) implement array_ptr like other cases

 public:
	~array_ptr() = default;

	constexpr array_ptr(array_ptr const&)     = default;
	constexpr array_ptr(array_ptr&&) noexcept = default;  // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor) TODO(correaa) change the implementation like the other cases

	constexpr explicit array_ptr(Ptr dat, typename multi::array_ref<T, 0, Ptr>::extensions_type extensions) : ref_(dat, extensions) {}
	constexpr explicit array_ptr(Ptr dat) : array_ptr(dat, typename multi::array_ref<T, 0, Ptr>::extensions_type{}) {}

	constexpr explicit operator bool() const { return ref_.base(); }  // cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr explicit operator Ptr() const { return ref_.base(); }

	auto operator=(array_ptr const&) -> array_ptr&     = default;
	auto operator=(array_ptr&&) noexcept -> array_ptr& = default;

	friend constexpr auto operator==(array_ptr const& self, array_ptr const& other) -> bool { return self.ref_.base() == other.ref_.base(); }
	friend constexpr auto operator!=(array_ptr const& self, array_ptr const& other) -> bool { return self.ref_.base() != other.ref_.base(); }

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto operator*() const -> multi::array_ref<T, 0, Ptr>& { return ref_; }  // moLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) make ref base class a mutable member

	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr auto operator->() const -> multi::array_ref<T, 0, Ptr>* { return &ref_; }  // moLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) make ref base class a mutable member
};

template<class TT, std::size_t N>
constexpr auto addressof(TT (&array)[N]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return detail::array_ptr<
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		std::decay_t<std::remove_all_extents_t<TT[N]>>, static_cast<dimensionality_type>(std::rank<TT[N]>{}), std::remove_all_extents_t<TT[N]>*>{&array};
}
}  // namespace detail

template<class T, dimensionality_type D, typename Ptr = T*>
using array_ptr [[deprecated]] = detail::array_ptr<T, D, Ptr>;

// template<dimensionality_type D, class P>
// constexpr auto make_array_ref(P data, multi::extents_t<D> extensions) {
// 	return array_ref<typename std::iterator_traits<P>::value_type, D, P>(data, extensions);
// }

// template<class P> auto make_array_ref(P data, extents_t<0> exts) { return make_array_ref<0>(data, exts); }
// template<class P> auto make_array_ref(P data, extents_t<1> exts) { return make_array_ref<1>(data, exts); }
// template<class P> auto make_array_ref(P data, extents_t<2> exts) { return make_array_ref<2>(data, exts); }
// template<class P> auto make_array_ref(P data, extents_t<3> exts) { return make_array_ref<3>(data, exts); }
// template<class P> auto make_array_ref(P data, extents_t<4> exts) { return make_array_ref<4>(data, exts); }
// template<class P> auto make_array_ref(P data, extents_t<5> exts) { return make_array_ref<5>(data, exts); }

#ifdef __cpp_deduction_guides
namespace detail {
template<class It, typename V = typename std::iterator_traits<It>::value_type>  // pointer_traits doesn't have ::value_type
array_ptr(It) -> array_ptr<V, 0, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>  // pointer_traits doesn't have ::value_type
array_ptr(It, index_extensions<0>) -> array_ptr<V, 0, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<1>) -> array_ptr<V, 1, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<2>) -> array_ptr<V, 2, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<3>) -> array_ptr<V, 3, It>;
#endif

#ifdef __cpp_deduction_guides
template<
	class T,
	std::size_t N,
	typename V = std::remove_all_extents_t<T[N]>, std::size_t D = std::rank_v<T[N]>  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	>
array_ptr(T (*)[N]) -> array_ptr<V, static_cast<multi::dimensionality_type>(D)>;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
#endif
}  // end namespace detail

#ifdef __cpp_deduction_guides
// template<typename It, class Tuple> array_ref(It, Tuple) -> array_ref<typename std::iterator_traits<It>::value_type, std::tuple_size_v<Tuple>, It>;

template<typename Ptr> array_ref(Ptr, extents_t<0>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 0, Ptr>;
template<typename Ptr> array_ref(Ptr, extents_t<1>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 1, Ptr>;
template<typename Ptr> array_ref(Ptr, extents_t<2>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 2, Ptr>;
template<typename Ptr> array_ref(Ptr, extents_t<3>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 3, Ptr>;
template<typename Ptr> array_ref(Ptr, extents_t<4>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 4, Ptr>;
template<typename Ptr> array_ref(Ptr, extents_t<5>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 5, Ptr>;
#endif

#ifdef __cpp_deduction_guides
template<class It> const_subarray(It, It) -> const_subarray<typename It::element, It::dimensionality + 1, typename It::element_ptr, layout_t<It::dimensionality + 1>>;

template<class T> const_subarray(std::initializer_list<T>) -> const_subarray<T, 1>;
template<class T> const_subarray(std::initializer_list<std::initializer_list<T>>) -> const_subarray<T, 2>;
#endif

// TODO(correaa) move to utility
template<class T, std::size_t N>
constexpr auto rotated(T const (&array)[N]) noexcept {                                                 // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(array))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
			   base(array), extensions(array)
	)
		.rotated();
}

template<class T, std::size_t N>
constexpr auto rotated(T (&array)[N]) noexcept {                                                       // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(array))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
			   base(array), extensions(array)
	)
		.rotated();
}

template<class RandomAccessIterator, dimensionality_type D>
constexpr auto operator/(RandomAccessIterator data, multi::extents_t<D> extensions)
	-> detail::array_ptr<typename std::iterator_traits<RandomAccessIterator>::value_type, D, RandomAccessIterator> {
	return {data, extensions};
}

namespace detail {
#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

// adl_uninitialized_copy's priority<3> overload finds it via unqualified lookup/ADL, the same path a user-provided uninitialized_copy override would take;
// this is what recurses into the N-1 dimension instead of letting std::uninitialized_copy placement-new a subarray proxy
// TODO(correaa) put this in adl.hpp with its own priority
template<class In, class T, dimensionality_type N, class TP, class = std::enable_if_t<(N > 1)>, class = decltype((void)adl_begin(*std::declval<In&>()), adl_end(*std::declval<In&>()))>
constexpr auto uninitialized_copy
	// require N>1 (this is important because it forces calling placement new on the pointer
	(In first, In last, multi::detail::array_iterator<T, N, TP> dest) {  // NOLINT(performance-unnecessary-value-param) TODO(correaa) inverstigate why I can't make this In const& last
	while(first != last) {                                               // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		adl_uninitialized_copy(adl_begin(*first), adl_end(*first), adl_begin(*dest));
		++first;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		++dest;
	}
	return dest;
}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif
}  // end namespace detail

// this multi::size has lower priority than std::size because of the T&&
template<class T> constexpr auto size(T&& rng) -> decltype(std::forward<T>(rng).size()) { return std::forward<T>(rng).size(); }

// begin and end for forwarding reference are needed in this namespace
// to overwrite the behavior of std::begin and std::end
// which take rvalue-references as const-references.

template<class T> constexpr auto begin(T&& rng) -> decltype(std::forward<T>(rng).begin()) { return std::forward<T>(rng).begin(); }  ///< Returns the beginning of the range (generic free function, usually return .begin())
template<class T> constexpr auto end(T&& rng) -> decltype(std::forward<T>(rng).end()) { return std::forward<T>(rng).end(); }        ///< Returns the end of the range (generic free function, usually return .begin())

// this has to take argument by forward reference to avoid collison with std::cbegin/std::cend
template<class T> constexpr auto cbegin(T&& rng) -> decltype(std::forward<T>(rng).begin()) { return std::forward<T>(rng).begin(); }  ///< Returns the beginning of the range (for constant access, generic free function, usually return .begin())
template<class T> constexpr auto cend(T&& rng) -> decltype(std::forward<T>(rng).end()) { return std::forward<T>(rng).end(); }        ///< Returns the end of the range (for constant access, generic free function, usually return .begin())

template<class T> constexpr auto                    extent(T const& rng) -> decltype(rng.extent()) { return rng.extent(); }
template<class T> [[deprecated("use extent")]] auto extension(T const& rng) -> decltype(rng.extent()) { return rng.extent(); }

template<class T> constexpr auto stride(T const& rng) -> decltype(rng.stride()) { return rng.stride(); }

template<class T, std::size_t N, std::size_t M>
auto transposed(T (&array)[N][M]) -> decltype(auto) { return ~multi::array_ref<T, 2>(array); }  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

}  // end namespace boost::multi

#ifndef BOOST_MULTI_SERIALIZATION_ARRAY_DEFAULT_VERSION
#define BOOST_MULTI_SERIALIZATION_ARRAY_DEFAULT_VERSION 0  // NOLINT(cppcoreguidelines-macro-usage) reserved for future use to select default format for serialization version //NOSONAR
#endif

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
namespace std::ranges {  // NOLINT(cert-dcl58-cpp) to enable borrowed, nvcc needs namespace
template<typename Element, ::boost::multi::dimensionality_type D, class... Rest>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::subarray<Element, D, Rest...>> = true;  // NOLINT(misc-definitions-in-headers)

template<typename Element, ::boost::multi::dimensionality_type D, class... Rest>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::const_subarray<Element, D, Rest...>> = true;  // NOLINT(misc-definitions-in-headers)

template<typename Element, ::boost::multi::dimensionality_type D, class... Rest>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::array_ref<Element, D, Rest...>> = true;  // NOLINT(misc-definitions-in-headers)

template<typename Ptr, class... Rest>
[[maybe_unused]] inline constexpr bool enable_borrowed_range<::boost::multi::detail::elements_range_t<Ptr, Rest...>> = true;  // NOLINT(misc-definitions-in-headers)
}  // end namespace std::ranges
#endif

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#undef BOOST_MULTI_HD
#undef BOOST_MULTI_NO_DANGLING

#endif  // BOOST_MULTI_ARRAY_REF_HPP
