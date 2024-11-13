// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ARRAY_REF_HPP_
#define BOOST_MULTI_ARRAY_REF_HPP_
#pragma once

#include <boost/multi/detail/tuple_zip.hpp>
#include <boost/multi/utility.hpp>  // IWYU pragma: export

namespace boost::multi {

template<class Element>
inline constexpr bool force_element_trivial = false;

template <class Element>
inline constexpr bool force_element_trivial_destruction = force_element_trivial<Element>;

template <class Element>
inline constexpr bool force_element_trivial_default_construction = force_element_trivial<Element>;

#ifdef _MULTI_FORCE_TRIVIAL_STD_COMPLEX
template<class T>
inline constexpr bool force_element_trivial<std::complex<T>> = std::is_trivial_v<T>;

template<class T>
inline constexpr bool force_element_trivial_destruction<std::complex<T>> = std::is_trivially_default_constructible_v<T>;

template<class T>
inline constexpr bool force_element_trivial_default_construction<std::complex<T>> = std::is_trivially_destructible_v<T>;

template<> inline constexpr bool force_element_trivial                     <std::complex<double>> = true;
template<> inline constexpr bool force_element_trivial_default_construction<std::complex<double>> = true;
template<> inline constexpr bool force_element_trivial_destruction         <std::complex<double>> = true;

template<> inline constexpr bool force_element_trivial                     <std::complex<float >> = true;
template<> inline constexpr bool force_element_trivial_default_construction<std::complex<float >> = true;
template<> inline constexpr bool force_element_trivial_destruction         <std::complex<float >> = true;
#endif

}  // end namespace boost::multi

#include <boost/multi/detail/adl.hpp>  // TODO(correaa) remove instantiation of force_element_trivial in this header
#include <boost/multi/detail/layout.hpp>          // IWYU pragma: export
#include <boost/multi/detail/memory.hpp>          // for pointer_traits
#include <boost/multi/detail/operators.hpp>       // for random_iterable
#include <boost/multi/detail/pointer_traits.hpp>  // IWYU pragma: export
#include <boost/multi/detail/serialization.hpp>
#include <boost/multi/detail/types.hpp>           // for dimensionality_type  // IWYU pragma: export

#include <boost/multi/detail/config/ASSERT.hpp>

#include <algorithm>   // fpr copy_n
#include <array>
#include <cstring>     // for std::memset in reinterpret_cast
#include <functional>  // for std::invoke
#include <iterator>    // for std::next
#include <memory>      // for std::pointer_traits
#include <new>         // for std::launder

#if __has_include(<span>)
#if !defined(_MSVC_LANG) || (_MSVC_LANG > 202002L)
#include <span>
#endif
#if defined(__cpp_lib_span) && __cpp_lib_span >= 202002L && !defined(_MSVC_LANG)
#define BOOST_MULTI_HAS_SPAN
#endif
#endif

#include <utility>     // for forward

#if !defined(__NVCC__)
	#define BOOST_MULTI_FRIEND_CONSTEXPR friend   constexpr  // this generates a problem with intel compiler 19 and v2021 "a constexpr function cannot have a nonliteral return type"
#else
	#define BOOST_MULTI_FRIEND_CONSTEXPR friend /*constexpr*/
#endif

#if defined(__NVCC__)
	#define BOOST_MULTI_HD __host__ __device__
#else
	#define BOOST_MULTI_HD
#endif

// #if defined(__APPLE__) && !defined(__clang__) /*&& defined(__aarch64__)*/ && defined(BOOST_UNIT_TEST_FRAMEWORK_DYN_LINK)
//  #include <boost/test/tools/detail/print_helper.hpp>
//  #include <boost/test/tree/test_unit.hpp>
//  #include <boost/test/unit_test_suite.hpp>
// auto boost::test_tools::tt_detail::print_log_value<char>::operator()(std::ostream&, char) -> void {}
// auto boost::unit_test::test_unit::full_name() const -> std::string { return std::string{}; }
// auto boost::unit_test::ut_detail::normalize_test_case_name(boost::unit_test::basic_cstring<char const> name) -> std::string {
//  return std::string(name.begin(), name.end());
// }
// namespace boost::unit_test {
// // void make_test_case(boost::function<void ()> const& f, boost::unit_test::basic_cstring<char const>, boost::unit_test::basic_cstring<char const>, unsigned long) {
// //  f();
// // }
// }
// #endif

namespace boost::multi {

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>>
struct const_subarray;

// template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>>
// using subarray = const_subarray<T, D, ElementPtr, Layout>;

template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>>
class subarray;

template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D, typename std::pointer_traits<ElementPtr>::difference_type>>
class move_subarray;

template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
constexpr auto is_subarray_aux(const_subarray<T, D, ElementPtr, Layout> const&) -> std::true_type;
constexpr auto is_subarray_aux(...                                            ) -> std::false_type;

template<class A> struct is_subarray: decltype(is_subarray_aux(std::declval<A>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<dimensionality_type D>
struct of_dim {
template<typename T, class ElementPtr, class Layout>
static constexpr auto is_subarray_of_dim_aux(subarray<T, D, ElementPtr, Layout> const&) -> std::true_type;
static constexpr auto is_subarray_of_dim_aux(...                                      ) -> std::false_type;

template<class A> struct is_subarray_of_dim: decltype(is_subarray_of_dim_aux(std::declval<A>())) {};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
};

template<typename T, dimensionality_type D, class A = std::allocator<T>> struct array;

template<
	typename T,
	dimensionality_type D,
	typename ElementPtr = T*,
	class Layout = layout_t<D, std::make_signed_t<typename std::pointer_traits<ElementPtr>::size_type> >
>
struct array_types : private Layout {  // cppcheck-suppress syntaxError ; false positive in cppcheck
	using element = T;
	using element_type = element;  // this follows more closely https://en.cppreference.com/w/cpp/memory/pointer_traits

	using element_ptr       = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element const>;
	using element_move_ptr  = multi::move_ptr<element, element_ptr>;

	using element_ref = typename std::iterator_traits<element_ptr>::reference;

	using layout_t = Layout;

	using rank = typename layout_t::rank;

	using layout_t::rank_v;

	using dimensionality_type = typename layout_t::dimensionality_type;
	using                                layout_t::dimensionality;

	using typename layout_t::stride_type;
	using          layout_t::stride     ;

	using layout_t::num_elements;
	using layout_t::offset;

	using layout_t::offsets;

	using typename layout_t::index;
	using typename layout_t::index_range;
	using typename layout_t::index_extension;

	using typename layout_t::strides_type;

	auto strides() const { return detail::convertible_tuple<strides_type>(layout_t::strides()); }

	using typename layout_t::difference_type;

	using typename layout_t::size_type;
	using          layout_t::size     ;

	using layout_t::nelems;

	using typename layout_t::extension_type;
	using          layout_t::extension;

	using typename layout_t::extensions_type;
	using          layout_t::extensions;

	constexpr auto extensions() const -> extensions_type {return static_cast<layout_t const&>(*this).extensions();}

	using layout_t::is_empty;
	using layout_t::   empty;

	using layout_t::sub;

	using typename layout_t::sizes_type;
	using          layout_t::sizes;

	[[deprecated("This is for compatiblity with Boost.MultiArray, you can use `rank` member type or `dimensionality` static member variable")]]
	static constexpr auto num_dimensions() { return dimensionality; }

	[[deprecated("This is for compatiblity with Boost.MultiArray, you can use `offsets` member function")]]
	auto index_bases() const -> std::ptrdiff_t const*;  // = delete;  this function is not implemented, it can give a linker error

	[[deprecated("This is for compatiblity with Boost.MultiArray, you can use `offsets` member function")]]
	constexpr auto shape() const { return detail::convertible_tuple(this->sizes()); }

	using layout_t::is_compact;

	friend constexpr auto size        (array_types const& self) noexcept -> size_type       {return self.size        ();}
	friend constexpr auto extension   (array_types const& self) noexcept -> extension_type  {return self.extension   ();}
	friend constexpr auto is_empty    (array_types const& self) noexcept -> bool            {return self.is_empty    ();}
	friend constexpr auto num_elements(array_types const& self) noexcept -> size_type       {return self.num_elements();}

	friend constexpr auto extensions  (array_types const& self) noexcept -> extensions_type {return self.extensions  ();}
	friend constexpr auto sizes       (array_types const& self) noexcept -> sizes_type      {return self.sizes       ();}

	// TODO(correaa) [[deprecated("use member syntax for non-salient properties")]]
	friend
	constexpr auto stride     (array_types const& self) noexcept -> stride_type         {return self.stride     ();}

	// TODO(correaa) [[deprecated("use member syntax for non-salient properties")]]
	friend
	constexpr auto strides    (array_types const& self) noexcept -> strides_type         {return self.strides     ();}

 protected:
	constexpr auto layout_mutable() -> layout_t& {return static_cast<layout_t&>(*this);}

 public:
	using value_type = typename std::conditional_t<
		(D > 1),
		array<element, D-1, typename multi::pointer_traits<element_ptr>::default_allocator_type>,
		element
	>;

	using reference = typename std::conditional_t<
		(D > 1),
		subarray<element, D-1, element_ptr>,
		typename std::iterator_traits<element_ptr>::reference
	>;

	using const_reference = typename std::conditional_t<
		(D > 1),
		const_subarray<element, D-1, element_ptr>,
		typename std::iterator_traits<element_const_ptr>::reference
	>;

	// BOOST_MULTI_HD constexpr auto base() &      -> element_ptr       {return base_;}
	// BOOST_MULTI_HD constexpr auto base() &&     -> element_ptr       {return base_;}
	// BOOST_MULTI_HD constexpr auto base() const& -> element_const_ptr {return base_;}
	BOOST_MULTI_HD constexpr auto base() const -> element_const_ptr { return base_; }

	BOOST_MULTI_HD constexpr auto  mutable_base() const -> element_ptr {return base_;}

	BOOST_MULTI_HD constexpr auto cbase() const  -> element_const_ptr {return base_;}
	BOOST_MULTI_HD constexpr auto mbase() const& -> element_ptr&      {return base_;}

	// friend /*constexpr*/ auto  base(array_types & self) -> element_ptr  {return self.base_;}
	// friend /*constexpr*/ auto  base(array_types && self) -> element_ptr  {return std::move(self).base_;}
	// friend /*constexpr*/ auto  base(array_types const& self) -> element_const_ptr  {return self.base_;}

	    BOOST_MULTI_HD constexpr auto layout()           const        -> layout_t const& {return *this;}
	friend constexpr auto layout(array_types const& self) -> layout_t const& {return self.layout();}

	       constexpr auto origin()           const&       -> decltype(auto) {return base_ + Layout::origin();}
	friend constexpr auto origin(array_types const& self) -> decltype(auto) {return self.origin();}

 protected:
	BOOST_MULTI_NO_UNIQUE_ADDRESS
	element_ptr base_;  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes,misc-non-private-member-variables-in-classes) : TODO(correaa) try to make it private, [static_]array needs mutation
	
	template<class, ::boost::multi::dimensionality_type, typename, bool, bool> friend struct array_iterator;

	using derived = subarray<T, D, ElementPtr, Layout>;
	BOOST_MULTI_HD constexpr explicit array_types(std::nullptr_t) : Layout{}, base_(nullptr) {}

 public:
	array_types() = default;

	BOOST_MULTI_HD constexpr array_types(layout_t const& lyt, element_ptr const& data)
	: Layout{lyt}, base_{data} {}

 protected:
	template<
		class ArrayTypes,
		typename = std::enable_if_t<! std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>
		, decltype(multi::detail::explicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	>
	// underlying pointers are explicitly convertible
	BOOST_MULTI_HD constexpr explicit array_types(ArrayTypes const& other)
	: Layout{other.layout()}, base_{other.base_} {}

	template<
		class ArrayTypes,
		typename = std::enable_if_t<!std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>,
		decltype(multi::detail::implicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointers are implicitly convertible
	BOOST_MULTI_HD constexpr /*implt*/ array_types(ArrayTypes const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : inherit behavior of underlying pointer
	: Layout{other.layout()}, base_{other.base_} {}

	template<
		typename ElementPtr2,
		typename = decltype(Layout{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().layout()}),
		typename = decltype(element_ptr{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().base_})
	>
	BOOST_MULTI_HD constexpr explicit array_types(array_types<T, D, ElementPtr2, Layout> const& other)
	: Layout{other.layout()}, base_{other.base_} {}

	template<class T2, ::boost::multi::dimensionality_type D2, class E2, class L2> friend struct array_types;
};

template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout, bool IsConst = false>
struct subarray_ptr;

template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = multi::layout_t<D> >
using const_subarray_ptr = subarray_ptr<T, D, ElementPtr, Layout, true>;

template<typename T, multi::dimensionality_type D, typename ElementPtr = T*, class Layout = multi::layout_t<D>, bool IsConst>
struct subarray_ptr  // NOLINT(fuchsia-multiple-inheritance) : to allow mixin CRTP
: boost::multi::iterator_facade<
	subarray_ptr<T, D, ElementPtr, Layout>, void, std::random_access_iterator_tag,
	subarray<T, D, ElementPtr, Layout> const&, typename Layout::difference_type
> {
 private:
	Layout layout_;
	ElementPtr base_;
	typename std::iterator_traits<ElementPtr>::difference_type offset_;

 public:
	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;
	template<typename, multi::dimensionality_type, typename, bool, bool> friend struct array_iterator;

	~subarray_ptr() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	using pointer         = subarray<T, D, ElementPtr, Layout>*;
	using element_type    = typename subarray<T, D, ElementPtr, Layout>::decay_type;
	using difference_type = typename Layout::difference_type;

	using value_type = element_type;

	using reference = std::conditional_t<
		IsConst,
		const_subarray<T, D, ElementPtr, Layout>,
		subarray<T, D, ElementPtr, Layout>>;

	using iterator_category = std::random_access_iterator_tag;

	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr subarray_ptr(std::nullptr_t nil) : layout_{}, base_{nil}, offset_{0} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) terse syntax and functionality by default
	
	subarray_ptr() = default;

	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;

	BOOST_MULTI_HD constexpr subarray_ptr(typename reference::element_ptr base, layout_t<typename reference::rank{} - 1> lyt) : layout_{lyt}, base_{base}, offset_{0} {}

	template<bool OtherIsConst, std::enable_if_t< ! OtherIsConst, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr/*mplct*/ subarray_ptr(subarray_ptr<T, D, ElementPtr, Layout, OtherIsConst> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : propagate implicitness of pointer
	: layout_{other.layout_}, base_{other.base_}, offset_{other.offset_} {}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherLayout, bool OtherIsConst,
		decltype(multi::detail::implicit_cast<typename reference::element_ptr>(std::declval<OtherEPtr>()))* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr/*mplct*/ subarray_ptr(subarray_ptr<OtherT, OtherD, OtherEPtr, OtherLayout, OtherIsConst> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : propagate implicitness of pointer
	: layout_{other.layout_}, base_{other.base_} {}

	template<
		class ElementPtr2, 
		std::enable_if_t<std::is_same_v<ElementPtr2, ElementPtr> && (D==0), int> =0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	BOOST_MULTI_HD constexpr explicit subarray_ptr(ElementPtr2 const& other) : layout_{}, base_{other} {}

	// template<class Array>
	// // cppcheck-suppress noExplicitConstructor ; no information loss, allows comparisons
	// BOOST_MULTI_HD constexpr subarray_ptr(Array* other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	// : subarray_ptr(other->data_elements(), other->layout()) {}
	// BOOST_MULTI_HD constexpr explicit subarray_ptr(subarray_ptr<T, D, ElementPtr, Layout>* other)
	// : subarray_ptr(other->base(), other->layout()) {}

	subarray_ptr(subarray_ptr const&) = default;
	subarray_ptr(subarray_ptr     &&) noexcept = default;  // TODO(correaa) remove inheritnace from reference to remove this move ctor

	auto operator=(subarray_ptr const&) -> subarray_ptr& = default;
	auto operator=(subarray_ptr     &&) noexcept -> subarray_ptr& = default;

	BOOST_MULTI_HD constexpr explicit operator bool() const { return static_cast<bool>(base()); }

	BOOST_MULTI_HD constexpr auto operator*() const -> reference { return reference(layout_, base_); }

	// [[deprecated("->() is experimental on msvc")]]
	BOOST_MULTI_HD constexpr auto operator->() const {
		// static_assert( sizeof(*this) == sizeof(reference) );
		class proxy {
			reference ref_;

		 public:
			BOOST_MULTI_HD constexpr explicit proxy(reference&& ref) : ref_{std::move(ref)} {}
			BOOST_MULTI_HD constexpr auto operator->() && -> reference* {return std::addressof(this->ref_);}
		};
		return proxy{operator*()};
	}

	// BOOST_MULTI_HD constexpr auto operator->() -> reference* { 
	//  return std::addressof(reinterpret_cast<reference&>(*this));
	// }

	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> reference { return *(*this + n); }

	BOOST_MULTI_HD constexpr auto operator<(subarray_ptr const& other) const -> bool { return distance_to(other) > 0; }

	BOOST_MULTI_HD constexpr subarray_ptr(typename reference::element_ptr base, Layout const& lyt) : layout_{lyt}, base_{base} {}

	template<typename, multi::dimensionality_type, typename, class> friend struct const_subarray;

	BOOST_MULTI_HD constexpr auto base() const -> typename reference::element_ptr {return base_;}

	friend BOOST_MULTI_HD constexpr auto base(subarray_ptr const& self) {return self.base();}

	template<class OtherSubarrayPtr, std::enable_if_t<!std::is_base_of_v<subarray_ptr, OtherSubarrayPtr>, int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator==(OtherSubarrayPtr const& other) const
	->decltype((base_ == other.base_) && (layout_ == other.layout_)) {
		return (base_ == other.base_) && (layout_ == other.layout_); }

	constexpr auto operator==(subarray_ptr const& other) const -> bool {
		return (base_ == other.base_) && (layout_ == other.layout_);
	}

	constexpr auto operator!=(subarray_ptr const& other) const -> bool {
		return (base_ != other.base_) || (layout_ != other.layout_);
	}

	template<
		typename OtherT, multi::dimensionality_type OtherD, typename OtherEPtr, class OtherL, bool OtherIsConst,
		std::enable_if_t<!std::is_base_of_v<subarray_ptr, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst>>, int> =0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	friend BOOST_MULTI_HD constexpr auto operator!=(subarray_ptr const& self, subarray_ptr<OtherT, OtherD, OtherEPtr, OtherL, OtherIsConst> const& other) -> bool {
		return self.base_ != other.base_ || self.layout_ != other.layout_;
	}

 protected:
	BOOST_MULTI_HD constexpr void increment() { base_ += layout_.nelems(); }
	BOOST_MULTI_HD constexpr void decrement() { base_ -= layout_.nelems(); }

	BOOST_MULTI_HD constexpr void advance(difference_type n) { base_ += layout_.nelems()*n; }
	BOOST_MULTI_HD constexpr auto distance_to(subarray_ptr const& other) const -> difference_type {
		assert( layout_.nelems() == other.layout_.nelems() );
		// assert( Ref::nelems() == other.Ref::nelems() && Ref::nelems() != 0 );
		// assert( (other.base() - base())%Ref::nelems() == 0);
		assert( layout_ == other.layout_ );
		return (other.base_ - base_)/layout_.nelems();
	}

 public:
	BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> subarray_ptr& { advance(n); return *this; }
};

template<class Element, dimensionality_type D, typename ElementPtr, bool IsConst = false, bool IsMove = false>
struct array_iterator;

template<class Element, ::boost::multi::dimensionality_type D, typename ElementPtr, bool IsConst, bool IsMove>
struct array_iterator  // NOLINT(fuchsia-multiple-inheritance)
: boost::multi::iterator_facade<
	array_iterator<Element, D, ElementPtr, IsConst, IsMove>, void, std::random_access_iterator_tag,
	subarray<Element, D-1, ElementPtr> const&, typename layout_t<D-1>::difference_type
>
, multi::decrementable<array_iterator<Element, D, ElementPtr, IsConst>>
, multi::incrementable<array_iterator<Element, D, ElementPtr, IsConst>>
, multi::affine<array_iterator<Element, D, ElementPtr, IsConst>, multi::difference_type>
, multi::totally_ordered2<array_iterator<Element, D, ElementPtr, IsConst>, void> {
	~array_iterator() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	constexpr auto operator=(array_iterator&&)  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	-> array_iterator& = default;

	array_iterator(array_iterator&&) noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	= default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	using difference_type = typename layout_t<D>::difference_type;
	using element = Element;
	using element_ptr = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element const>;
	using value_type = typename subarray<element, D-1, element_ptr>::decay_type;

	using pointer         = subarray<element, D - 1, element_ptr>*;
	using reference       = std::conditional_t<
		IsConst,
		const_subarray<element, D - 1, element_ptr>,
		subarray<element, D - 1, element_ptr>
	>;
	using const_reference = const_subarray<element, D - 1, element_ptr>;  // TODO(correaa) should be const_subarray (base of subarray)

	using iterator_category = std::random_access_iterator_tag;

	constexpr static dimensionality_type rank_v = D;
	using rank = std::integral_constant<dimensionality_type, D>;  // TODO(correaa) make rank a function for compat with mdspan?

	using ptr_type = subarray_ptr<element, D-1, element_ptr, layout_t<D-1>>;

	using stride_type = index;
	using layout_type = typename reference::layout_type;  // layout_t<D - 1>

	// BOOST_MULTI_HD constexpr explicit array_iterator(std::nullptr_t nil) : ptr_{nil} {}
	// BOOST_MULTI_HD constexpr array_iterator() : array_iterator{nullptr} {}
	BOOST_MULTI_HD constexpr array_iterator() : ptr_{}, stride_{} {}  // = default;  // TODO(correaa) make = default, now it is not compiling

	template<class, dimensionality_type, class, bool, bool> friend struct array_iterator;

	template<
		class EElement, typename PPtr,
		decltype(multi::detail::explicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().base()))* = nullptr
	>
	BOOST_MULTI_HD constexpr explicit array_iterator(array_iterator<EElement, D, PPtr> const& other)
	: ptr_{element_ptr{other.base()}, other.ptr_->layout()}, stride_{other.stride_} {}

	template<class EElement, typename PPtr,
		decltype(multi::detail::implicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().base()))* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr/*mplct*/ array_iterator(array_iterator<EElement, D, PPtr> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : propagate implicitness of pointer
	: ptr_(other.ptr_), stride_{other.stride_} {}

	array_iterator(array_iterator const&) = default;
	auto operator=(array_iterator const&) -> array_iterator& = default;

	BOOST_MULTI_HD constexpr explicit operator bool() const {return ptr_->base();}  // TODO(correaa) implement bool conversion for subarray_ptr
	BOOST_MULTI_HD constexpr auto operator*() const -> reference {return *ptr_;}

	BOOST_MULTI_HD constexpr auto operator->() const -> decltype(auto) {return ptr_;}

	BOOST_MULTI_HD constexpr auto operator+ (difference_type n) const -> array_iterator {array_iterator ret{*this}; ret += n; return ret;}
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> subarray<element, D-1, element_ptr> {return *((*this) + n);}

	template<bool OtherIsConst, 
		std::enable_if_t<(IsConst != OtherIsConst), int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	BOOST_MULTI_HD constexpr auto operator==(array_iterator<Element, D, ElementPtr, OtherIsConst> const& other) const -> bool {
		// assert( this->stride_ == other.stride_ );
		// assert( this->ptr_->layout() == other.ptr_->layout() );
		return (this->ptr_ == other.ptr_) && (this->stride_ == other.stride_) && ( (*(this->ptr_)).layout() == (*(other.ptr_)).layout() );
	}

	BOOST_MULTI_HD constexpr auto operator==(array_iterator const& other) const -> bool {
		// assert( this->stride_ == other.stride_ );
		// assert( this->ptr_->layout() == other.ptr_->layout() );
		return (this->ptr_ == other.ptr_) && (this->stride_ == other.stride_) && ( (*(this->ptr_)).layout() == (*(other.ptr_)).layout() );
	}

	BOOST_MULTI_HD constexpr auto operator!=(array_iterator const& other) const -> bool {
		return !operator==(other);
	}

	// BOOST_MULTI_HD constexpr auto operator==(array_iterator<Element, D, ElementPtr, true> const& other) const -> bool {
	//  return this->ptr_ == other.ptr_ && this->stride_== other.stride_ && this->ptr_->layout() == other.ptr_->layout();
	// }

	BOOST_MULTI_HD constexpr auto operator< (array_iterator const& other) const -> bool {
		assert((*ptr_).layout() == (*(other.ptr_)).layout());
		assert(stride_ != 0);
		return
			   ((0 < stride_) && (ptr_.base() - other.ptr_.base() < 0))
			|| ((stride_ < 0) && (0 < ptr_.base() - other.ptr_.base()));  // TODO(correaa) consider the case where stride_ is negative
	}

	BOOST_MULTI_HD constexpr explicit array_iterator(typename subarray<element, D-1, element_ptr>::element_ptr base, layout_t<D-1> const& lyt, index stride)
	: ptr_{base, lyt}, stride_{stride} {}

	template<class, dimensionality_type, class, class> friend struct const_subarray;

	template<class... As>
	BOOST_MULTI_HD constexpr auto operator()(index idx, As... args) const -> decltype(auto) {return this->operator[](idx)(args...); }
	BOOST_MULTI_HD constexpr auto operator()(index idx)             const -> decltype(auto) {return this->operator[](idx)         ; }

 private:
	template<class Self, typename Tuple, std::size_t ... I>
	static BOOST_MULTI_HD constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...>/*012*/) -> decltype(auto) {
		using std::get;
		return std::forward<Self>(self)(get<I>(tuple)...);
	}

 public:
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl) const -> decltype(auto) { return apply_impl_(          *this , tpl, std::make_index_sequence<std::tuple_size<Tuple>::value>()); }
	// template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl)     && -> decltype(auto) { return apply_impl_(std::move(*this), tpl, std::make_index_sequence<std::tuple_size<Tuple>::value>()); }
	// template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl)      & -> decltype(auto) { return apply_impl_(          *this , tpl, std::make_index_sequence<std::tuple_size<Tuple>::value>()); }

 private:
	ptr_type ptr_;
	stride_type stride_;  // = {1};  // nice non-zero default  // TODO(correaa) use INT_MAX?  // TODO(correaa) remove to make type trivial

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD constexpr void decrement_() { ptr_.base_ -= stride_; }
	BOOST_MULTI_HD constexpr void advance_(difference_type n) { ptr_.base_ += stride_*n; }

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	BOOST_MULTI_HD constexpr auto base() const -> element_ptr {return ptr_.base_;}
	BOOST_MULTI_HD constexpr auto stride() const -> stride_type {return stride_;}

	friend /*constexpr*/ auto base(array_iterator const& self) -> element_ptr {return self.base();}  // TODO(correaa) remove
	friend constexpr auto stride(array_iterator const& self) -> stride_type {return self.stride_;}  // TODO(correaa) remove

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	constexpr auto operator++() -> array_iterator& { ptr_.base_ += stride_; return *this; }
	constexpr auto operator--() -> array_iterator& { ptr_.base_ -= stride_; return *this; }

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	friend constexpr auto operator-(array_iterator const& self, array_iterator const& other) -> difference_type {
		assert(self.stride_ == other.stride_);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) normal in a constexpr function
		assert(self.stride_ != 0);
		return (self.ptr_.base() - other.ptr_.base())/self.stride_;
	}

	constexpr auto operator+=(difference_type n) -> array_iterator& { advance_(+n); return *this; }
	constexpr auto operator-=(difference_type n) -> array_iterator& { advance_(-n); return *this; }
};

template<typename ElementPtr, dimensionality_type D, class StridesType>
struct cursor_t {
	using difference_type = typename std::iterator_traits<ElementPtr>::difference_type;
	using strides_type    = StridesType;

	using element_ptr  = ElementPtr;
	using element_ref  = typename std::iterator_traits<element_ptr>::reference;
	using element_type = typename std::iterator_traits<element_ptr>::value_type;

	using pointer   = element_ptr;
	using reference = element_ref;

	using indices_type = typename extensions_t<D>::indices_type;

	cursor_t() = default;

 private:
	strides_type strides_;
	element_ptr  base_;

	template<class, dimensionality_type, class, class> friend struct const_subarray;
	template<class, dimensionality_type, class> friend struct cursor_t;

	BOOST_MULTI_HD constexpr cursor_t(element_ptr base, strides_type const& strides) : strides_{strides}, base_{base} {}

	template<class OtherCursor, class = decltype(multi::detail::implicit_cast<element_ptr>(std::declval<OtherCursor>().base()))>
	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr          cursor_t(OtherCursor const& other) : strides_{other.strides()}, base_{other.base()} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class OtherCursor>
	BOOST_MULTI_HD constexpr explicit cursor_t(OtherCursor const& other) : strides_{other.strides()}, base_{other.base()} {}

 public:
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> decltype(auto) {
		using std::get;
		if constexpr(D != 1) {
			return cursor_t<
				ElementPtr,
				D-1,
				std::decay_t<decltype(strides_.tail())>
			>{
				base_
				+ get<0>(strides_)*n,
				strides_.tail()
			};
		} else {
			return base_[get<0>(strides_)*n];
		}
	}

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
		using std::get;
		return ((get<I>(tup)*get<I>(strides_)) + ...);
	}

 public:
	template<class Tuple = indices_type>
	BOOST_MULTI_HD constexpr auto operator+=(Tuple const& tup) -> cursor_t& {
		base_ += apply_impl_(tup, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator* () const -> reference {return *base_;}
	BOOST_MULTI_HD constexpr auto operator->() const -> pointer   {return  base_;}

	BOOST_MULTI_HD constexpr auto base() const -> pointer { return base_; }
	BOOST_MULTI_HD constexpr auto strides() const -> strides_type { return strides_; }
	template<multi::dimensionality_type DD = 0>
	BOOST_MULTI_HD constexpr auto stride() const { using std::get; return get<DD>(strides_); }
};

template<typename Pointer, class LayoutType>
struct elements_iterator_t  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
: boost::multi::random_accessable<elements_iterator_t<Pointer, LayoutType>, typename std::iterator_traits<Pointer>::difference_type, typename std::iterator_traits<Pointer>::reference>
{
	using difference_type = typename std::iterator_traits<Pointer>::difference_type;
	using value_type = typename std::iterator_traits<Pointer>::value_type;
	using pointer = Pointer;
	using reference =  typename std::iterator_traits<Pointer>::reference;
	using iterator_category = std::random_access_iterator_tag;

	using const_pointer = typename std::pointer_traits<pointer>::template rebind<value_type const>;

	using layout_type = LayoutType;

 private:
	pointer base_;
	layout_type l_;
	difference_type n_ = 0;
	extensions_t<layout_type::dimensionality> xs_;

	using indices_type = typename extensions_t<layout_type::dimensionality>::indices_type;
	indices_type ns_ = {};

	template<typename, class> friend struct elements_iterator_t;
	template<typename, class> friend struct elements_range_t;

	constexpr elements_iterator_t(pointer base, layout_type const& lyt, difference_type n)
	: base_{base}, l_{lyt}, n_{n}, xs_{l_.extensions()}, ns_{lyt.is_empty()?indices_type{}:xs_.from_linear(n)} {}

 public:
	elements_iterator_t() = default;

	BOOST_MULTI_HD constexpr auto base()       ->       pointer {return base_;}
	BOOST_MULTI_HD constexpr auto base() const -> const_pointer {return base_;}

	BOOST_MULTI_HD constexpr auto layout() const -> layout_type {return l_;}

	template<class Other, decltype(multi::detail::implicit_cast<pointer>(std::declval<Other>().base_))* = nullptr>
	// cppcheck-suppress noExplicitConstructor
	BOOST_MULTI_HD constexpr /*impl*/ elements_iterator_t(Other const& other) : elements_iterator_t{other.base_, other.l_, other.n_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class Other>
	BOOST_MULTI_HD constexpr explicit elements_iterator_t(Other const& other) : elements_iterator_t{other.base_, other.l_, other.n_} {}

	elements_iterator_t(elements_iterator_t const&) = default;

	BOOST_MULTI_HD constexpr auto operator=(elements_iterator_t const& other) -> elements_iterator_t& {  // fixes (?) warning: definition of implicit copy assignment operator for 'elements_iterator_t<boost::multi::array<double, 3> *, boost::multi::layout_t<1>>' is deprecated because it has a user-declared copy constructor [-Wdeprecated-copy]
		if(&other == this) {return *this;}  // for cert-oop54-cpp
		base_ = other.base_;
		xs_ = other.xs_;
		n_ = other.n_;
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator++() -> elements_iterator_t& {
		std::apply([&xs = this->xs_](auto&... idxs) { return xs.next_canonical(idxs...); }, ns_);
		++n_;
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator--() -> elements_iterator_t& {
		std::apply([&xs = this->xs_](auto&... idxs) { return xs.prev_canonical(idxs...); }, ns_);
		--n_;
		return *this;
	}

	BOOST_MULTI_HD constexpr auto operator+=(difference_type n) -> elements_iterator_t& {
		auto const nn = std::apply(xs_, ns_);
		ns_ = xs_.from_linear(nn + n);
		n_ += n;
		return *this;
	}
	BOOST_MULTI_HD constexpr auto operator-=(difference_type n) -> elements_iterator_t& {
		auto const nn = std::apply(xs_, ns_);
		ns_ = xs_.from_linear(nn - n);
		n_ -= n;
		return *this;
	}

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD /*[[gnu::pure]]*/ constexpr auto operator-(elements_iterator_t const& other) const -> difference_type {
		assert(base_ == other.base_ && l_ == other.l_);
		return n_ - other.n_;
	}
	BOOST_MULTI_HD constexpr auto operator<(elements_iterator_t const& other) const -> difference_type {
		assert(base_ == other.base_ && l_ == other.l_);
		return n_ < other.n_;
	}

	constexpr auto current() const -> pointer {return base_ + std::apply(l_, ns_);}

	BOOST_MULTI_HD constexpr auto operator->() const -> pointer   {return base_ + std::apply(l_, ns_) ;}
	BOOST_MULTI_HD constexpr auto operator*()  const -> reference {return base_  [std::apply(l_, ns_)];}
	BOOST_MULTI_HD constexpr auto operator[](difference_type const& n) const -> reference {
		auto const nn = std::apply(xs_, ns_);
		return base_[std::apply(l_, xs_.from_linear(nn + n))];
	}  // explicit here is necessary for nvcc/thrust

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	BOOST_MULTI_HD constexpr auto operator+(difference_type n) const -> elements_iterator_t {auto ret{*this}; ret += n; return ret;}  // explicitly necessary for nvcc/thrust
	BOOST_MULTI_HD constexpr auto operator-(difference_type n) const -> elements_iterator_t {auto ret{*this}; ret -= n; return ret;}  // explicitly necessary for nvcc/thrust

	BOOST_MULTI_HD constexpr auto operator==(elements_iterator_t const& other) const -> bool {
		assert(base_ == other.base_ && l_ == other.l_);  // TODO(correaa) calling host function from host device
		return n_ == other.n_;  // and base_ == other.base_ and l_ == other.l_;
	}
	BOOST_MULTI_HD constexpr auto operator!=(elements_iterator_t const& other) const -> bool {
		assert(base_ == other.base_ && l_ == other.l_);  // TODO(correaa) calling host function from host device
		return n_ != other.n_;
	}
};

template<typename Pointer, class LayoutType>
struct elements_range_t {
	using       pointer = Pointer;
	using layout_type = LayoutType;

	using value_type      = typename std::iterator_traits<pointer>::value_type;
	using const_pointer   = typename std::pointer_traits<pointer>::template rebind<value_type const>;

	using       reference = typename std::iterator_traits<      pointer>::reference;
	using const_reference = typename std::iterator_traits<const_pointer>::reference;

	using size_type       = typename std::iterator_traits<pointer>::difference_type;
	using difference_type = typename std::iterator_traits<pointer>::difference_type;

	using       iterator = elements_iterator_t<pointer, layout_type>;
	using const_iterator = elements_iterator_t<const_pointer, layout_type>;

	// using element        = value_type;

 private:
	pointer base_;
	layout_type l_;

 public:
	template<class OtherRange, decltype(multi::detail::implicit_cast<pointer>(std::declval<OtherRange>().base_))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*impl*/ elements_range_t(OtherRange const& other) : base_{other.base}, l_{other.l_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) to reproduce the implicitness of the argument
	template<class OtherRange, decltype(multi::detail::explicit_cast<pointer>(std::declval<OtherRange>().base_))* = nullptr>
	constexpr explicit elements_range_t(OtherRange const& other) : elements_range_t{other} {}

	constexpr elements_range_t(pointer base, layout_type const& lyt) : base_{base}, l_{lyt} {}

	constexpr auto base()       ->       pointer {return base_;}
	constexpr auto base() const -> const_pointer {return base_;}

	constexpr auto layout() const -> layout_type {return l_;}

 private:

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	constexpr auto at_aux_(difference_type n) const -> reference {
		assert( ! is_empty() );
		return base_[std::apply(l_, l_.extensions().from_linear(n))];
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	BOOST_MULTI_HD constexpr auto operator[](difference_type n) const& -> const_reference {return at_aux_(n);}
	BOOST_MULTI_HD constexpr auto operator[](difference_type n)     && ->       reference {return at_aux_(n);}
	BOOST_MULTI_HD constexpr auto operator[](difference_type n)      & ->       reference {return at_aux_(n);}

	constexpr auto size() const -> size_type {return l_.num_elements();}

	[[nodiscard]]
	constexpr auto    empty() const -> bool {return l_.   empty();}
	constexpr auto is_empty() const -> bool {return l_.is_empty();}

	elements_range_t(elements_range_t const&) = delete;
	elements_range_t(elements_range_t     &&) = delete;

	template<typename OP, class OL> auto operator==(elements_range_t<OP, OL> const& other) const -> bool {
		if( is_empty() && other.is_empty()) {return true;}
		return size() == other.size() &&     adl_equal(other.begin(), other.end(), begin());
	}
	template<typename OP, class OL> auto operator!=(elements_range_t<OP, OL> const& other) const -> bool {
		if(is_empty() && other.is_empty()) {return false;}
		return size() != other.size() ||  ! adl_equal(other.begin(), other.end(), begin());
	}

	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&  other)  & noexcept {assert(size() == other.size()); adl_swap_ranges(begin(), end(), other.begin());}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&  other) && noexcept {assert(size() == other.size()); adl_swap_ranges(begin(), end(), other.begin());}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&& other)  & noexcept {assert(size() == other.size()); adl_swap_ranges(begin(), end(), std::move(other).begin());}
	template<typename OP, class OL> void swap(elements_range_t<OP, OL>&& other) && noexcept {assert(size() == other.size()); adl_swap_ranges(begin(), end(), std::move(other).begin());}

	~elements_range_t() = default;

 private:
	constexpr auto begin_aux_() const {return iterator{base_, l_, 0                };}
	constexpr auto end_aux_  () const {return iterator{base_, l_, l_.num_elements()};}

 public:
	constexpr auto begin() const& -> const_iterator {return begin_aux_();}
	constexpr auto end  () const& -> const_iterator {return end_aux_  ();}

	constexpr auto begin()     && ->       iterator {return begin_aux_();}
	constexpr auto end  ()     && ->       iterator {return end_aux_  ();}

	constexpr auto begin()      & ->       iterator {return begin_aux_();}
	constexpr auto end  ()      & ->       iterator {return end_aux_  ();}

	constexpr auto front() const& -> const_reference {return *begin();}
	constexpr auto back () const& -> const_reference {return *std::prev(end(), 1);}

	constexpr auto front()     && ->       reference {return *begin();}
	constexpr auto back ()     && ->       reference {return *std::prev(end(), 1);}

	constexpr auto front()      & ->       reference {return *begin();}
	constexpr auto back ()      & ->       reference {return *std::prev(end(), 1);}

	auto operator=(elements_range_t const&) -> elements_range_t& = delete;

	auto operator=(elements_range_t     && other) noexcept -> elements_range_t& {  // cannot be =delete in NVCC?
		if(! is_empty()) {adl_copy(other.begin(), other.end(), this->begin());}
		return *this;
	}

	template<class OtherElementRange, class = decltype(adl_copy(std::begin(std::declval<OtherElementRange&&>()), std::end(std::declval<OtherElementRange&&>()), std::declval<iterator>()))>
	auto operator=(OtherElementRange&& other)  & -> elements_range_t& {assert(size() == other.size());  // NOLINT(cppcoreguidelines-missing-std-forward) std::forward<OtherElementRange>(other) creates a problem with move-only elements
		if(! is_empty()) {adl_copy(std::begin(other), std::end(other), begin());}
		return *this;
	}

	template<class OtherElementRange, class = decltype(adl_copy(std::begin(std::declval<OtherElementRange&&>()), std::end(std::declval<OtherElementRange&&>()), std::declval<iterator>()))>
	constexpr auto operator=(OtherElementRange&& other) && -> elements_range_t& {assert(size() == other.size());  // NOLINT(cppcoreguidelines-missing-std-forward) std::forward<OtherElementRange>(other) creates a problem with move-only elements
		if(! is_empty()) {adl_copy(std::begin(other), std::end(other), begin());} 
		return *this;
	}

	auto operator=(std::initializer_list<value_type> values) && -> elements_range_t& {operator=(values); return *this;}
	auto operator=(std::initializer_list<value_type> values) &  -> elements_range_t& {
		assert(static_cast<size_type>(values.size()) == size());
		adl_copy_n(values.begin(), values.size(), begin());
		return *this;
	}
};

template<class It>
[[deprecated("remove")]] BOOST_MULTI_HD constexpr auto ref(It begin, It end)
->multi::subarray<typename It::element, It::rank_v, typename It::element_ptr> {
	return multi::subarray<typename It::element, It::rank_v, typename It::element_ptr>{begin, end};
}

template<typename, ::boost::multi::dimensionality_type, class Alloc> struct static_array;  // this might be needed by MSVC 14.3 in c++17 mode

template<typename T, ::boost::multi::dimensionality_type D, typename ElementPtr, class Layout>
struct const_subarray : array_types<T, D, ElementPtr, Layout> {
	using types = array_types<T, D, ElementPtr, Layout>;
	using ref_ = const_subarray;

	using array_types<T, D, ElementPtr, Layout>::rank_v;

	friend struct const_subarray<typename types::element, D + 1, typename types::element_ptr >;

	using types::layout;
	using typename types::element_type;

	using layout_type = Layout;

	BOOST_MULTI_HD constexpr auto layout() const -> decltype(auto) {return array_types<T, D, ElementPtr, Layout>::layout();}

	using basic_const_array = subarray<T, D, typename std::pointer_traits<ElementPtr>::template rebind<element_type const>, Layout>;

	const_subarray() = default;
	auto operator=(const_subarray const&) -> const_subarray& = delete;
	auto operator=(const_subarray     &&) -> const_subarray& = delete;

	BOOST_MULTI_HD constexpr const_subarray(layout_type const& layout, ElementPtr const& base)
	: array_types<T, D, ElementPtr, Layout>{layout, base} {}

 protected:
	// using types::types;
	BOOST_MULTI_HD constexpr explicit const_subarray(std::nullptr_t nil) : types{nil} {}

	template<typename, ::boost::multi::dimensionality_type, class Alloc> friend struct static_array;

	// TODO(correaa) vvv consider making it explicit (seems that in C++23 it can prevent auto s = a[0];)
	const_subarray(const_subarray const&) = default;  // NOTE: reference type cannot be copied. perhaps you want to return by std::move or std::forward if you got the object from a universal reference argument

	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;

 public:
	using element           = typename types::element;
	using element_ptr       = typename types::element_ptr;
	using element_const_ptr = typename types::element_const_ptr;
	using element_ref       = typename types::element_ref;
	using element_cref      = typename std::iterator_traits<element_const_ptr>::reference;

	using  elements_iterator = elements_iterator_t<element_ptr      , layout_type>;
	using celements_iterator = elements_iterator_t<element_const_ptr, layout_type>;

	using       elements_range = elements_range_t<element_ptr      , layout_type>;
	using const_elements_range = elements_range_t<element_const_ptr, layout_type>;

	using index_gen [[deprecated("here to fulfill MultiArray concept")]] = char*;
	using extent_gen [[deprecated("here to fulfill MultiArray concept")]] = void ;
	using extent_range [[deprecated("here to fulfill MultiArray concept")]] = void;

 private:
	constexpr auto elements_aux_() const {return elements_range(this->base_, this->layout());}

 public:
	const_subarray(const_subarray&&) noexcept = default;  // lints(readability-redundant-access-specifiers)

	constexpr auto       elements()      & ->       elements_range { return elements_aux_(); }
	constexpr auto       elements()     && ->       elements_range { return elements_aux_(); }
	constexpr auto       elements() const&                         { return const_elements_range(this->base(), this->layout()); }
	constexpr auto const_elements() const  -> const_elements_range { return elements_aux_(); }

	constexpr auto hull() const -> std::pair<element_const_ptr, size_type> {
		return {this->base(), std::abs(this->hull_size())};
	}

	~const_subarray() = default;  // this lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	BOOST_MULTI_FRIEND_CONSTEXPR auto sizes(const_subarray const& self) noexcept -> typename const_subarray::sizes_type {return self.sizes();}  // needed by nvcc
	BOOST_MULTI_FRIEND_CONSTEXPR auto size (const_subarray const& self) noexcept -> typename const_subarray::size_type  {return self.size ();}  // needed by nvcc

	//  template<class T2> friend constexpr auto reinterpret_array_cast(const_subarray     && self) {return std::move(self).template reinterpret_array_cast<T2, typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>>();}
	//  template<class T2> friend constexpr auto reinterpret_array_cast(const_subarray const& self) {return           self .template reinterpret_array_cast<T2, typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>>();}

	friend constexpr auto dimensionality(const_subarray const& /*self*/) {return D;}

	using typename types::reference;

	using default_allocator_type = typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type;

	constexpr auto get_allocator() const -> default_allocator_type {
		using multi::get_allocator;
		return get_allocator(this->base());
	}

	BOOST_MULTI_FRIEND_CONSTEXPR auto get_allocator(const_subarray const& self) -> default_allocator_type {return self.get_allocator();}

	using decay_type = array<typename types::element_type, D, typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type>;

	friend constexpr auto decay(const_subarray const& self) -> decay_type {return self.decay();}
	       constexpr auto decay()           const&    -> decay_type {
		decay_type ret{*this};
		return ret;
	}

	constexpr auto operator+() const -> decay_type {return decay();}
	using typename types::const_reference;

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	BOOST_MULTI_HD constexpr auto at_aux_(index idx) const {
		#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
		#endif
		return const_reference(
			this->layout().sub(),
			this->base_ + (idx*this->layout().stride() - this->layout().offset())
		);  // cppcheck-suppress syntaxError ; bug in cppcheck 2.5
		#if defined(__clang__)
		#pragma clang diagnostic pop
		#endif
	}

 public:

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD constexpr auto operator[](index idx) const& -> const_reference {
		return const_reference(
			this->layout().sub(),
			this->base_ + (idx*this->layout().stride() - this->layout().offset())
		);  // cppcheck-suppress syntaxError ; bug in cppcheck 2.5
		// return at_aux_(idx);
	}  // TODO(correaa) use return type to cast

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	//  BOOST_MULTI_HD constexpr auto operator[](index idx)     && ->     reference { return                              at_aux_(idx) ; }
	//  BOOST_MULTI_HD constexpr auto operator[](index idx)      & ->     reference { return                              at_aux_(idx) ; }

	template<class Tuple = std::array<index, static_cast<std::size_t>(D)>, 
		typename = std::enable_if_t<(std::tuple_size<Tuple>::value > 1)>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	BOOST_MULTI_HD constexpr auto operator[](Tuple const& tup) const
	->decltype(operator[](detail::head(tup))[detail::tuple_tail(tup)]) {
		return operator[](detail::head(tup))[detail::tuple_tail(tup)]; }

	template<class Tuple, typename = std::enable_if_t<(std::tuple_size<Tuple>::value == 1)> >  // NOLINT(modernize-use-constraints) TODO(correaa)
	BOOST_MULTI_HD constexpr auto operator[](Tuple const& tup) const
	->decltype(operator[](detail::head(tup))) {
		return operator[](detail::head(tup)); }

	constexpr auto front() const& -> const_reference { return *begin(); }
	constexpr auto back()  const& -> const_reference { return *(end() - 1); }  // std::prev(end(), 1);}

	constexpr auto front()     && ->       reference { return *begin(); }
	constexpr auto back()      && ->       reference { return *(end() - 1); }  // std::prev(end(), 1);}

	constexpr auto front()      & ->       reference {return *begin();}
	constexpr auto back()       & ->       reference {return *(end() - 1);}  // std::prev(end(), 1);}

	using typename types::index;

	constexpr auto reindexed(index first) const& {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return const_subarray(new_layout, types::base_);
	}
	constexpr auto reindexed(index first)& {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return const_subarray(new_layout, types::base_);
	}
	constexpr auto reindexed(index first)&& -> const_subarray {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}

	// TODO(correaa) : implement reindexed_aux
	template<class... Indexes>
	constexpr auto reindexed(index first, Indexes... idxs) const& -> const_subarray {
		return ((reindexed(first).rotated()).reindexed(idxs...)).unrotated();
	}

 private:
	constexpr auto taked_aux_(difference_type n) const {
		assert( n <= this->size() );
		typename types::layout_t const new_layout(
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride()*n
		);
		return const_subarray(new_layout, this->base_);
	}

 public:
	constexpr auto taked(difference_type n) const& -> basic_const_array { return taked_aux_(n); }
	// constexpr auto taked(difference_type n)     && -> const_subarray    { return taked_aux_(n); }
	// constexpr auto taked(difference_type n)      & -> const_subarray    { return taked_aux_(n); }

 private:
	constexpr auto dropped_aux_(difference_type n) const {
		assert( n <= this->size() );
		typename types::layout_t const new_layout{
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride()*(this->size() - n)
		};

		#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
		#endif

		return const_subarray(new_layout, this->base_ + n*this->layout().stride() - this->layout().offset());

		#if defined(__clang__)
		#pragma clang diagnostic pop
		#endif
	}

 public:
	constexpr auto dropped(difference_type n) const& -> basic_const_array { return dropped_aux_(n); }
	constexpr auto dropped(difference_type n)     && ->    const_subarray { return dropped_aux_(n); }
	constexpr auto dropped(difference_type n)      & ->    const_subarray { return dropped_aux_(n); }

 private:

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD constexpr auto sliced_aux_(index first, index last) const {
		// TODO(correaa) remove first == last condition
		BOOST_MULTI_ACCESS_ASSERT(((first==last) || this->extension().contains(first   ))&&"sliced first out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		BOOST_MULTI_ACCESS_ASSERT(((first==last) || this->extension().contains(last - 1))&&"sliced last  out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		typename types::layout_t new_layout = this->layout();
		new_layout.nelems() = this->stride()*(last - first);  // TODO(correaa) : reconstruct layout instead of mutating it
		BOOST_MULTI_ACCESS_ASSERT(this->base_ || ((first*this->layout().stride() - this->layout().offset()) == 0) );  // it is UB to offset a nullptr
		return const_subarray{new_layout, this->base_ + (first*this->layout().stride() - this->layout().offset())};
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) const& -> const_subarray { return sliced_aux_(first, last); }

	constexpr auto blocked(index first, index last) const& -> basic_const_array { return sliced(first, last).reindexed(first); }
	constexpr auto blocked(index first, index last)      & -> const_subarray { return sliced(first, last).reindexed(first); }

	using iextension = typename const_subarray::index_extension;

	constexpr auto stenciled(iextension iex)                                                         & -> const_subarray { return blocked(iex.first(), iex.last()); }
	constexpr auto stenciled(iextension iex, iextension iex1)                                        & -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2)                       & -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1, iex2)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3)      & -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3)).unrotated(); }
	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs)     & -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3, iexs...)).unrotated(); }

	constexpr auto stenciled(iextension iex)                                                        && -> const_subarray { return blocked(iex.first(), iex.last()); }
	constexpr auto stenciled(iextension iex, iextension iex1)                                       && -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2)                      && -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1, iex2)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3)     && -> const_subarray { return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3)).unrotated(); }
	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs)    && -> const_subarray{ return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3, iexs...)).unrotated(); }

	constexpr auto stenciled(iextension iex)                                                    const& -> basic_const_array { return blocked(iex.first(), iex.last()); }
	constexpr auto stenciled(iextension iex, iextension iex1)                                   const& -> basic_const_array { return ((stenciled(iex).rotated()).stenciled(iex1)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2)                  const& -> basic_const_array { return ((stenciled(iex).rotated()).stenciled(iex1, iex2)).unrotated(); }
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3) const& -> basic_const_array { return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3)).unrotated(); }

	template<class... Xs>
	constexpr auto stenciled(iextension iex, iextension iex1, iextension iex2, iextension iex3, Xs... iexs) const& -> basic_const_array {
		return ((stenciled(iex).rotated()).stenciled(iex1, iex2, iex3, iexs...)).unrotated();
	}

	constexpr auto elements_at(size_type idx) const& -> decltype(auto) {
		assert(idx < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](idx / sub_num_elements).elements_at(idx % sub_num_elements);
	}
	constexpr auto elements_at(size_type idx) && -> decltype(auto) {
		assert(idx < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](idx / sub_num_elements).elements_at(idx % sub_num_elements);
	}
	constexpr auto elements_at(size_type idx) & -> decltype(auto) {
		assert(idx < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](idx / sub_num_elements).elements_at(idx % sub_num_elements);
	}

 private:
	constexpr auto strided_aux_(difference_type diff) const {
		typename types::layout_t const new_layout{this->layout().sub(), this->layout().stride()*diff, this->layout().offset(), this->layout().nelems()};
		return const_subarray(new_layout, types::base_);
	}

 public:
	constexpr auto strided(difference_type diff) const& -> basic_const_array { return strided_aux_(diff); }
	constexpr auto strided(difference_type diff)     && ->    const_subarray { return strided_aux_(diff); }
	constexpr auto strided(difference_type diff)      & ->    const_subarray { return strided_aux_(diff); }

	constexpr auto sliced(
		typename types::index first, typename types::index last, typename types::index stride_
	) const& -> const_subarray {
		return sliced(first, last).strided(stride_);
	}

	using index_range = typename const_subarray::index_range;

	BOOST_MULTI_HD constexpr auto range(index_range irng) const& -> decltype(auto) {return                  sliced(irng.front(), irng.front() + irng.size());}
	// constexpr auto range(index_range irng)     && -> decltype(auto) {return std::move(*this).sliced(irng.front(), irng.front() + irng.size());}
	// constexpr auto range(index_range irng)      & -> decltype(auto) {return                  sliced(irng.front(), irng.front() + irng.size());}

	[[deprecated("is_flattable will be a property of the layout soon")]]
	constexpr auto is_flattable() const -> bool{
		return
			   (this->size() <= 1)
			|| (this->stride() == this->layout().sub().nelems())
		;
	}

	// friend constexpr auto flatted(const_subarray const& self) {return self.flatted();}
	constexpr auto flatted()              const& {
		// assert(is_flattable() && "flatted doesn't work for all layouts!");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		multi::layout_t<D-1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return const_subarray<T, D-1, ElementPtr>{new_layout, types::base_};
	}

	void flattened() const = delete;
	// {
	//  multi::biiterator<std::decay_t<decltype(this->begin())>> biit{this->begin(), 0, size(*(this->begin()))};
	//  return basic_array<T, D-1, decltype(biit)>(this->layout().sub, biit);
	// }

	constexpr auto broadcasted() const& {
		multi::layout_t<D + 1> const new_layout{layout(), 0, 0, (std::numeric_limits<size_type>::max)()};  // paren for MSVC macros
		return const_subarray<T, D+1, typename const_subarray::element_const_ptr>{new_layout, types::base_};
	}

	// TODO(correaa) : define a diagonal_aux
	constexpr auto diagonal()    && {return this->diagonal();}

	constexpr auto diagonal()     & -> const_subarray<T, D-1, typename const_subarray::element_ptr> {
		using boost::multi::detail::get;
		auto square_size = (std::min)(get<0>(this->sizes()), get<1>(this->sizes()));  // paren for MSVC macros
		multi::layout_t<D-1> new_layout{(*this)({0, square_size}, {0, square_size}).layout().sub()};
		new_layout.nelems() += (*this)({0, square_size}, {0, square_size}).layout().nelems();  // TODO(correaa) : don't use mutation
		new_layout.stride() += (*this)({0, square_size}, {0, square_size}).layout().stride();  // TODO(correaa) : don't use mutation
		return {new_layout, types::base_};
	}

	template<class Dummy = void, std::enable_if_t<(D > 1) && sizeof(Dummy*), int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	constexpr auto diagonal() const& -> const_subarray<T, D-1, typename const_subarray::element_const_ptr> {
		using std::get;
		auto const square_size = (std::min)(get<0>(this->sizes()), get<1>(this->sizes()));  // parenthesis min for MSVC macros
		multi::layout_t<D-1> new_layout{(*this)({0, square_size}, {0, square_size}).layout().sub()};
		new_layout.nelems() += (*this)({0, square_size}, {0, square_size}).layout().nelems();
		new_layout.stride() += (*this)({0, square_size}, {0, square_size}).layout().stride();  // cppcheck-suppress arithOperationsOnVoidPointer ; false positive D == 1 doesn't happen here
		return {new_layout, types::base_};
	}

	friend constexpr auto diagonal(const_subarray const& self) {return           self .diagonal();}
	friend constexpr auto diagonal(const_subarray&       self) {return           self .diagonal();}
	friend constexpr auto diagonal(const_subarray&&      self) {return std::move(self).diagonal();}

	// using partitioned_type       = const_subarray<T, D+1, element_ptr      >;
	// using partitioned_const_type = const_subarray<T, D+1, element_const_ptr>;

 private:
	BOOST_MULTI_HD constexpr auto partitioned_aux_(size_type n) const {
		assert(n != 0);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		// vvv TODO(correaa) should be size() here?
		// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) normal in a constexpr function
		assert( (this->layout().nelems() % n) == 0);  // if you get an assertion here it means that you are partitioning an array with an incommunsurate partition
		multi::layout_t<D+1> new_layout{this->layout(), this->layout().nelems()/n, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= n;
		return subarray<T, D+1, element_ptr>(new_layout, types::base_);
	}

 public:
	BOOST_MULTI_HD constexpr auto partitioned(size_type n) const& -> const_subarray<T, D+1, element_ptr> {return partitioned_aux_(n);}

	// friend BOOST_MULTI_HD constexpr auto partitioned(const_subarray const& self, size_type n) -> partitioned_const_type {return           self .partitioned(n);}
	// friend BOOST_MULTI_HD constexpr auto partitioned(const_subarray      & self, size_type n) -> partitioned_type       {return           self .partitioned(n);}
	// friend BOOST_MULTI_HD constexpr auto partitioned(const_subarray     && self, size_type n) -> partitioned_type       {return std::move(self).partitioned(n);}

 private:
	BOOST_MULTI_HD constexpr auto chunked_aux_(size_type count) const {
		assert( this->size() % count == 0 );
		return partitioned_aux_(this->size()/count);
	}

 public:  // in Mathematica this is called Partition https://reference.wolfram.com/language/ref/Partition.html in RangesV3 it is called chunk
	BOOST_MULTI_HD constexpr auto chunked(size_type count) const& -> const_subarray<T, D+1, element_ptr> {return chunked_aux_(count);}
	//  BOOST_MULTI_HD constexpr auto chunked(size_type count)      & -> partitioned_type       {return chunked_aux_(count);}
	//  BOOST_MULTI_HD constexpr auto chunked(size_type count)     && -> partitioned_type       {return chunked_aux_(count);}

	constexpr auto tiled(size_type count) const & {
		assert(count != 0);
		struct divided_type {
			const_subarray<T, D+1, element_ptr> quotient;
			const_subarray<T, D, element_ptr> remainder;
		};
		return divided_type{
			this->taked(this->size() - (this->size() % count)).chunked(count),
			this->dropped(this->size() - (this->size() % count))
		};
	}

 private:
	constexpr auto reversed_aux_() const -> const_subarray {
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}

 public:
	       constexpr auto reversed()           const&    -> basic_const_array { return reversed_aux_(); }
	       constexpr auto reversed()                &    ->          const_subarray { return reversed_aux_(); }
	       constexpr auto reversed()               &&    ->          const_subarray { return reversed_aux_(); }
	friend constexpr auto reversed(const_subarray const& self) -> basic_const_array { return           self .reversed(); }
	friend constexpr auto reversed(const_subarray      & self) ->          const_subarray { return           self .reversed(); }
	friend constexpr auto reversed(const_subarray     && self) ->          const_subarray { return std::move(self).reversed(); }

 private:
	BOOST_MULTI_HD constexpr auto transposed_aux_() const {
		auto new_layout = this->layout();
		new_layout.transpose();
		return const_subarray(new_layout, types::base_);
	}

 public:
	BOOST_MULTI_HD constexpr auto transposed() const& -> const_subarray { return transposed_aux_(); }
	// BOOST_MULTI_HD constexpr auto transposed()      & -> const_subarray { return transposed_aux_(); }
	// BOOST_MULTI_HD constexpr auto transposed()     && -> const_subarray { return transposed_aux_(); }

	// friend BOOST_MULTI_HD /*constexpr*/ auto transposed(const_subarray const& self) -> basic_const_array { return           self .transposed(); }
	// friend BOOST_MULTI_HD /*constexpr*/ auto transposed(const_subarray      & self) ->    const_subarray { return           self .transposed(); }
	// friend BOOST_MULTI_HD /*constexpr*/ auto transposed(const_subarray     && self) ->    const_subarray { return std::move(self).transposed(); }

	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	auto operator~ (const_subarray const& self) -> const_subarray {return self.transposed();}
	// BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	// auto operator~ (const_subarray& self) -> const_subarray {return self.transposed();}
	// BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	// auto operator~ (const_subarray&& self) -> const_subarray {return std::move(self).transposed();}

 private:
	BOOST_MULTI_HD constexpr auto rotated_aux_() const {
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate();
		return const_subarray(new_layout, types::base_);
	}
	BOOST_MULTI_HD constexpr auto unrotated_aux_() const {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return const_subarray(new_layout, types::base_);
	}

 public:
	BOOST_MULTI_HD constexpr auto   rotated() const& -> const_subarray { return   rotated_aux_(); }
	BOOST_MULTI_HD constexpr auto unrotated() const& -> const_subarray { return unrotated_aux_(); }

	// BOOST_MULTI_FRIEND_CONSTEXPR auto unrotated(const_subarray const& self) { return           self .unrotated(); }
	// BOOST_MULTI_FRIEND_CONSTEXPR auto unrotated(const_subarray      & self) { return           self .unrotated(); }
	// BOOST_MULTI_FRIEND_CONSTEXPR auto unrotated(const_subarray     && self) { return std::move(self).unrotated(); }

	// constexpr auto operator|(typename const_subarray::size_type n)      & -> decltype(auto) { return                  partitioned(n); }
	// constexpr auto operator|(typename const_subarray::size_type n)     && -> decltype(auto) { return std::move(*this).partitioned(n); }
	// constexpr auto operator|(typename const_subarray::size_type n) const& -> decltype(auto) { return                  partitioned(n); }

 private:
	template<typename, ::boost::multi::dimensionality_type, typename, class> friend struct const_subarray;

	BOOST_MULTI_HD constexpr auto paren_aux_() const& {return const_subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_);}

 public:
	BOOST_MULTI_HD constexpr auto operator()() const& -> const_subarray {return paren_aux_();}  // NOLINT(readability-redundant-access-specifiers,readability-const-return-type)

 private:
	template<class... As>  BOOST_MULTI_HD  constexpr auto paren_aux_(index_range rng, As... args) const& {return range(rng).rotated().paren_aux_(args...).unrotated();}

	// template<class... As>    constexpr auto paren_aux_(intersecting_range<index> inr, As... args)      & -> decltype(auto) {return paren_aux_(intersection(this->extension(), inr), args...);}
	// template<class... As>    constexpr auto paren_aux_(intersecting_range<index> inr, As... args)     && -> decltype(auto) {return paren_aux_(intersection(this->extension(), inr), args...);}
	template<class... As> BOOST_MULTI_HD   constexpr auto paren_aux_(intersecting_range<index> inr, As... args) const& -> decltype(auto) {return paren_aux_(intersection(this->extension(), inr), args...);}

	// template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args)      & -> decltype(auto) {return operator[](idx).paren_aux_(args...);}
	// template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args)     && -> decltype(auto) {return operator[](idx).paren_aux_(args...);}
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args) const& -> decltype(auto) {return operator[](idx).paren_aux_(args...);}

 public:
	// vvv DO NOT remove default parameter `= irange` : the default template parameters below help interpret the expression `{first, last}` syntax as index ranges
	template<class A1 = irange>                                                                        BOOST_MULTI_HD constexpr auto operator()(A1 arg1)                                        const& -> decltype(auto) {return paren_aux_(arg1);}  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange>                                                     BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2)                               const& -> decltype(auto) {return paren_aux_(arg1, arg2);}  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange>                                  BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3)                      const& -> decltype(auto) {return paren_aux_(arg1, arg2, arg3);}  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As>  BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) const& -> decltype(auto) {return paren_aux_(arg1, arg2, arg3, arg4, args...);}

 private:
	template<typename Tuple, std::size_t ... I> BOOST_MULTI_HD constexpr auto apply_impl_(Tuple const& tuple, std::index_sequence<I...>/*012*/) const& -> decltype(auto) {return this->operator()(std::get<I>(tuple)...);}

 public:
	template<typename Tuple> constexpr auto apply(Tuple const& tuple) const& -> decltype(auto) {return apply_impl_(tuple, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	// template<typename Tuple> constexpr auto apply(Tuple const& tuple)     && -> decltype(auto) {return apply_impl_(tuple, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	// template<typename Tuple> constexpr auto apply(Tuple const& tuple)      & -> decltype(auto) {return apply_impl_(tuple, std::make_index_sequence<std::tuple_size<Tuple>::value>());}

	using       iterator = array_iterator<element, D, element_ptr      >;
	using const_iterator = array_iterator<element, D, element_ptr, true>;
	using move_iterator  = array_iterator<element, D, element_ptr, false, true>;

	// using  move_iterator = array_iterator<element, D, element_move_ptr >;

	// using       reverse_iterator [[deprecated]] = std::reverse_iterator<      iterator>;
	// using const_reverse_iterator [[deprecated]] = std::reverse_iterator<const_iterator>;

	const_subarray(iterator first, iterator last)
	: const_subarray(layout_type(first->layout(), first.stride(), 0, (last - first)*first->size()), first.base()) {
		assert(first->layout() == last->layout());
	}

 private:
	// [[deprecated("remove")]] BOOST_MULTI_HD constexpr explicit const_subarray(iterator begin, iterator end)
	// : const_subarray(
	//  layout_type{begin->layout(), begin.stride(), 0, begin.stride() * (end - begin)},
	//  begin.base()
	// ) {
	//  assert(begin.stride() == end.stride());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  assert(begin->layout() == end->layout());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	// }

	friend BOOST_MULTI_HD constexpr auto ref<iterator>(iterator begin, iterator end) -> multi::subarray<typename iterator::element, iterator::rank_v, typename iterator::element_ptr>;

 public:
	using       ptr =       subarray_ptr<T, D, ElementPtr, Layout>;
	using const_ptr = const_subarray_ptr<T, D, ElementPtr, Layout>;  // TODO(correaa) add const_subarray_ptr

	using pointer = ptr;
	using const_pointer = const_ptr;

 private:
	constexpr auto addressof_aux_() const {return ptr(this->base_, this->layout());}

 public:
	constexpr auto addressof()     && ->       ptr { return addressof_aux_(); }
	constexpr auto addressof()      & ->       ptr { return addressof_aux_(); }
	constexpr auto addressof() const& -> const_ptr { return addressof_aux_(); }

	// NOLINTBEGIN(google-runtime-operator) //NOSONAR
	// operator& is not defined for r-values anyway
	constexpr auto operator&()     && { return addressof(); }  // NOLINT(runtime/operator) //NOSONAR
	// [[deprecated("controversial")]]
	constexpr auto operator&()      & { return addressof(); }  // NOLINT(runtime/operator) //NOSONAR
	// [[deprecated("controversial")]]
	constexpr auto operator&() const& { return addressof(); }  // NOLINT(runtime/operator) //NOSONAR
	// NOLINTEND(google-runtime-operator)
        
 private:
	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD constexpr auto begin_aux_() const {return iterator(types::base_                 , this->sub(), this->stride());}
	               constexpr auto end_aux_  () const {return iterator(types::base_ + this->nelems(), this->sub(), this->stride());}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	       BOOST_MULTI_HD constexpr     auto begin()       &       {return begin_aux_();}
	                      constexpr     auto end  ()       &       {return end_aux_()  ;}
	friend BOOST_MULTI_HD /*constexpr*/ auto begin(const_subarray& self) {return self.begin();}
	friend constexpr                    auto end  (const_subarray& self) {return self.end  ();}

	       constexpr     auto begin()       &&       {return begin();}
	       constexpr     auto end  ()       &&       {return end()  ;}
	friend /*constexpr*/ auto begin(const_subarray&& self) {return std::move(self).begin();}
	friend /*constexpr*/ auto end  (const_subarray&& self) {return std::move(self).end()  ;}

	       constexpr     auto  begin()        const&       -> const_iterator { return begin_aux_(); }
	       constexpr     auto  end  ()        const&       -> const_iterator { return end_aux_()  ; }
	friend /*constexpr*/ auto  begin(const_subarray const& self) -> const_iterator { return self.begin(); }  // NOLINT(whitespace/indent) constexpr doesn't work with nvcc friend
	friend /*constexpr*/ auto  end  (const_subarray const& self) -> const_iterator { return self.end()  ; }  // NOLINT(whitespace/indent) constexpr doesn't work with nvcc friend

	BOOST_MULTI_HD     constexpr auto cbegin()           const& {return begin();}
	/*fd*/ constexpr auto cend()             const& {return end()  ;}
	friend constexpr auto cbegin(const_subarray const& self) {return self.cbegin();}
	friend constexpr auto cend  (const_subarray const& self) {return self.cend()  ;}

	using       cursor = cursor_t<typename const_subarray::element_ptr      , D, typename const_subarray::strides_type>;
	using const_cursor = cursor_t<typename const_subarray::element_const_ptr, D, typename const_subarray::strides_type>;
 private:
	constexpr auto home_aux_() const {return cursor(this->base_, this->strides());}

 public:
	constexpr auto home() const& -> const_cursor { return home_aux_(); }

	// const_subarray(cursor home, typename const_subarray::sizes_type szs) {}

	template<class It> constexpr auto assign(It first) & -> It { adl_copy_n(first, this->size(), begin()); std::advance(first, this->size()); return first; }
	template<class It> constexpr auto assign(It first)&& -> It { return assign(first);}

	template<
		class Range,
		std::enable_if_t<!has_extensions<std::decay_t<Range>>::value, int> = 0,
		//  std::enable_if_t<not multi::is_implicitly_convertible_v<subarray, Range>, int> =0,
		class = decltype(Range(std::declval<typename const_subarray::const_iterator>(), std::declval<typename const_subarray::const_iterator>()))>
	constexpr explicit operator Range() const { return Range(begin(), end()); }  // NOLINT(fuchsia-default-arguments-calls) for example std::vector(it, ti, alloc = {})

	template<class TT, class... As>
	friend constexpr auto operator==(const_subarray const& self, const_subarray<TT, D, As...> const& other) -> bool {
		return (self.extension() == other.extension()) && (self.elements() == other.elements());
	}
	template<class TT, class... As>
	friend constexpr auto operator!=(const_subarray const& self, const_subarray<TT, D, As...> const& other) -> bool {
		return (self.extension() != other.extension()) ||  (self.elements() != other.elements());
	}

	constexpr auto operator==(const_subarray const& other) const -> bool {
		return (this->extension() == other.extension()) && (this->elements() == other.elements());
	}
	constexpr auto operator!=(const_subarray const& other) const -> bool {
		return (this->extension() != other.extension()) ||  (this->elements() != other.elements());
	}

	friend constexpr auto lexicographical_compare(const_subarray const& self, const_subarray const& other) -> bool {
		if(self.extension().first() > other.extension().first()) {return true ;}
		if(self.extension().first() < other.extension().first()) {return false;}
		return adl_lexicographical_compare(
			self.begin(), self.end(),
			other.begin(), other.end()
		);
	}

	constexpr auto operator< (const_subarray const& other) const& -> bool {return lexicographical_compare(*this, other);}
	constexpr auto operator<=(const_subarray const& other) const& -> bool {return *this == other || lexicographical_compare(*this, other);}
	constexpr auto operator> (const_subarray const& other) const& -> bool {return other < *this;}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		std::enable_if_t<  std::is_const_v<typename std::pointer_traits<P2>::element_type>,int> =0  // NOLINT(modernize-use-constraints) TODO(correaa)
	>
	constexpr auto static_array_cast() const & {  // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), static_cast<P2>(this->base_));  // TODO(correaa) might violate constness
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		std::enable_if_t<! std::is_const_v<typename std::pointer_traits<P2>::element_type>,int> =0  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	[[deprecated("violates constness")]]
	constexpr auto static_array_cast() const & {  // name taken from std::static_pointer_cast
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
	constexpr auto static_array_cast_(Args&&... args) const & {  // name taken from std::static_pointer_cast
		return subarray<T2, D, P2>(this->layout(), P2{this->base_, std::forward<Args>(args)...});
	}

 public:
	template<class UF>
	constexpr auto element_transformed(UF&& fun) const& {
		return static_array_cast_<
		//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
			std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
			transform_ptr<
			//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
				std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
				UF, element_const_ptr, std::invoke_result_t<UF const&, element_cref>
			>
		>(std::forward<UF>(fun));
	}
	template<class UF>
	constexpr auto element_transformed(UF&& fun)  & {
		return static_array_cast_<
			std::decay_t<std::invoke_result_t<UF const&, element_ref >>,
			transform_ptr<
				std::decay_t<std::invoke_result_t<UF const&, element_ref >>,
				UF, element_ptr      , std::invoke_result_t<UF const&, element_ref >
			>
		>(std::forward<UF>(fun));
	}
	template<class UF>
	constexpr auto element_transformed(UF&& fun) && {return element_transformed(std::forward<UF>(fun));}

	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2 const>,
		class Element = typename const_subarray::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM member) const& -> subarray<T2, D, P2> {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
			"Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements."
		);

		return subarray<T2, D, P2>{this->layout().scale(sizeof(T), sizeof(T2)), static_cast<P2>(&(this->base_->*member))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM member) & -> subarray<T2, D, P2> {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
			"Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements"
		);

		return subarray<T2, D, P2>{this->layout().scale(sizeof(T), sizeof(T2)), static_cast<P2>(&(this->base_->*member))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM member) && -> subarray<T2, D, P2> {
		return this->member_cast<T2, P2, Element, PM>(member);
	}

	template<class T2, class P2 = typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>>
	using rebind = subarray<std::decay_t<T2>, D, P2>;

	template<
		class T2 = std::remove_const_t<T>,
		class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		std::enable_if_t<  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
			std::is_same_v<  // check that pointer family is not changed
				typename std::pointer_traits<P2>::template rebind<T2>,
				typename std::pointer_traits<element_ptr>::template rebind<T2>
			>
			&& 
			std::is_same_v<  // check that only constness is changed
				std::remove_const_t<typename std::pointer_traits<P2>::element_type>,
				std::remove_const_t<typename const_subarray::element_type>
			>
		, int> =0
	>
	constexpr auto const_array_cast() const {
		if constexpr(std::is_pointer_v<P2>) {
			return rebind<T2, P2>(this->layout(), const_cast      <P2       >(this->base_));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		} else {
			return rebind<T2, P2>(this->layout(), reinterpret_cast<P2 const&>(this->base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)  //NOSONAR
		}
	}

	constexpr auto as_const() const {
		return rebind<element, element_const_ptr>{this->layout(), this->base()};
	}

 private:
	template<class T2, class P2>
	constexpr auto reinterpret_array_cast_aux_() const -> rebind<T2, P2> {
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return {
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)  // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		};
	}

 public:
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2 const>>
	constexpr auto reinterpret_array_cast() const& {return reinterpret_array_cast_aux_<T2, P2>().as_const();}

	template<
		class T2,
		class P2 =
			std::conditional_t<
				std::is_const_v<typename std::pointer_traits<typename const_subarray::element_ptr>::element_type>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2 const>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<T2>
			>
	>
	constexpr auto reinterpret_array_cast(size_type count) const& {
		static_assert(sizeof(T) % sizeof(T2) == 0, "error: reinterpret_array_cast is limited to integral stride values");

		assert(sizeof(T) == sizeof(T2) * static_cast<std::size_t>(count));  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : checck implicit size compatibility

		if constexpr(std::is_pointer_v<ElementPtr>) {
			using void_ptr_like = std::conditional_t<
				std::is_const_v<typename std::pointer_traits<decltype(this->base_)>::element_type>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<void const>,
				typename std::pointer_traits<typename const_subarray::element_ptr>::template rebind<void>
			>;
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

	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int version) {
		using AT = multi::archive_traits<Archive>;
		if(version == 0) {
			std::for_each(this->begin(), this->end(), [&](reference&& item) {arxiv & AT    ::make_nvp("item", std::move(item));});
		} else {
			std::for_each(this->elements().begin(), this->elements().end(), [&](element& elem) {arxiv & AT    ::make_nvp("elem", elem);});
		}
	//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv & cereal::make_nvp("item", item);});
	//  std::for_each(this->begin(), this->end(), [&](auto&& item) {arxiv &                          item ;});
	}
};

template<class T>
BOOST_MULTI_HD constexpr auto move(T&& val) noexcept -> decltype(auto) {
	if constexpr(has_member_move<T>::value) {
		return std::forward<T>(val).move();
	} else {
		return std::move(std::forward<T>(val));
	}
}

template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout>
class move_subarray : public subarray<T, D, ElementPtr, Layout> {
	BOOST_MULTI_HD constexpr move_subarray(subarray<T, D, ElementPtr, Layout>& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  TODO(correa) check if this is necessary
	: subarray<T, D, ElementPtr, Layout>(other.layout(), other.mutable_base()) {}

	friend class subarray<T, D, ElementPtr, Layout>;

 public:
	using subarray<T, D, ElementPtr, Layout>::operator[];
	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> decltype(auto) {
		return multi::move(subarray<T, D, ElementPtr, Layout>::operator[](idx));
	}

	using subarray<T, D, ElementPtr, Layout>::begin;
	using subarray<T, D, ElementPtr, Layout>::end;

	auto begin() && { return this->mbegin(); }
	auto end  () && { return this->mend  (); }
};

template<typename T, multi::dimensionality_type D, typename ElementPtr, class Layout>
class subarray : public const_subarray<T, D, ElementPtr, Layout> {
	BOOST_MULTI_HD constexpr subarray(const_subarray<T, D, ElementPtr, Layout> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  TODO(correa) check if this is necessary
	: subarray(other.layout(), other.mutable_base()) {}

	template<typename, multi::dimensionality_type, typename, class> friend class subarray;
	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;

	template<class, multi::dimensionality_type, class, bool, bool> friend struct array_iterator;

	subarray(subarray const&) = default;

 public:
	BOOST_MULTI_HD constexpr auto move() {return move_subarray<T, D, ElementPtr, Layout>(*this);}
	friend BOOST_MULTI_HD constexpr auto move(subarray& self) { return self.move(); }
	friend BOOST_MULTI_HD constexpr auto move(subarray&& self) { return std::move(self).move(); }

	using move_iterator  = array_iterator<T, D, ElementPtr, false, true>;

	subarray(subarray&&) noexcept = default;
	~subarray() = default;

	using ptr = subarray_ptr<T, D, ElementPtr, Layout>;

	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&()     && { return subarray_ptr<T, D, ElementPtr, Layout>(this->base_, this->layout()); }  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR
	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&()      & { return subarray_ptr<T, D, ElementPtr, Layout>(this->base_, this->layout()); } // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR

	using const_subarray<T, D, ElementPtr, Layout>::operator&;
	// NOLINTNEXTLINE(runtime/operator)
	// BOOST_MULTI_HD constexpr auto operator&() const& {return subarray_ptr<const_subarray, Layout>{this->base_, this->layout()};}  // NOLINT(google-runtime-operator) extend semantics  //NOSONAR

	using const_subarray<T, D, ElementPtr, Layout>::const_subarray;

	constexpr auto mbegin() { return move_iterator{this->begin()}; }
	constexpr auto mend()   { return move_iterator{this->end()  }; }

	using const_subarray<T, D, ElementPtr, Layout>::home;
	constexpr auto home()     && { return this->home_aux_(); }
	constexpr auto home()      & { return this->home_aux_(); }

	using const_subarray<T, D, ElementPtr, Layout>::strided;
	constexpr auto strided(difference_type diff) && -> subarray { return this->strided_aux_(diff);}
	constexpr auto strided(difference_type diff)  & -> subarray { return this->strided_aux_(diff);}

	using const_subarray<T, D, ElementPtr, Layout>::taked;
	constexpr auto taked(difference_type count)  && -> subarray {return this->taked_aux_(count);}
	constexpr auto taked(difference_type count)   & -> subarray {return this->taked_aux_(count);}

	using const_subarray<T, D, ElementPtr, Layout>::dropped;
	constexpr auto dropped(difference_type count) && -> subarray { return this->dropped_aux_(count); }
	constexpr auto dropped(difference_type count)  & -> subarray { return this->dropped_aux_(count); }

	using const_subarray<T, D, ElementPtr, Layout>::rotated;
	BOOST_MULTI_HD constexpr auto rotated() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::rotated(); }
	BOOST_MULTI_HD constexpr auto rotated() &  -> subarray { return const_subarray<T, D, ElementPtr, Layout>::rotated(); }

	using const_subarray<T, D, ElementPtr, Layout>::unrotated;
	BOOST_MULTI_HD constexpr auto unrotated() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unrotated(); }
	BOOST_MULTI_HD constexpr auto unrotated() &  -> subarray { return const_subarray<T, D, ElementPtr, Layout>::unrotated(); }

	using const_subarray<T, D, ElementPtr, Layout>::transposed;
	BOOST_MULTI_HD constexpr auto transposed() && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::transposed(); }
	BOOST_MULTI_HD constexpr auto transposed() &  -> subarray { return const_subarray<T, D, ElementPtr, Layout>::transposed(); }

	// BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	// auto operator~ (subarray const& self) { return self.transposed(); }
	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	auto operator~ (subarray& self) { return self.transposed(); }
	BOOST_MULTI_FRIEND_CONSTEXPR BOOST_MULTI_HD
	auto operator~ (subarray&& self) { return std::move(self).transposed(); }

	using const_subarray<T, D, ElementPtr, Layout>::reindexed;

	template<class... Indexes>
	constexpr auto reindexed(index first, Indexes... idxs) & -> subarray {
		return const_subarray<T, D, ElementPtr, Layout>::reindexed(first, idxs...);
		// return ((this->reindexed(first).rotated()).reindexed(idxs...)).unrotated();
	}
	template<class... Indexes>
	constexpr auto reindexed(index first, Indexes... idxs)&& -> subarray {
		return const_subarray<T, D, ElementPtr, Layout>::reindexed(first, idxs...);
		// return ((std::move(*this).reindexed(first).rotated()).reindexed(idxs...)).unrotated();
	}

	BOOST_MULTI_HD constexpr auto base() const& -> typename subarray::element_const_ptr { return this->base_; }
	BOOST_MULTI_HD constexpr auto base() &  -> ElementPtr { return this->base_; }
	BOOST_MULTI_HD constexpr auto base() && -> ElementPtr { return this->base_; }
	// BOOST_MULTI_HD constexpr auto base() const& -> element_const_ptr {return base_;}

	constexpr auto operator=(const_subarray<T, D, ElementPtr, Layout> const& other)      & -> subarray& {
		if(this == std::addressof(other)) {return *this;}
		assert(this->extension() == other.extension());
		this->elements() = other.elements();
		return *this;
	}

	constexpr void swap(subarray&& other) && noexcept {
		assert(this->extension() == other.extension());
		adl_swap_ranges(this->elements().begin(), this->elements().end(), std::move(other).elements().begin());
	}
	friend constexpr void swap(subarray&& self, subarray&& other) noexcept { std::move(self).swap(std::move(other)); }

	// template<class A, typename = std::enable_if_t<!std::is_base_of_v<subarray, std::decay_t<A>>>> friend constexpr void swap(subarray&& self, A&& other) noexcept { std::move(self).swap(std::forward<A>(other)); }
	// template<class A, typename = std::enable_if_t<!std::is_base_of_v<subarray, std::decay_t<A>>>> friend constexpr void swap(A&& other, subarray&& self) noexcept { std::move(self).swap(std::forward<A>(other)); }

	// template<class Array> constexpr void swap(Array&& other) && noexcept {
	//  assert( std::move(*this).extension() == other.extension() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  this->elements().swap(std::forward<Array>(other).elements());
	// //  adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	// }
	// template<class A> constexpr void swap(A&& other) & noexcept {return swap(std::forward<A>(other));}

	// friend constexpr void swap(subarray&& self, subarray&& other) noexcept {std::move(self).swap(std::move(other));}

	// template<class Array> constexpr void swap(subarray& self, Array&& other) noexcept {self.swap(std::forward<Array>(other));}  // TODO(correaa) remove
	// template<class Array> constexpr void swap(Array&& other, subarray& self) noexcept {self.swap(std::forward<Array>(other));}

	// fix mutation
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...> const& other) && -> decltype(auto) {operator=(          other ); return *this;}
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...> const& other)  & -> subarray& {
		assert(other.extensions() == this->extensions());
		this->elements() = other.elements();
		return *this;
	}

	// fix mutation
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...>     && other) && -> subarray& {operator=(std::move(other)); return *this;}
	template<class TT, class... As> constexpr auto operator=(const_subarray<TT, D, As...>     && other)  & -> subarray& {
		assert(this->extensions() == other.extensions());
		this->elements() = std::move(other).elements();
		return *this;
	}

	template<
		class Range,
		class = std::enable_if_t<! std::is_base_of_v<subarray, Range> >,
		class = std::enable_if_t<! is_subarray<Range>::value>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	constexpr auto operator=(Range const& rng) &  // TODO(correaa) check that you LHS is not read-only?
	-> subarray& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		assert(this->size() == static_cast<size_type>(adl_size(rng)));  // TODO(correaa) or use std::cmp_equal?
		adl_copy_n(adl_begin(rng), adl_size(rng), this->begin());
		return *this;
	}

	template<
		class Range,
		class = std::enable_if_t<! std::is_base_of_v<subarray, Range> >,  // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
		class = std::enable_if_t<! is_subarray<Range>::value          >  // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
	>
	constexpr auto operator=(Range const& rng) && -> subarray& {operator=(rng); return *this;}

	// template<class TT, class... As>
	// constexpr auto operator=(const_subarray<TT, D, As...> const& other) && -> subarray& {operator=(other); return *this;}

	// template<class TT, class... As>
	// constexpr auto operator=(subarray<TT, D, As...>&& other) && -> subarray& {operator=(std::move(other)); return *this;}

	// template<class TT, class... As>
	// constexpr
	// auto operator=(const_subarray<TT, D, As...> const& other) & -> subarray& {
	//  assert(this->extension() == other.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  this->elements() = other.elements();
	//  return *this;
	// }

	template<class TT, class... As>
	constexpr
	auto operator=(const_subarray<TT, D, As...> const& other) && -> subarray& {
		assert(this->extension() == other.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		this->elements() = other.elements();
		return *this;
	}

	template<class TT, class... As>
	constexpr
	auto operator=(subarray<TT, D, As...>&& other) & -> subarray& {
		assert(this->extension() == other.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		this->elements() = std::move(other).elements();
		return *this;
	}

	// template<class TT, class... As>
	// constexpr
	// auto operator=(array<TT, D, As...>&& other) & -> subarray& {
	//  operator=(static_cast<move_subarray<TT, D, As...>&&>(std::move(other)));
	//  other.clear();  // TODO(correaa) is this a good idea?
	//  return *this;
	// }

	constexpr auto operator=(const_subarray<T, D, ElementPtr, Layout> const& other) const&& -> subarray&;  // for std::indirectly_writable
	// {
	//  assert(this->extension() == other.extension());
	//  this->elements() = other.elements();
	//  return std::move(*this);
	// }

	constexpr auto operator=(subarray const& other) & -> subarray& {
		if(this == std::addressof(other)) { return *this; }
		assert(this->extension() == other.extension());
		this->elements() = other.elements();
		return *this;
	}
	constexpr auto operator=(subarray&& other) & noexcept -> subarray& {  // TODO(correaa) make conditionally noexcept
		if(this == std::addressof(other)) { return *this; }
		assert(this->extension() == other.extension());
		this->elements() = other.elements();
		return *this;
	}

	auto operator=(std::initializer_list<typename subarray::value_type> values) && -> subarray& {operator=(values); return *this; }
	auto operator=(std::initializer_list<typename subarray::value_type> values) &  -> subarray& {
		assert( static_cast<size_type>(values.size()) == this->size() );
		adl_copy_n(values.begin(), values.size(), this->begin());
		return *this;
	}

//  constexpr
//  auto operator=(const_subarray               const& other) & -> const_subarray& {
//      if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
//      assert(this->extension() == other.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
//      // MULTI_MARK_SCOPE("multi::operator= [D="+std::to_string(D)+"] from "+typeid(T).name()+" to "+typeid(T).name() );
//      elements() = other.elements();
// //      if(this->num_elements() == this->nelems() and o.num_elements() == this->nelems() and this->layout() == o.layout()) {
// //          adl_copy_n(o.base(), o.num_elements(), this->base());
// //      } else if(o.stride() < (~o).stride()) {
// //          adl_copy_n( (~o).begin(), (~o).size(), (~(*this)).begin() );
// //      } else {
// //          assign(o.begin());
// //      }
//      return *this;
//  }

	// constexpr auto operator=(const_subarray const& other) &&
	// -> const_subarray& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	//  if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
	//  operator=(other);
	//  return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	// }

	// template<
	//  class Range,
	//  class = decltype( static_cast<Range>(std::declval<const_subarray<T, D, ElementPtr, Layout> const&>()) )
	// >
	// constexpr explicit operator Range() const {
	//  return static_cast<const_subarray<T, D, ElementPtr, Layout> const&>(*this).operator Range();
	// }

	// BOOST_MULTI_HD constexpr auto operator()() const& -> const_subarray<T, D, ElementPtr, Layout> /*const*/ {return this->paren_aux_();}  // NOLINT(readability-redundant-access-specifiers,readability-const-return-type)

	// using reference = typename std::conditional_t<
	//  (D > 1),
	//  subarray<typename subarray::element, D-1, typename subarray::element_ptr>,
	//  typename std::iterator_traits<typename subarray::element_ptr>::reference
	// >;

	// using const_reference = typename std::conditional_t<
	//  (D > 1),
	//  const_subarray<typename subarray::element, D-1, element_const_ptr>,
	//  typename std::iterator_traits<element_const_ptr>::reference
	// >;

	// BOOST_MULTI_HD constexpr auto operator[](index idx) const&    { return static_cast<typename subarray::const_reference>(this->at_aux_(idx)); }  // TODO(correaa) use return type to cast
	using const_subarray<T, D, ElementPtr, Layout>::operator[];
	// BOOST_MULTI_HD constexpr auto operator[](index idx) const& { return const_subarray<T, D, ElementPtr, Layout>::operator[](idx); }
	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> typename subarray::reference { return this->at_aux_(idx) ; }
	BOOST_MULTI_HD constexpr auto operator[](index idx)  & -> typename subarray::reference { return this->at_aux_(idx) ; }

	using const_subarray<T, D, ElementPtr, Layout>::sliced;
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) && -> subarray { return const_subarray<T, D, ElementPtr, Layout>::sliced(first, last) ; }
	BOOST_MULTI_HD constexpr auto sliced(index first, index last)  & -> subarray { return const_subarray<T, D, ElementPtr, Layout>::sliced(first, last) ; }

	using const_subarray<T, D, ElementPtr, Layout>::range;
	BOOST_MULTI_HD constexpr auto range(index_range irng)     && -> decltype(auto) {return std::move(*this).sliced(irng.front(), irng.front() + irng.size());}
	BOOST_MULTI_HD constexpr auto range(index_range irng)      & -> decltype(auto) {return                  sliced(irng.front(), irng.front() + irng.size());}

 private:
	using const_subarray<T, D, ElementPtr, Layout>::paren_aux_;

	BOOST_MULTI_HD constexpr auto paren_aux_() &  {return subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_);}
	BOOST_MULTI_HD constexpr auto paren_aux_() && {return subarray<T, D, ElementPtr, Layout>(this->layout(), this->base_);}

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx)  & -> decltype(auto) { return operator[](idx); }
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx) && -> decltype(auto) { return operator[](idx); }

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args)  & -> decltype(auto) {return operator[](idx).paren_aux_(args...);}
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(index idx, As... args) && -> decltype(auto) {return operator[](idx).paren_aux_(args...);}

	template<class... As>
	BOOST_MULTI_HD constexpr auto paren_aux_(index_range irng, As... args)  & {
		return this->range(irng).rotated().paren_aux_(args...).unrotated();
	}
	template<class... As>
	BOOST_MULTI_HD constexpr auto paren_aux_(index_range irng, As... args) && {
		return std::move(*this).range(irng).rotated().paren_aux_(args...).unrotated();
	}

	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(intersecting_range<index> inr, As... args)      & -> decltype(auto) {return paren_aux_(intersection(this->extension(), inr), args...);}
	template<class... As> BOOST_MULTI_HD constexpr auto paren_aux_(intersecting_range<index> inr, As... args)     && -> decltype(auto) {return paren_aux_(intersection(this->extension(), inr), args...);}

 public:
	using const_subarray<T, D, ElementPtr, Layout>::operator();

	BOOST_MULTI_HD constexpr auto operator()()      & ->    subarray           {return this->paren_aux_();}
	BOOST_MULTI_HD constexpr auto operator()()     && ->    subarray           {return this->paren_aux_();}

	template<class A1 = irange>                                                                       BOOST_MULTI_HD constexpr auto operator()(A1 arg1)                                        &  -> decltype(auto) {return this->paren_aux_(arg1);}  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange>                                                    BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2)                               &  -> decltype(auto) {return this->paren_aux_(arg1, arg2);}
	template<class A1 = irange, class A2 = irange, class A3 = irange>                                 BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3)                      &  -> decltype(auto) {return this->paren_aux_(arg1, arg2, arg3);}  // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) &  -> decltype(auto) {return this->paren_aux_(arg1, arg2, arg3, arg4, args...);}  // NOLINT(whitespace/line_length) pattern line

	template<class A1 = irange>                                                                       BOOST_MULTI_HD constexpr auto operator()(A1 arg1)                                        && -> decltype(auto) {return std::move(*this).paren_aux_(arg1);}                             // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange>                                                    BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2)                               && -> decltype(auto) {return std::move(*this).paren_aux_(arg1, arg2);}                       // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange>                                 BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3)                      && -> decltype(auto) {return std::move(*this).paren_aux_(arg1, arg2, arg3);}                 // NOLINT(whitespace/line_length) pattern line
	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> BOOST_MULTI_HD constexpr auto operator()(A1 arg1, A2 arg2, A3 arg3, A4 arg4, As... args) && -> decltype(auto) {return std::move(*this).paren_aux_(arg1, arg2, arg3, arg4, args...);}  // NOLINT(whitespace/line_length) pattern line


 private:
	template<class Self, typename Tuple, std::size_t ... I>
	static BOOST_MULTI_HD constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...>/*012*/) -> decltype(auto) {
		using std::get;
		return std::forward<Self>(self)(get<I>(tuple)...);
	}

 public:
	using const_subarray<T, D, ElementPtr, Layout>::apply;
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl) && -> decltype(auto) { return apply_impl_(std::move(*this), tpl, std::make_index_sequence<std::tuple_size<Tuple>::value>()); }
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tpl)  & -> decltype(auto) { return apply_impl_(          *this , tpl, std::make_index_sequence<std::tuple_size<Tuple>::value>()); }

	using const_subarray<T, D, ElementPtr, Layout>::partitioned;
	BOOST_MULTI_HD constexpr auto partitioned(size_type size)      & -> subarray<T, D+1, typename subarray::element_ptr> { return this->partitioned_aux_(size); }
	BOOST_MULTI_HD constexpr auto partitioned(size_type size)     && -> subarray<T, D+1, typename subarray::element_ptr> { return this->partitioned_aux_(size); }

	using const_subarray<T, D, ElementPtr, Layout>::flatted;
	constexpr auto flatted() & {
		// assert(is_flattable() && "flatted doesn't work for all layouts!");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		multi::layout_t<D-1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return subarray<T, D-1, ElementPtr>(new_layout, this->base_);
	}
	constexpr auto flatted() && {return this->flatted();}

	using const_subarray<T, D, ElementPtr, Layout>::reinterpret_array_cast;
	// template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	// constexpr auto reinterpret_array_cast() & {
	//  assert(this->layout().stride() * static_cast<size_type>(sizeof(T)) % static_cast<size_type>(sizeof(T2)) == 0);

	//  return subarray<T2, 1, P2>{
	//    typename subarray::layout_type{this->layout().sub(), this->layout().stride() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2)), this->layout().offset() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2)), this->layout().nelems() * static_cast<size_type>(sizeof(T)) / static_cast<size_type>(sizeof(T2))},
	//    reinterpret_pointer_cast<P2>(this->base_)
    //     };
	// }

	// template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	// constexpr auto reinterpret_array_cast() && { return this->reinterpret_array_cast(); }

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast() & {
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return subarray<T2, D, P2>(
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)  // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		);
	}

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast() && {
		// static_assert( sizeof(T)%sizeof(T2) == 0,
		//  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return subarray<T2, D, P2>(
			this->layout().scale(sizeof(T), sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_pointer_cast<P2>(this->base_)  // if ADL gets confused here (e.g. multi:: and thrust::) then adl_reinterpret_pointer_cast will be necessary
		);
	}

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(size_type count) & {
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");

		assert( sizeof(T) == sizeof(T2)*static_cast<std::size_t>(count) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : checck implicit size compatibility
		if constexpr(std::is_pointer_v<ElementPtr>) {
			return subarray<T2, D + 1, P2>(
				layout_t<D+1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				static_cast<P2>(static_cast<void*>(this->base_))  // NOLINT(bugprone-casting-through-void) direct reinterepret_cast doesn't work here
			);
		} else {
			return subarray<T2, D + 1, P2>(
				layout_t<D+1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				reinterpret_cast<P2 const&>(this->base_)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
			);
		}
	}

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(size_type count) && {
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");

		assert( sizeof(T) == sizeof(T2)*static_cast<std::size_t>(count) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : checck implicit size compatibility
		if constexpr(std::is_pointer_v<ElementPtr>) {
			return subarray<T2, D + 1, P2>(
				layout_t<D+1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				static_cast<P2>(static_cast<void*>(this->base_))  // NOLINT(bugprone-casting-through-void) direct reinterepret_cast doesn't work here
			);
		} else {
			return subarray<T2, D + 1, P2>(
				layout_t<D+1>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, count).rotate(),
				reinterpret_cast<P2 const&>(this->base_)  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-casting-through-void) direct reinterepret_cast doesn't work here
			);
		}
	}

	using element_move_ptr  = multi::move_ptr<typename subarray::element, typename subarray::element_ptr>;

	constexpr auto element_moved()  & {return subarray<T, D, typename subarray::element_move_ptr, Layout>(this->layout(), element_move_ptr{this->base_});}
	constexpr auto element_moved() && {return element_moved();}
};

template<class Element, typename Ptr> struct array_iterator<Element, 0, Ptr>{};

template<class Element, typename Ptr, bool IsConst, bool IsMove>
struct array_iterator<Element, 1, Ptr, IsConst, IsMove>  // NOLINT(fuchsia-multiple-inheritance,cppcoreguidelines-pro-type-member-init,hicpp-member-init) stride_ is not initialized in some constructors
	: boost::multi::iterator_facade<
		array_iterator<Element, 1, Ptr, IsConst, IsMove>,
		Element, std::random_access_iterator_tag,
		std::conditional_t<
			IsConst,
			typename std::iterator_traits<typename std::pointer_traits<Ptr>::template rebind<Element const> >::reference,
			typename std::iterator_traits<Ptr>::reference
		>, multi::difference_type
	>
	, multi::affine          <array_iterator<Element, 1, Ptr, IsConst>, multi::difference_type>
	, multi::decrementable   <array_iterator<Element, 1, Ptr, IsConst>>
	, multi::incrementable   <array_iterator<Element, 1, Ptr, IsConst>>
	, multi::totally_ordered2<array_iterator<Element, 1, Ptr, IsConst>, void> {
	using affine = multi::affine<array_iterator<Element, 1, Ptr, IsConst>, multi::difference_type>;

	using pointer = std::conditional_t<
		IsConst,
		typename std::pointer_traits<Ptr>::template rebind<Element const>,
		Ptr
	>;

 private:
	using reference_aux = std::conditional_t<
		IsConst,
		typename std::iterator_traits<typename std::pointer_traits<Ptr>::template rebind<Element const> >::reference,
		typename std::iterator_traits<Ptr>::reference
	>;

 public:
	using reference = std::conditional_t<
		IsMove,
		std::add_rvalue_reference_t<reference_aux>,
		reference_aux
	>;

	using difference_type = typename affine::difference_type;
	static constexpr dimensionality_type dimensionality = 1;

	array_iterator() = default;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
	using layout_type = multi::layout_t<0>;

	template<
		bool OtherIsConst, std::enable_if_t< ! OtherIsConst, int> =0  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	> 
	BOOST_MULTI_HD constexpr explicit array_iterator(array_iterator<Element, 1, Ptr, OtherIsConst> const& other)
	: ptr_{other.base()}, stride_{other.stride()} {}

	template<
		class Other,
		decltype(multi::detail::implicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().base())* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	BOOST_MULTI_HD constexpr/*mplct*/ array_iterator(Other const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to reproduce the implicitness of the argument
	: ptr_{other.base()}, stride_{other.stride()} {}

	template<
		class Other,
		decltype(multi::detail::explicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().data_)* = nullptr
	>
	constexpr explicit array_iterator(Other const& other)
	: ptr_{other.data_}, stride_{other.stride_} {}

	template<class, dimensionality_type, class, bool, bool> friend struct array_iterator;

	template<
		class EElement, typename PPtr,
		typename = decltype(multi::detail::implicit_cast<Ptr>(std::declval<array_iterator<EElement, 1, PPtr>>().data_))
	>
	BOOST_MULTI_HD constexpr /*impl*/ array_iterator(array_iterator<EElement, 1, PPtr> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to reproduce the implicitness of original pointer
	: ptr_{other.base()}, stride_{other.stride_} {}

	constexpr explicit operator bool() const {return static_cast<bool>(this->ptr_);}

	BOOST_MULTI_HD constexpr auto operator[](typename array_iterator::difference_type n) const -> decltype(auto) {
		return *((*this) + n);
	}

	constexpr auto operator->() const {return static_cast<pointer>(ptr_);}

	using element = Element;
	using element_ptr = Ptr;
	// using pointer = element_ptr;
	using stride_type = multi::index;

	static constexpr dimensionality_type rank_v = 1;
	using rank = std::integral_constant<dimensionality_type, rank_v>;

	constexpr auto operator<(array_iterator const& other) const -> bool { 
		assert(other.stride_ == stride_);
		assert((ptr_ - other.ptr_)%stride_ == 0);
		return (ptr_ - other.ptr_)/stride_ < 0;
		// return distance_to_(other) > 0;
	}

	BOOST_MULTI_HD explicit constexpr array_iterator(Ptr ptr, typename const_subarray<Element, 1, Ptr>::index stride)
	: ptr_{ptr}, stride_{stride} {}

 private:
	friend struct const_subarray<Element, 1, Ptr>;

	element_ptr ptr_;  // {nullptr};  // TODO(correaa) : consider uninitialized pointer
	stride_type stride_;  // = {1};  // = {0};  // TODO(correaa) change to make it trivially default constructible

	// constexpr auto distance_to_(array_iterator const& other) const -> difference_type {
	//  assert(stride_==other.stride_ && (other.data_-data_)%stride_ == 0);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  return (other.data_ - data_)/stride_;  // with struct-overflow=3 error: assuming signed overflow does not occur when simplifying `X - Y > 0` to `X > Y` [-Werror=strict-overflow]
	// }

 public:
	BOOST_MULTI_HD constexpr auto operator+(difference_type n) const -> array_iterator { array_iterator ret{*this}; ret+=n; return ret; }
	BOOST_MULTI_HD constexpr auto operator-(difference_type n) const -> array_iterator { array_iterator ret{*this}; ret-=n; return ret; }

	BOOST_MULTI_HD constexpr auto base() const {return static_cast<pointer>(ptr_);}

	[[deprecated("use base() for iterator")]]
	BOOST_MULTI_HD constexpr auto data() const {return base();}

	BOOST_MULTI_FRIEND_CONSTEXPR
	auto base(array_iterator const& self) { return self.base(); }

	BOOST_MULTI_HD constexpr auto stride() const -> stride_type {return      stride_;}
	friend    constexpr auto stride(array_iterator const& self) -> stride_type { return self.stride_; }

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	constexpr auto operator++() -> array_iterator& {ptr_ += stride_; return *this;}
	constexpr auto operator--() -> array_iterator& {ptr_ -= stride_; return *this;}

	constexpr auto operator+=(difference_type n) -> array_iterator& {assert(stride_ != 0); ptr_ += stride_*n; return *this;}
	constexpr auto operator-=(difference_type n) -> array_iterator& {assert(stride_ != 0); ptr_ -= stride_*n; return *this;}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	constexpr auto operator-(array_iterator const& other) const -> difference_type {
		assert(stride_==other.stride_ && (ptr_ - other.ptr_)%stride_ == 0);
		return (ptr_ - other.ptr_)/stride_;  // with struct-overflow=3 error: assuming signed overflow does not occur when simplifying `X - Y > 0` to `X > Y` [-Werror=strict-overflow]
		// return -distance_to_(other);
	}

	constexpr auto operator==(array_iterator const& other) const -> bool {
		assert(this->stride_ == other.stride_);
		return this->ptr_ == other.ptr_;
	}

	template<bool OtherIsConst, std::enable_if_t<OtherIsConst != IsConst, int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator==(array_iterator<Element, 1, Ptr, OtherIsConst> const& other) const -> bool {
		assert(this->stride_ == other.stride_);
		return this->ptr_ == other.ptr_;
	}

	// constexpr auto operator!=(array_iterator<Element, 1, Ptr, true> const& other) const -> bool {
	//  assert(this->stride_ == other.stride_);
	//  return this->ptr_ != other.ptr_;
	// }

//  friend constexpr auto operator==(array_iterator const& self, array_iterator const& other) -> bool {return self.ptr_ == other.ptr_;}

	BOOST_MULTI_HD constexpr auto operator*() const -> decltype(auto) {
		if constexpr(IsMove) {
			return multi::move(*ptr_);
		} else {
			return static_cast<reference>(*ptr_);
		}
	}
};

template<class Element, dimensionality_type D, typename... Ts>
using iterator = array_iterator<Element, D, Ts...>;

template<typename T, typename ElementPtr, class Layout>
class const_subarray<T, 0, ElementPtr, Layout>
: public array_types<T, 0, ElementPtr, Layout> {
 public:
	using types = array_types<T, 0, ElementPtr, Layout>;
	using types::types;

	using element      = typename types::element;
	using element_ref  = typename std::iterator_traits<typename const_subarray::element_ptr>::reference;
	using element_cref = typename std::iterator_traits<typename const_subarray::element_const_ptr>::reference;
	using iterator = array_iterator<T, 0, ElementPtr>;

	using layout_type = Layout;

	constexpr auto operator= (element const& elem)  & -> const_subarray& {
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D=0 from "}+typeid(T).name()+" to "+typeid(T).name() );
		adl_copy_n(&elem, 1, this->base_);
		return *this;
	}
	constexpr auto operator= (element const& elem) && -> const_subarray& {
		operator=(elem);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	constexpr auto operator==(element const& elem) const -> bool {
		assert(this->num_elements() == 1);
		return adl_equal(&elem, std::next(&elem, this->num_elements()), this->base());
	}
	constexpr auto operator!=(element const& elem) const {return ! operator==(elem);}

	template<class Range0>
	constexpr
	auto operator=(Range0 const& rng) & -> const_subarray& {
		adl_copy_n(&rng, 1, this->base_);
		return *this;
	}

	constexpr auto elements_at(size_type idx [[maybe_unused]]) const& -> element_cref {assert(idx < this->num_elements()); return *(this->base_);}
	constexpr auto elements_at(size_type idx [[maybe_unused]])     && -> element_ref  {assert(idx < this->num_elements()); return *(this->base_);}
	constexpr auto elements_at(size_type idx [[maybe_unused]])      & -> element_ref  {assert(idx < this->num_elements()); return *(this->base_);}

	constexpr auto operator!=(const_subarray const& other) const {return ! adl_equal(other.base_, other.base_ + 1, this->base_);}
	constexpr auto operator==(const_subarray const& other) const {return   adl_equal(other.base_, other.base_ + 1, this->base_);}

	constexpr auto operator<(const_subarray const& other) const {
		return adl_lexicographical_compare(
			this->base_, this->base_ + this->num_elements(),
			other.base_, other.base_ + other.num_elements()
		);
	}

	using decay_type = typename types::element;

	BOOST_MULTI_HD constexpr auto operator()() const& -> element_ref {return *(this->base_);}  // NOLINT(hicpp-explicit-conversions)

	constexpr operator element_ref ()     && noexcept {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	constexpr operator element_ref ()      & noexcept {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	constexpr operator element_cref() const& noexcept {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax

	template<class IndexType>
	constexpr auto operator[](IndexType const&) const& = delete;
	constexpr auto sliced() const& = delete;
	constexpr auto partitioned() const& = delete;

	constexpr auto strided(difference_type) const& = delete;

	constexpr auto taked(difference_type) const& = delete;
	constexpr auto dropped(difference_type) const& = delete;

	BOOST_MULTI_HD constexpr auto reindexed() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto   rotated() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto unrotated() const& { return operator()(); }

	auto transposed() const& = delete;
	auto flatted() const& = delete;
	auto range() const& -> const_subarray = delete;

	using       cursor = cursor_t<typename const_subarray::element_ptr      , 0, typename const_subarray::strides_type>;
	using const_cursor = cursor_t<typename const_subarray::element_const_ptr, 0, typename const_subarray::strides_type>;
 private:
	constexpr auto home_aux_() const {return cursor(this->base_, this->strides());}

 public:
	constexpr auto home() const& -> const_cursor {return home_aux_();}

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	auto paren_aux_() const& {return operator()();}

 public:
	template<class Tuple>
	BOOST_MULTI_HD constexpr auto apply(Tuple const& /*unused*/) const {
		static_assert(std::tuple_size_v<Tuple> == 0);
		return operator()();
	}

	BOOST_MULTI_HD constexpr auto operator&() const& {  // NOLINT(google-runtime-operator)
		return /*TODO(correaa) add const*/ subarray_ptr<T, 0, ElementPtr, Layout>(this->base_, this->layout());
	}  // NOLINT(google-runtime-operator) extend semantics  //NOSONAR

	template<class T2, class P2 = typename std::pointer_traits<ElementPtr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast()  const& {
		return const_subarray<T2, 0, P2>{
			typename const_subarray::layout_type{this->layout()},
			reinterpret_pointer_cast<P2>(this->base_)
		};
	}

	constexpr auto broadcasted() const& {
		multi::layout_t<1> const new_layout{this->layout(), 0, 0, (std::numeric_limits<size_type>::max)()};  // paren for MSVC macros
		return subarray<T, 1, typename const_subarray::element_const_ptr>{new_layout, types::base_};
	}

	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int const /*version*/) {
		using AT = multi::archive_traits<Archive>;
		auto& element_ = *(this->base_);
		arxiv &     AT::make_nvp("element", element_);
	//  arxiv & cereal::make_nvp("element", element_);
	//  arxiv &                             element_ ;
	}
};

template<typename T, typename ElementPtr, class Layout>
struct const_subarray<T, 1, ElementPtr, Layout>  // NOLINT(fuchsia-multiple-inheritance) : to define operators via CRTP
	: multi::random_iterable<subarray<T, 1, ElementPtr, Layout> >  // paren for msvc 19.14?
	, array_types<T, ::boost::multi::dimensionality_type{1}, ElementPtr, Layout> {  // paren for msvc 19.14?
	~const_subarray() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	// boost serialization needs `delete`. void boost::serialization::extended_type_info_typeid<T>::destroy(const void*) const [with T = boost::multi::subarray<double, 1, double*, boost::multi::layout_t<1> >]
	// void operator delete(void* ptr) noexcept = delete;
	// void operator delete(void* ptr, void* place ) noexcept = doperator=(elete;  // NOLINT(bugprone-easily-swappable-parameters)

	static constexpr dimensionality_type rank_v = 1;
	using rank = std::integral_constant<dimensionality_type, rank_v>;

	using types = array_types<T, dimensionality_type{1}, ElementPtr, Layout>;
	using types::types;
	using layout_type = Layout;
	using ref_ = const_subarray;

	using element_type = T;

	using element_ptr       = typename types::element_ptr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element_type const>;
	using element_move_ptr  = multi::move_ptr<element_type, element_ptr>;
	using element_ref       = typename types::element_ref;
	using element_cref      = typename std::iterator_traits<element_const_ptr>::reference;

	using const_pointer     = element_const_ptr;
	using       pointer     = element_ptr;
	using const_reference   = typename array_types<T, dimensionality_type{1}, ElementPtr, Layout>::const_reference;
	using       reference   = typename array_types<T, dimensionality_type{1}, ElementPtr, Layout>::      reference;

	using default_allocator_type = typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type;

	constexpr auto get_allocator() const -> default_allocator_type {return default_allocator_of(const_subarray::base());}
	BOOST_MULTI_FRIEND_CONSTEXPR
	auto get_allocator(const_subarray const& self) -> default_allocator_type {return self.get_allocator();}

	using decay_type = array<typename types::element, dimensionality_type{1}, typename multi::pointer_traits<typename const_subarray::element_ptr>::default_allocator_type>;

	       constexpr auto decay()           const        -> decay_type {return decay_type{*this};}
	BOOST_MULTI_FRIEND_CONSTEXPR auto decay(const_subarray const& self) -> decay_type {return self.decay();}

	using basic_const_array = const_subarray<
		T, 1,
		typename std::pointer_traits<ElementPtr>::template rebind<typename const_subarray::element_type const>,
		Layout
	>;

 protected:
	template<class A> constexpr void intersection_assign(A&& other)&& {intersection_assign(std::forward<A>(other));}
	template<class A> constexpr void intersection_assign(A&& other)&  {  // NOLINT(cppcoreguidelines-rvalue-reference-param-not-moved,cppcoreguidelines-missing-std-forward) false positive clang-tidy 17
		std::for_each(
			intersection(types::extension(), extension(other)).begin(),
			intersection(types::extension(), extension(other)).end()  ,
			[&](auto const idx) {operator[](idx) = std::forward<A>(other)[idx];}
		);
	}

	const_subarray(const_subarray const&) = default;

	template<typename, ::boost::multi::dimensionality_type, typename EP, class LLayout> friend struct const_subarray;
	template<typename, ::boost::multi::dimensionality_type, class Alloc>                friend struct static_array;  // TODO(correaa) check if this is necessary

	template<class T2, class P2, class TT, dimensionality_type DD, class PP>
	friend constexpr auto static_array_cast(subarray<TT, DD, PP> const&) -> decltype(auto);

	// template<class T2>
	// friend constexpr auto reinterpret_array_cast(const_subarray&& self) {
	//  return std::move(self).template reinterpret_array_cast<T2, typename std::pointer_traits<element_ptr>::template rebind<T2>>();
	// }
	// template<class T2>
	// friend constexpr auto reinterpret_array_cast(const_subarray const& self) {
	//  return self.template reinterpret_array_cast<T2, typename std::pointer_traits<element_ptr>::template rebind<T2>>();
	// }

 public:
	friend constexpr auto sizes(const_subarray const& self) noexcept -> typename const_subarray::sizes_type {return self.sizes();}  // needed by nvcc
	friend constexpr auto size (const_subarray const& self) noexcept -> typename const_subarray::size_type  {return self.size ();}  // needed by nvcc

	constexpr auto operator+() const { return decay(); }

	const_subarray(const_subarray&&) noexcept = default;  // in C++ 14 this is necessary to return array references from functions
	// in c++17 things changed and non-moveable non-copyable types can be returned from functions and captured by auto

 protected:
	template<typename, multi::dimensionality_type, typename, class, bool> friend struct subarray_ptr;
	template<class, dimensionality_type D, class, bool, bool> friend struct array_iterator;

 public:
	friend constexpr auto dimensionality(const_subarray const& /*self*/) -> dimensionality_type {return 1;}

	// // NOLINTNEXTLINE(runtime/operator)
	// BOOST_MULTI_HD constexpr auto operator&()     && { return subarray_ptr<const_subarray, Layout>{this->base_, this->layout()}; }  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR
	// // NOLINTNEXTLINE(runtime/operator)
	// BOOST_MULTI_HD constexpr auto operator&()      & { return subarray_ptr<const_subarray, Layout>{this->base_, this->layout()}; } // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed  //NOSONAR
	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() const& {return const_subarray_ptr<T, 1, ElementPtr, Layout>{this->base_, this->layout()};}  // NOLINT(google-runtime-operator) extend semantics  //NOSONAR

	BOOST_MULTI_HD constexpr void assign(std::initializer_list<typename const_subarray::value_type> values) const {assert( values.size() == static_cast<std::size_t>(this->size()) );
		assign(values.begin(), values.end());
	}
	template<class It>
	constexpr auto assign(It first) & -> It {adl_copy_n(first, this->size(), this->begin()); std::advance(first, this->size()); return first;}
	template<class It>
	constexpr auto assign(It first)&& -> It {return assign(first);}
	template<class It>
	constexpr void assign(It first, It last) & {
		assert( std::distance(first, last) == this->size() ); (void)last;  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		assign(first);
	}
	template<class It>
	constexpr void assign(It first, It last)&& {assign(first, last);}

	// constexpr auto operator=(const_subarray&& other) & noexcept(std::is_nothrow_copy_assignable_v<T>) // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor,cppcoreguidelines-noexcept-move-operations) //NOSONAR
	// -> const_subarray& {  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	//  operator=(other);
	//  return *this;  // lints([cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	// }

	constexpr auto operator=(const_subarray const&) const& -> const_subarray const& = delete;
	constexpr auto operator=(const_subarray     &&) const& -> const_subarray const& = delete;

	// constexpr auto operator=(const_subarray const& other) && -> const_subarray& {
	//  if(this == std::addressof(other)) {return *this;}  // lints cert-oop54-cpp
	//  operator=(other); return *this;
	// }

	// [[deprecated("for compatibility with ranges")]] constexpr auto operator=(const_subarray const& other) const&& -> const_subarray const&&;  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) //NOSONAR this is needed to satify the std::indirectly_writable concept
	// {  // something like this will fail
	//  if(this == std::addressof(other)) {return static_cast<subarray const&&>(*this);}  // lints cert-oop54-cpp
	//  const_cast<subarray&&>(*this).operator=(other);
	//  return static_cast<subarray const&&>(*this);
	// }

	template<
		class ECPtr,
		class = std::enable_if_t< std::is_same_v<element_const_ptr, ECPtr> && !std::is_same_v<element_const_ptr, element_ptr> >  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	>
	constexpr auto operator=(const_subarray<T, 1L, ECPtr, Layout> const& other) const && -> const_subarray& {assert(0); operator=(          other ); return *this;}  // required by https://en.cppreference.com/w/cpp/iterator/indirectly_writable for std::ranges::copy_n

	using       cursor = cursor_t<typename const_subarray::element_ptr      , 1, typename const_subarray::strides_type>;
	using const_cursor = cursor_t<typename const_subarray::element_const_ptr, 1, typename const_subarray::strides_type>;
 private:
	constexpr auto home_aux_() const {return cursor(this->base_, this->strides());}

 public:
	constexpr auto home() const& -> const_cursor {return home_aux_();}

 private:
	template<typename, multi::dimensionality_type, typename, class> friend class subarray;

	BOOST_MULTI_HD constexpr auto at_aux_(index idx) const -> typename const_subarray::reference {  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment
	//  MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
		#endif
		auto ba = this->base_;  // NOLINT(llvm-qualified-auto,readability-qualified-auto)
		auto of = (idx*this->stride() - this->offset());  // NOLINT(llvm-qualified-auto,readability-qualified-auto)
		auto pt = ba + of;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		return *pt;  // in C++17 this is allowed even with syntethic references
	//  return *(this->base() + (idx*this->stride() - this->offset()));  // TODO(correaa) use this->base()[(i*this->stride() - this->offset())]
		#if defined(__clang__)
		#pragma clang diagnostic pop
		#endif
	}

 public:
	constexpr auto broadcasted() const& {
		multi::layout_t<2> const new_layout{this->layout(), 0, 0, (std::numeric_limits<size_type>::max)()};
		return const_subarray<T, 2, typename const_subarray::element_const_ptr>{new_layout, types::base_};
	}

	BOOST_MULTI_HD constexpr auto operator[](index idx) const& -> typename const_subarray::const_reference { return at_aux_(idx); }  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment
	//  BOOST_MULTI_HD constexpr auto operator[](index idx)      & -> typename const_subarray::      reference { return at_aux_(idx); }  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment
	//  BOOST_MULTI_HD constexpr auto operator[](index idx)     && -> typename const_subarray::      reference { return at_aux_(idx); }  // NOLINT(readability-const-return-type) fancy pointers can deref into const values to avoid assignment

	constexpr auto front() const& -> const_reference {return *begin();}
	constexpr auto back()  const& -> const_reference {return *std::prev(end(), 1);}

	constexpr auto front()     && ->       reference {return *begin();}
	constexpr auto back()      && ->       reference {return *std::prev(end(), 1);}

	constexpr auto front()      & ->       reference {return *begin();}
	constexpr auto back()       & ->       reference {return *std::prev(end(), 1);}

	// template<class ElementPtr2,
	//  std::enable_if_t<std::is_same_v<ElementPtr2, typename const_subarray::element_const_ptr>, int> = 0
	// >
	// constexpr operator subarray<T, 1, ElementPtr2, Layout>&& () const & {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) this is needed by std::ranges, TODO(correaa) think if this can be solved by inheritance from subarray<T, D, const ptr>
	//  return std::move(reinterpret_cast<subarray<T, 1, ElementPtr2, Layout> const&>(*this));  // NOLINT([ppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-reinterpret-cast)  think if this can be solved by inheritance from subarray<T, D, const ptr>
	// }

 private: 
	template<class Self, typename Tuple, std::size_t ... I, const_subarray* = nullptr>
	static constexpr auto apply_impl_(Self&& self, Tuple const& tuple, std::index_sequence<I...> /*012*/) -> decltype(auto) {
		return std::forward<Self>(self)(std::get<I>(tuple)...);
	}

 public:
	template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tuple) const& -> decltype(auto) {return apply_impl_(          *this , tuple, std::make_index_sequence<std::tuple_size_v<Tuple>>());}
	// template<typename Tuple> BOOST_MULTI_HD constexpr auto apply(Tuple const& tuple)     && -> decltype(auto) {return apply_impl_(std::move(*this), tuple, std::make_index_sequence<std::tuple_size_v<Tuple>>());}
	// template<typename Tuple>                constexpr auto apply(Tuple const& tuple)      & -> decltype(auto) {return apply_impl_(          *this , tuple, std::make_index_sequence<std::tuple_size_v<Tuple>>());}

	template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value == 0), int> = 0> BOOST_MULTI_HD constexpr auto operator[](Tuple const& /*empty*/) const& -> decltype(auto) { return *this; }  // NOLINT(modernize-use-constraints) TODO(correaa)
	template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value == 1), int> = 0> BOOST_MULTI_HD constexpr auto operator[](Tuple const& indices  ) const& -> decltype(auto) { return operator[](std::get<0>(indices)); }  // NOLINT(modernize-use-constraints) TODO(correaa)
	template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>::value >  1), int> = 0> BOOST_MULTI_HD constexpr auto operator[](Tuple const& indices  ) const&  // NOLINT(modernize-use-constraints) TODO(correaa)
	->decltype(operator[](std::get<0>(indices))[detail::tuple_tail(indices)]) {
		return operator[](std::get<0>(indices))[detail::tuple_tail(indices)]; }

    // Warning C4459 comes from boost::multi_array having a namespace indices which collides with the variable name?
    #ifdef _MSC_VER
    #pragma warning( push )
    #pragma warning( disable : 4459 )
    #endif

    [[deprecated("BMA compat, finish impl")]] BOOST_MULTI_HD constexpr auto operator[](std::tuple<irange> const& indices) const& { return (*this)({std::get<0>(indices).front(), std::get<0>(indices).back() + 1}); }

    #ifdef _MSC_VER
    #pragma warning( pop )
    #endif

	BOOST_MULTI_HD constexpr auto elements_at(size_type idx) const& -> decltype(auto) {assert(idx < this->num_elements()); return operator[](idx);}
	BOOST_MULTI_HD constexpr auto elements_at(size_type idx)     && -> decltype(auto) {assert(idx < this->num_elements()); return operator[](idx);}
	BOOST_MULTI_HD constexpr auto elements_at(size_type idx)      & -> decltype(auto) {assert(idx < this->num_elements()); return operator[](idx);}

	constexpr auto reindexed(index first) && {return reindexed(first);}
	constexpr auto reindexed(index first)  & {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return const_subarray{new_layout, types::base_};
	}

 private:
	BOOST_MULTI_HD constexpr auto taked_aux_(difference_type count) const {
		assert( count <= this->size() );  // calculating size is expensive that is why
		typename types::layout_t const new_layout{
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride()*count
		};
		return const_subarray{new_layout, this->base_};
	}

 public:
	constexpr auto taked(difference_type count) const& -> const_subarray<T, 1, ElementPtr, Layout> {return taked_aux_(count);}

 private:
	BOOST_MULTI_HD constexpr auto dropped_aux_(difference_type count) const -> const_subarray {
		assert( count <= this->size() );
		typename types::layout_t const new_layout{
			this->layout().sub(),
			this->layout().stride(),
			this->layout().offset(),
			this->stride()*(this->size() - count)
		};

		#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
		#endif

		return const_subarray{new_layout, this->base_ + (count*this->layout().stride() - this->layout().offset())};

		#if defined(__clang__)
		#pragma clang diagnostic pop
		#endif
	}

 public:
	constexpr auto dropped(difference_type count) const& -> const_subarray { return dropped_aux_(count); }

 private:

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	BOOST_MULTI_HD constexpr auto sliced_aux_(index first, index last) const {
		typename types::layout_t new_layout = this->layout();
		if(this->is_empty()) {
			assert(first == last);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
			new_layout.nelems() = 0;  // TODO(correaa) : don't use mutation
		} else {
			(new_layout.nelems() /= this->size())*=(last - first);
		}

		return const_subarray{new_layout, this->base_ + (first*this->layout().stride() - this->layout().offset())};
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	BOOST_MULTI_HD constexpr auto sliced(index first, index last) const& -> basic_const_array /*const*/ { return basic_const_array{sliced_aux_(first, last)};}  // NOLINT(readability-const-return-type)
	BOOST_MULTI_HD constexpr auto sliced(index first, index last)      & ->          const_subarray           { return                   sliced_aux_(first, last) ;}
	BOOST_MULTI_HD constexpr auto sliced(index first, index last)     && ->          const_subarray           { return                   sliced_aux_(first, last) ;}

	using  elements_iterator = elements_iterator_t<element_ptr      , layout_type>;
	using celements_iterator = elements_iterator_t<element_const_ptr, layout_type>;

	using       elements_range = elements_range_t<element_ptr      , layout_type>;
	using const_elements_range = elements_range_t<element_const_ptr, layout_type>;

 private:
	constexpr auto elements_aux_() const {return elements_range{this->base_, this->layout()};}

 public:
	constexpr auto  elements()      & ->       elements_range {return elements_aux_();}
	constexpr auto  elements()     && ->       elements_range {return elements_aux_();}
	constexpr auto  elements() const& -> const_elements_range {return const_elements_range{this->base(), this->layout()};}  // TODO(correaa) simplify

	constexpr auto celements() const  -> const_elements_range {return elements_aux_();}

	constexpr auto hull() const -> std::pair<element_const_ptr, size_type> {
		return {(std::min)(this->base(), this->base() + this->hull_size()), std::abs(this->hull_size())};  // paren for MSVC macros
	}

	/*[[gnu::pure]]*/ constexpr auto blocked(index first, index last)& -> const_subarray {
		return sliced(first, last).reindexed(first);
	}
	/*[[gnu::pure]]*/ constexpr auto stenciled(typename const_subarray::index_extension ext) -> const_subarray {
		return blocked(ext.first(), ext.last());
	}

 private:
	constexpr auto strided_aux_(difference_type diff) const -> const_subarray {
		auto const new_layout = typename types::layout_t{this->layout().sub(), this->layout().stride()*diff, this->layout().offset(), this->layout().nelems()};
		return {new_layout, types::base_};
	}

 public:
	constexpr auto strided(difference_type diff) const& -> basic_const_array { return strided_aux_(diff);}

	BOOST_MULTI_HD constexpr auto sliced(index first, index last, difference_type stride) const& -> basic_const_array { return sliced(first, last).strided(stride); }
	// BOOST_MULTI_HD constexpr auto sliced(index first, index last, difference_type stride)     && -> const_subarray { return sliced(first, last).strided(stride); }
	// BOOST_MULTI_HD constexpr auto sliced(index first, index last, difference_type stride)      & -> const_subarray { return sliced(first, last).strided(stride); }

	// BOOST_MULTI_HD constexpr auto range(index_range const& rng)      & {return                  sliced(rng.front(), rng.last());}
	// BOOST_MULTI_HD constexpr auto range(index_range const& rng)     && {return std::move(*this).sliced(rng.front(), rng.last());}
	BOOST_MULTI_HD constexpr auto range(index_range const& rng) const& {return                  sliced(rng.front(), rng.last());}

	BOOST_MULTI_HD constexpr auto operator()() const& -> const_subarray {return *this;}  // const_subarray(this->base(), this->layout());}
	// BOOST_MULTI_HD constexpr auto operator()()     && -> const_subarray       {return *this;}
	// BOOST_MULTI_HD constexpr auto operator()()      & -> const_subarray       {return *this;}

	// BOOST_MULTI_HD constexpr auto operator()(index_range const& rng)      & {return                  range(rng);}
	// BOOST_MULTI_HD constexpr auto operator()(index_range const& rng)     && {return std::move(*this).range(rng);}
	BOOST_MULTI_HD constexpr auto operator()(index_range const& rng) const& {return                  range(rng);}

	// BOOST_MULTI_HD constexpr auto operator()(index idx)      & -> decltype(auto) {return                  operator[](idx);}
	// BOOST_MULTI_HD constexpr auto operator()(index idx)     && -> decltype(auto) {return std::move(*this).operator[](idx);}
	BOOST_MULTI_HD constexpr auto operator()(index idx) const -> decltype(auto) {return operator[](idx);}

 private:
	// BOOST_MULTI_HD constexpr auto paren_aux_()      & {return operator()();}
	// BOOST_MULTI_HD constexpr auto paren_aux_()     && {return operator()();}
	BOOST_MULTI_HD constexpr auto paren_aux_() const& {return operator()();}

	// BOOST_MULTI_HD constexpr auto paren_aux_(index_range const& rng)      & {return range(rng);}
	// BOOST_MULTI_HD constexpr auto paren_aux_(index_range const& rng)     && {return range(rng);}
	BOOST_MULTI_HD constexpr auto paren_aux_(index_range const& rng) const& {return range(rng);}

	// BOOST_MULTI_HD constexpr auto paren_aux_(index idx)      & -> decltype(auto) {return operator[](idx);}
	// BOOST_MULTI_HD constexpr auto paren_aux_(index idx)     && -> decltype(auto) {return operator[](idx);}
	BOOST_MULTI_HD constexpr auto paren_aux_(index idx) const& -> decltype(auto) {return operator[](idx);}

	// constexpr auto paren_aux_(intersecting_range<index> const& rng)      & -> decltype(auto) {return                  paren_aux_(intersection(this->extension(), rng));}
	// constexpr auto paren_aux_(intersecting_range<index> const& rng)     && -> decltype(auto) {return std::move(*this).paren_aux_(intersection(this->extension(), rng));}
	constexpr auto paren_aux_(intersecting_range<index> const& rng) const& -> decltype(auto) {return                  paren_aux_(intersection(this->extension(), rng));}

 public:
	// constexpr auto operator()(intersecting_range<index> const& isrange)      & -> decltype(auto) {return                  paren_aux_(isrange);}
	// constexpr auto operator()(intersecting_range<index> const& isrange)     && -> decltype(auto) {return std::move(*this).paren_aux_(isrange);}
	BOOST_MULTI_HD constexpr auto operator()(intersecting_range<index> const& isrange) const& -> decltype(auto) {return                  paren_aux_(isrange);}

	// template<class... Args>
	// constexpr auto operator()(Args&&... args) &
	// ->decltype(paren_(*this, std::forward<Args>(args)...)) {
	//  return paren_(*this, std::forward<Args>(args)...); }

	// template<class... Args>
	// constexpr auto operator()(Args&&... args) &&
	// ->decltype(paren_(std::move(*this), std::forward<Args>(args)...)) {
	//  return paren_(std::move(*this), std::forward<Args>(args)...); }

	template<class... Args>
	BOOST_MULTI_HD constexpr auto operator()(Args&&... args) const&
	->decltype(paren_(*this, std::forward<Args>(args)...)) {
		return paren_(*this, std::forward<Args>(args)...); }

 private:
	BOOST_MULTI_HD constexpr auto partitioned_aux_(size_type size) const {
		assert( size != 0 );
		assert( (this->layout().nelems() % size) == 0 );  // TODO(correaa) remove assert? truncate left over? (like mathematica) // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		multi::layout_t<2> new_layout{this->layout(), this->layout().nelems()/size, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= size;  // TODO(correaa) : don't use mutation
		return subarray<T, 2, element_ptr>(new_layout, types::base_);
	}

 public:
	BOOST_MULTI_HD constexpr auto partitioned(size_type size) const& -> const_subarray<T, 2, element_ptr> {return partitioned_aux_(size);}

 private:
	BOOST_MULTI_HD constexpr auto chunked_aux_(size_type size) const {
		assert( this->size() % size == 0 );
		return partitioned_aux_(this->size()/size);
	}

 public:  // in Mathematica this is called Partition https://reference.wolfram.com/language/ref/Partition.html in RangesV3 it is called chunk
	BOOST_MULTI_HD constexpr auto chunked(size_type size) const& -> const_subarray<T, 2, element_ptr> {return chunked_aux_(size);}
	// BOOST_MULTI_HD constexpr auto chunked(size_type size)      & -> partitioned_type       {return chunked_aux_(size);}
	// BOOST_MULTI_HD constexpr auto chunked(size_type size)     && -> partitioned_type       {return chunked_aux_(size);}

	constexpr auto tiled(size_type count) const & {
		assert(count != 0);
		struct divided_type {
			const_subarray<T, 2, element_ptr> quotient;
			const_subarray<T, 1, element_ptr> remainder;
		};
		return divided_type{
			this->taked(this->size() - (this->size() % count)).chunked(count),
			this->dropped(this->size() - (this->size() % count))
		};
	}

 private:
	constexpr auto reversed_aux_() const -> const_subarray {
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}

 public:
	constexpr auto reversed() const& -> basic_const_array {return reversed_aux_();}
	constexpr auto reversed()      & -> const_subarray {return reversed_aux_();}
	constexpr auto reversed()     && -> const_subarray {return reversed_aux_();}

	friend constexpr auto reversed(const_subarray const& self) -> basic_const_array {return           self .reversed();}
	friend constexpr auto reversed(const_subarray      & self) -> const_subarray {return           self .reversed();}
	friend constexpr auto reversed(const_subarray     && self) -> const_subarray {return std::move(self).reversed();}

	// friend constexpr auto   rotated(const_subarray const& self) -> decltype(auto) {return self.  rotated();}
	// friend constexpr auto unrotated(const_subarray const& self) -> decltype(auto) {return self.unrotated();}

	// constexpr auto   rotated()      & -> decltype(auto) {return operator()();}
	// constexpr auto   rotated()     && -> decltype(auto) {return operator()();}
	BOOST_MULTI_HD constexpr auto   rotated() const& { return operator()(); }
	BOOST_MULTI_HD constexpr auto unrotated() const& { return operator()(); }

	auto transposed() const& = delete;
	auto flatted() const& = delete;

	using         iterator = typename multi::array_iterator<element_type, 1, typename types::element_ptr      >;
	using   const_iterator = typename multi::array_iterator<element_type, 1, typename types::element_ptr, true>;
	using    move_iterator = typename multi::array_iterator<element_type, 1, typename types::element_ptr, false, true>;

	using       reverse_iterator [[deprecated]] = std::reverse_iterator<      iterator>;
	using const_reverse_iterator [[deprecated]] = std::reverse_iterator<const_iterator>;

	struct [[deprecated("BMA compatibility")]] index_gen {auto operator[](irange const& rng) const {return std::make_tuple(rng);}};
	using extent_gen [[deprecated("BMA compatibility")]] = std::array<irange, 1>;
	using extent_range [[deprecated("BMA compatibility")]] = irange;

	template<
		class Range,
		std::enable_if_t<! has_extensions<std::decay_t<Range>>::value, int> =0,
		std::enable_if_t<! is_subarray<std::decay_t<Range>>::value, int> =0,
		class = decltype((void)std::declval<Range>().begin(), std::declval<Range>().end() ),
		class = decltype(Range{std::declval<typename const_subarray::const_iterator>(), std::declval<typename const_subarray::const_iterator>()})
	>
	constexpr explicit operator Range() const {
		// vvv Range{...} needed by Windows GCC?
		return Range{begin(), end()};  // NOLINT(fuchsia-default-arguments-calls) e.g. std::vector(it, it, alloc = {})
	}

 private:
	// [[deprecated("remove")]] BOOST_MULTI_HD constexpr explicit const_subarray(iterator begin, iterator end)
	// : const_subarray {
	//  layout_type{ {}/*begin->layout()*/, begin.stride(), 0, begin.stride()*(end - begin)},
	//  begin.base()
	// } {
	//  assert(begin.stride()  == end.stride() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	// //  assert(begin->layout() == end->layout());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	// }
	// friend constexpr auto ref<iterator>(iterator begin, iterator end) -> multi::subarray<typename iterator::element, iterator::rank_v, typename iterator::element_ptr>;

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	constexpr BOOST_MULTI_HD auto begin_aux_() const {return iterator{this->base_                  , this->stride()};}
	constexpr                auto end_aux_  () const {return iterator{this->base_ + types::nelems(), this->stride()};}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

 public:
	BOOST_MULTI_HD constexpr auto  begin() const& -> const_iterator {return begin_aux_();}
	               constexpr auto  begin()      & ->       iterator {return begin_aux_();}
	               constexpr auto  begin()     && ->       iterator {return begin_aux_();}

	// constexpr auto mbegin()      & {return move_iterator{begin()};}
	// constexpr auto mend  ()      & {return move_iterator{end  ()};}

	// constexpr auto mbegin()     && {return move_iterator{begin()};}
	// constexpr auto mend  ()     && {return move_iterator{end  ()};}

	constexpr auto  end  () const& -> const_iterator {return end_aux_();}
	constexpr auto  end  ()      & ->       iterator {return end_aux_();}
	constexpr auto  end  ()     && ->       iterator {return end_aux_();}

	[[deprecated("implement as negative stride")]] constexpr auto rbegin() const& {return const_reverse_iterator(end  ());}  // TODO(correaa) implement as negative stride?
	[[deprecated("implement as negative stride")]] constexpr auto rend  () const& {return const_reverse_iterator(begin());}  // TODO(correaa) implement as negative stride?

	BOOST_MULTI_FRIEND_CONSTEXPR auto begin(const_subarray const& self) -> const_iterator {return           self .begin();}
	BOOST_MULTI_FRIEND_CONSTEXPR auto begin(const_subarray      & self) ->       iterator {return           self .begin();}
	BOOST_MULTI_FRIEND_CONSTEXPR auto begin(const_subarray     && self) ->       iterator {return std::move(self).begin();}

	BOOST_MULTI_FRIEND_CONSTEXPR auto end  (const_subarray const& self) -> const_iterator {return           self .end()  ;}
	BOOST_MULTI_FRIEND_CONSTEXPR auto end  (const_subarray      & self) ->       iterator {return           self .end()  ;}
	BOOST_MULTI_FRIEND_CONSTEXPR auto end  (const_subarray     && self) ->       iterator {return std::move(self).end()  ;}

	BOOST_MULTI_HD constexpr auto cbegin()           const& -> const_iterator {return begin();}
	   constexpr auto cend  ()           const& -> const_iterator {return end()  ;}

	friend BOOST_MULTI_HD /*constexpr*/ auto cbegin(const_subarray const& self) {return self.cbegin();}
	BOOST_MULTI_FRIEND_CONSTEXPR auto cend  (const_subarray const& self) {return self.cend()  ;}

	// // fix mutation
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...> const& other) && -> const_subarray& {operator=(          other ); return *this;}
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...> const& other)  & -> const_subarray& {
	//  assert(other.extensions() == this->extensions());
	//  elements() = other.elements();
	//  return *this;
	// }

	// // fix mutation
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...>     && other) && -> const_subarray& {operator=(std::move(other)); return *this;}
	// template<class TT, class... As> constexpr auto operator=(const_subarray<TT, 1L, As...>     && other)  & -> const_subarray& {
	//  assert(this->extensions() == other.extensions());
	//  elements() = std::move(other).elements();
	//  return *this;
	// }

	// template<
	//  class Range,
	//  class = std::enable_if_t<! std::is_base_of_v<const_subarray, Range> >,
	//  class = std::enable_if_t<! is_subarray<Range>::value>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	// >
	// constexpr auto operator=(Range const& rng) &  // TODO(correaa) check that you LHS is not read-only?
	// -> const_subarray& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	//  assert(this->size() == static_cast<size_type>(adl_size(rng)));  // TODO(correaa) or use std::cmp_equal?
	//  adl_copy_n(adl_begin(rng), adl_size(rng), begin());
	//  return *this;
	// }
	// template<
	//  class Range,
	//  class = std::enable_if_t<! std::is_base_of_v<const_subarray, Range>>,
	//  class = std::enable_if_t<! is_subarray<Range>::value>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	// >
	// constexpr auto operator=(Range const& rng) && -> const_subarray& {operator=(rng); return *this;}

	template<class It> constexpr auto assign(It first) &&
	->decltype(adl_copy_n(first, std::declval<size_type>(), std::declval<iterator>()), void()) {
		return adl_copy_n(first, this->       size()      , std::move(*this).begin()), void(); }

	friend constexpr auto operator==(const_subarray const& self, const_subarray const& other) -> bool {
		return
			self.extension() == other.extension()
			&& self.elements() == other.elements()
		;
	}

	friend constexpr auto operator!=(const_subarray const& self, const_subarray const& other) -> bool {
		return
			self.extension() != other.extension()
			|| self.elements() != other.elements()
		;
	}

	template<class OtherT, typename OtherEP, class OtherLayout>
	friend constexpr auto operator==(const_subarray const& self, const_subarray<OtherT, 1, OtherEP, OtherLayout> const& other) -> bool {
		return
			self.extension() == other.extension()
			&& self.elements() == other.elements()
		;
	}

	template<class TT, typename EEPP, class LL>
	friend constexpr auto operator!=(const_subarray const& self, const_subarray<TT, 1, EEPP, LL> const& other) -> bool {
		return
			self.extension() != other.extension()
			|| self.elements() != other.elements()
		;
	}

	friend constexpr auto operator<(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(self, other); }
	friend constexpr auto operator>(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(other, self); }  // NOLINT(readability-suspicious-call-argument)

	friend constexpr auto operator<=(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(self, other) || self == other; }
	friend constexpr auto operator>=(const_subarray const& self, const_subarray const& other) -> bool { return lexicographical_compare_(other, self) || self == other; }  // NOLINT(readability-suspicious-call-argument)

 private:
	template<class A1, class A2>
	static constexpr auto lexicographical_compare_(A1 const& self, A2 const& other) -> bool {  // NOLINT(readability-suspicious-call-argument)
		if(self.extension().first() > other.extension().first()) {
			return true;
		}
		if(self.extension().first() < other.extension().first()) {
			return false;
		}
		return adl_lexicographical_compare(adl_begin(self), adl_end(self), adl_begin(other), adl_end(other));
	}

 public:
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() const -> subarray<T2, 1, P2> {  // name taken from std::static_pointer_cast
		return {this->layout(), static_cast<P2>(this->base_)};
	}
	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>, class... Args>
	constexpr auto static_array_cast(Args&&... args) const -> subarray<T2, 1, P2> {  // name taken from std::static_pointer_cast
		return {this->layout(), P2{this->base_, std::forward<Args>(args)...}};
	}

	template<class UF>
	constexpr auto element_transformed(UF&& fun) const& {
		return static_array_cast<
		//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
			std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
			transform_ptr<
			//  std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<UF const&, element_cref>>>,
				std::decay_t<std::invoke_result_t<UF const&, element_cref>>,
				UF, element_const_ptr, std::invoke_result_t<UF const&, element_cref>
			>
		>(std::forward<UF>(fun));
	}
	template<class UF>
	constexpr auto element_transformed(UF&& fun)  & {
		return static_array_cast<
			std::decay_t<std::invoke_result_t<UF const&, element_ref >>,
			transform_ptr<
				std::decay_t<std::invoke_result_t<UF const&, element_ref >>,
				UF, element_ptr      , std::invoke_result_t<UF const&, element_ref >
			>
		>(std::forward<UF>(fun));
	}
	template<class UF>
	constexpr auto element_transformed(UF&& fun) && {return element_transformed(std::forward<UF>(fun));}

	template<
		class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>,
		class Element = typename const_subarray::element,
		class PM = T2 std::decay_t<Element>::*
	>
	constexpr auto member_cast(PM member) const {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. "
			"Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements"
		);

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) reinterpret is what the function does. alternative for GCC/NVCC
		auto&& r1 = (*(reinterpret_cast<typename const_subarray::element_type* const&>(const_subarray::base_))).*member;  // ->*pm;
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa) find a better way
		auto* p1 = &r1; P2 p2 = reinterpret_cast<P2&>(p1);  //NOSONAR
#else
		auto p2 = static_cast<P2>(&(this->base_->*member));  // this crashes nvcc 11.2-11.4 and some? gcc compiler
#endif
		return subarray<T2, 1, P2>(this->layout().scale(sizeof(T), sizeof(T2)), p2);
	}

	// constexpr auto element_moved()  & {return subarray<typename subarray::element, D, element_move_ptr, Layout>(this->layout(), element_move_ptr{this->base_});}
	// constexpr auto element_moved() && {return element_moved();}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2>>
	constexpr auto reinterpret_array_cast()  const& {
		assert( this->layout().stride()*static_cast<size_type>(sizeof(T)) % static_cast<size_type>(sizeof(T2)) == 0 );

		return const_subarray<T2, 1, P2>{
			layout_type{this->layout().sub(), this->layout().stride()*static_cast<size_type>(sizeof(T))/static_cast<size_type>(sizeof(T2)), this->layout().offset()*static_cast<size_type>(sizeof(T))/static_cast<size_type>(sizeof(T2)), this->layout().nelems()*static_cast<size_type>(sizeof(T))/static_cast<size_type>(sizeof(T2))},
			reinterpret_pointer_cast<P2>(this->base_)
		};
	}

	template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2 const> >
	constexpr auto reinterpret_array_cast(size_type n) const& -> subarray<std::decay_t<T2>, 2, P2> {  // TODO(correaa) : use rebind for return type
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

		return subarray<std::decay_t<T2>, 2, P2>{
			layout_t<2>{this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, n},
			reinterpret_pointer_cast<P2>(this->base())
		}.rotated();
	}

	// // TODO(correaa) : rename to reinterpret_pointer_cast?
	// template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2> >
	// constexpr auto reinterpret_array_cast(size_type n)& {  // -> subarray<std::decay_t<T2>, 2, P2> {
	//  // static_assert( sizeof(T)%sizeof(T2)== 0,
	//  //  "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

	//  return subarray<std::decay_t<T2>, 2, P2>(
	//    layout_t<2>(this->layout().scale(sizeof(T), sizeof(T2)), 1, 0, n),
	//    reinterpret_pointer_cast<P2>(this->base())
	//  ).rotated();
	// }
	// template<class T2, class P2 = typename std::pointer_traits<element_ptr>::template rebind<T2> >
	// constexpr auto reinterpret_array_cast(size_type n)&& -> subarray<std::decay_t<T2>, 2, P2> {
	//  return this->reinterpret_array_cast<T2, P2>(n);
	// }

	template<class TT = typename const_subarray::element_type>
	constexpr auto fill(TT const& value) & -> decltype(auto) {
		return adl_fill_n(this->begin(), this->size(), value), *this;
	}
	constexpr auto fill()& -> decltype(auto) {return fill(typename const_subarray::element_type{});}

	template<class TT = typename const_subarray::element_type>
	[[deprecated]] constexpr auto fill(TT const& value) && -> decltype(auto) {return std::move(this->fill(value));}
	[[deprecated]] constexpr auto fill() && -> decltype(auto) {
		return std::move(*this).fill(typename const_subarray::element_type{});
	}

	template<class Archive>
	void serialize(Archive& arxiv, unsigned /*version*/) {
		using AT = multi::archive_traits<Archive>;
		std::for_each(this->begin(), this->end(), [&](reference& item) {arxiv & AT    ::make_nvp("item", item);});
	//  std::for_each(this->begin(), this->end(), [&](auto&&     item) {arxiv & cereal::make_nvp("item", item);});
	//  std::for_each(this->begin(), this->end(), [&](auto&&     item) {arxiv &                          item ;});
	}
};

template<class T2, class P2, class Array, class... Args>
constexpr auto static_array_cast(Array&& self, Args&&... args) -> decltype(auto) {
	return std::forward<Array>(self).template static_array_cast<T2, P2>(std::forward<Args>(args)...);
}

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref  // TODO(correaa) : inheredit from multi::partially_ordered2<array_ref<T, D, ElementPtr>, void>?
: subarray<T, D, ElementPtr>
{
	~array_ref() = default;  // lints(cppcoreguidelines-special-member-functions)

	using layout_type = typename array_ref::types::layout_t;

	using iterator = typename subarray<T, D, ElementPtr>::iterator;

 public:
	constexpr  // attempt for MSVC
	array_ref() = delete;  // because reference cannot be unbound

	array_ref(iterator, iterator) = delete;

	// return type removed for MSVC
	friend constexpr auto sizes(array_ref const& self) noexcept /*-> typename array_ref::sizes_type*/ {return self.sizes();}  // needed by nvcc
	friend constexpr auto size (array_ref const& self) noexcept /*-> typename array_ref::size_type*/  {return self.size ();}  // needed by nvcc

	[[deprecated("references are not copyable, use auto&&")]]
	array_ref(array_ref const&) = default;  // don't try to use `auto` for references, use `auto&&` or explicit value type

	#if defined(__NVCC__)
	array_ref(array_ref&&) noexcept = default;  // this needs to be public in nvcc c++17
	#else
	array_ref(array_ref&&) = delete;
	#endif

	#if defined(BOOST_MULTI_HAS_SPAN) && !defined(__NVCC__)
	template<class Tconst = const typename array_ref::element_type,
		std::enable_if_t<std::is_convertible_v<typename array_ref::element_const_ptr, Tconst*> && (D == 1), int> = 0  // NOLINT(modernize-use-constraints) TODO(correaa)
	>
	constexpr explicit operator std::span<Tconst>() const& {return std::span<Tconst>(this->data_elements(), this->size());}
	#endif

	template<class OtherPtr, class=std::enable_if_t<! std::is_same<OtherPtr, ElementPtr>{}>, decltype(multi::detail::explicit_cast<ElementPtr>(std::declval<OtherPtr>()))* = nullptr>
	constexpr explicit array_ref(array_ref<T, D, OtherPtr>&& other)
	: subarray<T, D, ElementPtr>{other.layout(), ElementPtr{std::move(other).base()}} {}  // cppcheck-suppress internalAstError ; bug in cppcheck 2.13.0

	template<class OtherPtr, class=std::enable_if_t<! std::is_same<OtherPtr, ElementPtr>{}>, decltype(multi::detail::implicit_cast<ElementPtr>(std::declval<OtherPtr>()))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	constexpr /*implicit*/ array_ref(array_ref<T, D, OtherPtr>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: subarray<T, D, ElementPtr>{other.layout(), ElementPtr{std::move(other).base()}} {}

	constexpr array_ref(ElementPtr dat, ::boost::multi::extensions_t<D> const& xs) /*noexcept*/  // TODO(correa) eliminate this ctor
	: subarray<T, D, ElementPtr>{typename subarray<T, D, ElementPtr>::types::layout_t(xs), dat} {}

	// constexpr array_ref(typename array_ref::extensions_type extensions, typename array_ref::element_ptr dat) noexcept
	// : subarray<T, D, ElementPtr>{typename array_ref::types::layout_t{extensions}, dat} {}

	constexpr array_ref(::boost::multi::extensions_t<D> exts, ElementPtr dat) noexcept
	: subarray<T, D, ElementPtr>{typename array_ref::types::layout_t(exts), dat} {}

	template<
		class Array,
		std::enable_if_t<! std::is_array_v<Array> > =0,
		std::enable_if_t<! std::is_base_of_v<array_ref, std::decay_t<Array>>, int> =0,
		std::enable_if_t<std::is_convertible_v<decltype(multi::data_elements(std::declval<Array&>())), ElementPtr>, int> =0  // NOLINT(modernize-use-constraints,ppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) support legacy c-arrays
	>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax and because a reference to c-array can be represented as an array_ref
	constexpr array_ref(  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax and because a reference to c-array can be represented as an array_ref
		Array& array  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	)
	: array_ref(
	 multi::data_elements(array),
	 extensions(array)
	) {}

	template<class TT = void, std::enable_if_t<sizeof(TT*) && D == 0, int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax and because a reference to c-array can be represented as an array_ref
	constexpr array_ref(  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax and because a reference to c-array can be represented as an array_ref
		T& elem // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	)
	: array_ref(&elem, {}) {}

	template<class TT, std::size_t N>
	constexpr array_ref(TT (&arr)[N])  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-explicit-constructor,hicpp-explicit-conversions) : for backward compatibility // NOSONAR
	: array_ref(
		::boost::multi::extensions(arr),
		::boost::multi::data_elements(arr)
	)
	{}

	template<
		class TT, std::size_t N//,
		//std::enable_if_t<std::is_convertible_v<decltype(data_elements(std::declval<TT(&)[N]>())), ElementPtr>, int> =0  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) support legacy c-arrays
	>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr array_ref(std::array<TT, N>& arr) : array_ref(::boost::multi::extensions(arr), ::boost::multi::data_elements(arr)) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) array_ptr is more general than pointer c-array support legacy c-arrays  // NOSONAR

//  this ctor makes memcheck complain about memmory used after scope
	template<class TT, std::enable_if_t<std::is_same_v<typename array_ref::value_type, TT>, int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor
	array_ref(std::initializer_list<TT> il) : array_ref(il.begin(), typename array_ref::extensions_type{static_cast<typename array_ref::size_type>(il.size())}) {}

	// template<class TT, std::enable_if_t<std::is_same_v<typename array_ref::value_type, TT>, int> =0>
	// array_ref(std::initializer_list<TT>&& il) = delete;

	using subarray<T, D, ElementPtr>::operator=;

 private:
	template<class It> constexpr auto copy_elements_(It first) {
		return adl_copy_n(first, this->num_elements(), this->data_elements());
	}

 public:
	BOOST_MULTI_HD constexpr auto data_elements() const& -> typename array_ref::element_const_ptr { return array_ref::base_; }

	template<class TT, class... As, std::enable_if_t<! std::is_base_of_v<array_ref, array_ref<TT, D, As...>> ,int> =0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr auto operator=(array_ref<TT, D, As...> const& other) && -> array_ref& {
		assert(this->extensions() == other.extensions());
		array_ref::copy_elements_(other.data_elements());
		return *this;
	}

	constexpr auto operator=(array_ref const& other) & -> array_ref& {
		if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
		// TODO(correaa) assert on extensions, not on num elements
		assert(this->num_elements() == other.num_elements());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		array_ref::copy_elements_(other.data_elements());
		return *this;
	}

	constexpr auto operator=(array_ref const& other) && -> array_ref& {
		if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
		operator=(other);
		return *this;
	}

	constexpr auto operator=(array_ref&& other) & noexcept(std::is_nothrow_copy_assignable_v<T>) // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor,cppcoreguidelines-noexcept-move-operations)  //NOSONAR(cppS5018)
	-> array_ref& {
		if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
		operator=(std::as_const(other));
		return *this;
	}
	constexpr auto operator=(array_ref&& other) && noexcept(std::is_nothrow_copy_assignable_v<T>) // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor,cppcoreguidelines-noexcept-move-operations)
	-> array_ref& {
		if(this == std::addressof(other)) {return *this;}  // lints(cert-oop54-cpp)
		operator=(std::as_const(other));
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
//  constexpr
	auto operator=(array_ref<TT, DD, As...> const& other)& -> array_ref& {
		assert( this->extensions() == other.extensions() );
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D="}+std::to_string(D)+" from "+typeid(TT).name()+" to "+typeid(T).name() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr auto operator=(array_ref<TT, DD, As...> const& other) && -> array_ref& {
		this->operator=(other);
		return *this;  // lints (cppcoreguidelines-c-copy-assignment-signature)
	}

	using  elements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_ptr      >;
	using celements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_const_ptr>;

 private:
	constexpr auto elements_aux_() const {
		return elements_type{
			this->base_,
			typename elements_type::extensions_type{multi::iextension{this->num_elements()}}
		};
	}

 public:
	       constexpr auto  elements()        const&       -> celements_type {return elements_aux_();}
	       constexpr auto  elements()             &       ->  elements_type {return elements_aux_();}
	       constexpr auto  elements()            &&       ->  elements_type {return elements_aux_();}

	friend constexpr auto elements(array_ref      & self) ->  elements_type {return           self . elements();}
	friend constexpr auto elements(array_ref     && self) ->  elements_type {return std::move(self). elements();}
	friend constexpr auto elements(array_ref const& self) -> celements_type {return           self . elements();}

	       constexpr auto celements()         const&       {return celements_type{array_ref::data_elements(), array_ref::num_elements()};}
	friend constexpr auto celements(array_ref const& self) {return self.celements();}

	template<typename TT, class... As>
	friend constexpr auto operator==(array_ref const& self, array_ref<TT, D, As...> const& other) -> bool {
		if(self.extensions() != other.extensions()) { return false; }

		#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
		#endif

		return adl_equal(
			other.data_elements(), other.data_elements() + other.num_elements(),
			self .data_elements()
		);

		#if defined(__clang__)
		#pragma clang diagnostic pop
		#endif

	}

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	template<typename TT, class... As>
	friend constexpr auto operator!=(array_ref const& self, array_ref<TT, D, As...> const& other) -> bool {
		if(self.extensions() != other.extensions()) { return true; }
		return !adl_equal(
			other.data_elements(), other.data_elements() + other.num_elements(),
			self .data_elements()
		);
		// return ! operator==(self, other);  // commented due to bug in nvcc 22.11
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif


	    BOOST_MULTI_HD constexpr auto data_elements() &      -> typename array_ref::element_ptr       { return array_ref::base_; }
	    BOOST_MULTI_HD constexpr auto data_elements() &&     -> typename array_ref::element_ptr       { return array_ref::base_; }
	//  BOOST_MULTI_HD constexpr auto data_elements() const& -> typename array_ref::element_const_ptr { return array_ref::base_; }

	friend constexpr auto data_elements(array_ref&& self) -> typename array_ref::element_ptr {return std::move(self).data_elements();}

	// data() is here for compatibility with std::vector
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data() const& { return data_elements(); }  // NOLINT(modernize-use-constraints) TODO(correaa)
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data()     && { return data_elements(); }  // NOLINT(modernize-use-constraints) TODO(correaa)
	template<class Dummy = void, std::enable_if_t<(D == 1) && sizeof(Dummy*), int> = 0> constexpr auto data()      & { return data_elements(); }  // NOLINT(modernize-use-constraints) TODO(correaa)

	// TODO(correaa) : find a way to use [[deprecated("use data_elements()")]] for friend functions
	friend constexpr auto data(array_ref const& self) -> typename array_ref::element_ptr { return           self .data_elements(); }
	friend constexpr auto data(array_ref      & self) -> typename array_ref::element_ptr { return           self .data_elements(); }
	friend constexpr auto data(array_ref     && self) -> typename array_ref::element_ptr { return std::move(self).data_elements(); }

	using decay_type = typename array_ref::decay_type;

	       constexpr auto decay()         const&       -> decay_type const& {return static_cast<decay_type const&>(*this);}
	friend constexpr auto decay(array_ref const& self) -> decay_type const& {return self.decay();}

 private:
	template<class TTN, std::size_t DD = 0>
	void check_sizes_() const {
		using std::get;
		if(size_type{get<DD>(this->sizes())} != size_type{std::extent<TTN, unsigned{DD}>::value}) {
			throw std::bad_cast{};
		}
		if constexpr(DD + 1 != D) {
			check_sizes_<TTN, DD + 1>();
		}
	}

	template<class TT> static auto launder_(TT* pointer) -> TT* {
	#if(defined(__cpp_lib_launder) && ( __cpp_lib_launder >= 201606L)) 
		return std::launder(pointer);
	#else
		return              pointer ;
	#endif
	}

	template<class, ::boost::multi::dimensionality_type, class> friend struct array;

	template<class TTN>
	constexpr auto to_carray_()& -> TTN& {
		check_sizes_<TTN>();
		return *launder_(reinterpret_cast<TTN*>(array_ref::base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	template<class TTN>
	constexpr auto to_carray_() const& -> TTN const& {
		check_sizes_<TTN>();
		return *launder_(reinterpret_cast<TTN const*>(array_ref::base_));  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

 public:
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr explicit operator TTN const&() const& { return to_carray_<TTN>(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20 
	constexpr explicit operator TTN&() && { return to_carray_<TTN>(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> = 0>  // NOLINT(modernize-use-constraints)  TODO(correaa) for C++20
	constexpr explicit operator TTN&() & { return to_carray_<TTN>(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)

 private:
	template<class Ar>
	auto serialize_structured_(Ar& arxiv, unsigned int const version) {
		subarray<T, D, ElementPtr>::serialize(arxiv, version);
	}
	template<class Archive>
	auto serialize_flat_(Archive& arxiv, unsigned int const /*version*/) {
		using AT = multi::archive_traits<Archive>;
		arxiv & AT::make_nvp("elements", AT::make_array(this->data_elements(), static_cast<std::size_t>(this->num_elements())));
	}
//  template<class Ar, class AT = multi::archive_traits<Ar>>
//  auto serialize_binary_if(std::true_type, Ar& ar) {
//      ar & AT::make_nvp("binary_data", AT::make_binary_object(this->data_elements(), static_cast<std::size_t>(this->num_elements())*sizeof(typename array_ref::element)));
//  }
//  template<class Ar>
//  auto serialize_binary_if(std::false_type, Ar& ar) {return serialize_flat(ar);}

 public:
	template<class Archive>
	auto serialize(Archive& arxiv, unsigned int const version) {
		serialize_flat(arxiv, version);
//      serialize_structured(ar, version);
//      switch(version) {
//          case static_cast<unsigned int>( 0): return serialize_flat(arxiv);
//          case static_cast<unsigned int>(-1): return serialize_structured(arxiv, version);
//      //  case 2: return serialize_binary_if(std::is_trivially_copy_assignable<typename array_ref::element>{}, arxiv);
//          default:
//              if( this->num_elements() <= version ){serialize_structured(arxiv, version);}
//              else                                 {serialize_flat      (arxiv         );}
//      }
	}
};

template<class T, dimensionality_type D, class Ptr = T*>
using array_cref = array_ref<
	std::decay_t<T>, D,
	typename std::pointer_traits<Ptr>::template rebind<T const>
>;

template<class T, dimensionality_type D, class Ptr = T*>
using array_mref = array_ref<
	std::decay_t<T>, D,
	std::move_iterator<Ptr>
>;

template<class TT, std::size_t N>
constexpr auto ref(
	TT(&arr)[N]  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) interact with legacy  // NOSONAR
) {
	return array_ref<std::remove_all_extents_t<TT[N]>, std::rank_v<TT[N]>>(arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) interact with legacy
}

template<class T, dimensionality_type D, typename Ptr = T*>
struct array_ptr
: subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_t> {
	using basic_ptr = subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_t>;

	constexpr array_ptr(Ptr data, multi::extensions_t<D> extensions)
	: basic_ptr{data, multi::layout_t<D>{extensions}} {}

	constexpr explicit array_ptr(std::nullptr_t nil) : array_ptr{nil, multi::extensions_t<D>{}} {}

	template<typename CArray>
	// cppcheck-suppress constParameterPointer ;  workaround cppcheck 2.11
	constexpr explicit array_ptr(CArray* data) : array_ptr{data_elements(*data), extensions(*data)} {}

	template<
		class TT, std::size_t N,
		std::enable_if_t<std::is_convertible_v<decltype(data_elements(std::declval<TT(&)[N]>())), Ptr>,int> =0  // NOLINT(modernize-use-constraints,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) support legacy c-arrays TODO(correaa) for C++20
	>
	// cppcheck-suppress noExplicitConstructor ;  // NOLINTNEXTLINE(runtime/explicit)
	constexpr array_ptr(TT(*array)[N]) : array_ptr{data_elements(*array), extensions(*array)} {}  // NOLINT(modernize-use-constraints,google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) array_ptr is more general than pointer c-array support legacy c-arrays  TODO(correaa) for C++20  // NOSONAR

	constexpr auto operator*() const -> array_ref<T, D, Ptr> {
		return array_ref<T, D, Ptr>((*static_cast<subarray_ptr<T, D, Ptr, typename array_ref<T, D, Ptr>::layout_t> const&>(*this)).extensions(), this->base());
	}
};

template<class T, typename Ptr>
class [[deprecated("no good uses found")]] array_ptr<T, 0, Ptr> {  // TODO(correaa) make it private mutable member
	mutable multi::array_ref<T, 0, Ptr> ref_;  // TODO(correaa) implement array_ptr like other cases

 public:
	~array_ptr()                              = default;
	constexpr array_ptr(array_ptr const&)     = default;
	constexpr array_ptr(array_ptr&&) noexcept = default;  // NOLINT(hicpp-noexcept-move,performance-noexcept-move-constructor) TODO(correaa) change the implementation like the other cases

	constexpr explicit array_ptr(Ptr dat, typename multi::array_ref<T, 0, Ptr>::extensions_type extensions) : ref_(dat, extensions) {}
	constexpr explicit array_ptr(Ptr dat) : array_ptr(dat, typename multi::array_ref<T, 0, Ptr>::extensions_type{}) {}

	constexpr explicit operator bool() const {return ref_.base();}
	constexpr explicit operator Ptr () const {return ref_.base();}

	auto operator=(array_ptr const&) -> array_ptr& = default;
	auto operator=(array_ptr     &&) noexcept -> array_ptr& = default;

	friend constexpr auto operator==(array_ptr const& self, array_ptr const& other) -> bool {return self.ref_.base() == other.ref_.base();}
	friend constexpr auto operator!=(array_ptr const& self, array_ptr const& other) -> bool {return self.ref_.base() != other.ref_.base();}

	constexpr auto operator* () const -> multi::array_ref<T, 0, Ptr>& {return  ref_;}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) make ref base class a mutable member
	constexpr auto operator->() const -> multi::array_ref<T, 0, Ptr>* {return &ref_;}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) make ref base class a mutable member
};

template<class TT, std::size_t N>
constexpr auto addressof(TT(&array)[N]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return array_ptr<
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		std::decay_t<std::remove_all_extents_t<TT[N]>>, static_cast<dimensionality_type>(std::rank<TT[N]>{}), std::remove_all_extents_t<TT[N]>*
	>{&array};
}

template<class T, dimensionality_type D, typename Ptr = T*>
using array_cptr = array_ptr<T, D,  typename std::pointer_traits<Ptr>::template rebind<T const>>;

template<dimensionality_type D, class P>
constexpr auto make_array_ref(P data, multi::extensions_t<D> extensions) {
	return array_ref<typename std::iterator_traits<P>::value_type, D, P>(data, extensions);
}

template<class P> auto make_array_ref(P data, extensions_t<0> exts) {return make_array_ref<0>(data, exts);}
template<class P> auto make_array_ref(P data, extensions_t<1> exts) {return make_array_ref<1>(data, exts);}
template<class P> auto make_array_ref(P data, extensions_t<2> exts) {return make_array_ref<2>(data, exts);}
template<class P> auto make_array_ref(P data, extensions_t<3> exts) {return make_array_ref<3>(data, exts);}
template<class P> auto make_array_ref(P data, extensions_t<4> exts) {return make_array_ref<4>(data, exts);}
template<class P> auto make_array_ref(P data, extensions_t<5> exts) {return make_array_ref<5>(data, exts);}

#if defined(__cpp_deduction_guides)

template<class It, typename V = typename std::iterator_traits<It>::value_type>  // pointer_traits doesn't have ::value_type
array_ptr(It)->array_ptr<V, 0, It>;

template<class It, typename V = typename std::iterator_traits<It>::value_type>  // pointer_traits doesn't have ::value_type
array_ptr(It, index_extensions<0>)->array_ptr<V, 0, It>;

template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<1>)->array_ptr<V, 1, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<2>)->array_ptr<V, 2, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<3>)->array_ptr<V, 3, It>;

template<
	class T, std::size_t N,
	typename V = std::remove_all_extents_t<T[N]>, std::size_t D = std::rank_v<T[N]>  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
>
array_ptr(T(*)[N])->array_ptr<V, static_cast<multi::dimensionality_type>(D)>;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility

template<class Ptr> array_ref(Ptr, index_extensions<0>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 0, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<1>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 1, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<2>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 2, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<3>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 3, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<4>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 4, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<5>) -> array_ref<typename std::iterator_traits<Ptr>::value_type, 5, Ptr>;

template<class It, class Tuple> array_ref(It, Tuple)->array_ref<typename std::iterator_traits<It>::value_type, std::tuple_size<Tuple>::value, It>;
#endif

// TODO(correaa) move to utility
template<class T, std::size_t N>
constexpr auto rotated(const T(&array)[N]) noexcept {                                                  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(array))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		base(array), extensions(array)
	).rotated();
}
template<class T, std::size_t N>
constexpr auto rotated(T(&array)[N]) noexcept {                                                        // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(array))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		base(array), extensions(array)
	).rotated();
}

template<class RandomAccessIterator, dimensionality_type D>
constexpr auto operator/(RandomAccessIterator data, multi::extensions_t<D> extensions)
-> multi::array_ptr<typename std::iterator_traits<RandomAccessIterator>::value_type, D, RandomAccessIterator>
{return {data, extensions};}

template<class In, class T, dimensionality_type N, class TP, class = std::enable_if_t<(N > 1)>, class = decltype((void)adl_begin(*In{}), adl_end(*In{}))>
constexpr auto uninitialized_copy
// require N>1 (this is important because it forces calling placement new on the pointer
(In first, In last, multi::array_iterator<T, N, TP> dest) {

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	while(first != last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		adl_uninitialized_copy(adl_begin(*first), adl_end(*first), adl_begin(*dest));
		++first;
		++dest;
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	return dest;
}

// begin and end for forwarding reference are needed in this namespace
// to overwrite the behavior of std::begin and std::end
// which take rvalue-references as const-references.

template<class T> auto begin(T&& rng) -> decltype(std::forward<T>(rng).begin()) {return std::forward<T>(rng).begin();}
template<class T> auto end  (T&& rng) -> decltype(std::forward<T>(rng).end()  ) {return std::forward<T>(rng).end()  ;}

template<class T, std::size_t N, std::size_t M>
auto transposed(T(&array)[N][M]) -> decltype(auto) {return ~multi::array_ref<T, 2>(array);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

template<class T, dimensionality_type D, class TPtr = T const*>
using array_const_view = array_ref<T, D, TPtr> const&;

template<class T, dimensionality_type D, class TPtr = T*>
using array_view = array_ref<T, D, TPtr>&;

}  // end namespace boost::multi

#ifndef BOOST_MULTI_SERIALIZATION_ARRAY_VERSION
#define BOOST_MULTI_SERIALIZATION_ARRAY_VERSION 0  // NOLINT(cppcoreguidelines-macro-usage) gives user opportunity to select serialization version //NOSONAR
// #define BOOST_MULTI_SERIALIZATION_ARRAY_VERSION  0 // save data as flat array
// #define BOOST_MULTI_SERIALIZATION_ARRAY_VERSION -1 // save data as structured nested labels array
// #define BOOST_MULTI_SERIALIZATION_ARRAY_VERSION 16 // any other value, structure for N <= 16, flat otherwise N > 16

namespace boost::multi {
	constexpr inline int serialization_array_version = BOOST_MULTI_SERIALIZATION_ARRAY_VERSION;
}  // end namespace boost::multi
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_ARRAY_REF_HPP_
