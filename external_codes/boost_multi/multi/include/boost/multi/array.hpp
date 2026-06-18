// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ARRAY_HPP
#define BOOST_MULTI_ARRAY_HPP
// #pragma once

#include "boost/multi/array_ref.hpp"  // IWYU pragma: export
#include "boost/multi/detail/adl.hpp"
#include "boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp"
#include "boost/multi/detail/memory.hpp"
#include "boost/multi/detail/static_allocator.hpp"  // IWYU pragma: export
#include "boost/multi/restriction.hpp"

#include <iterator>  // for std::sentinel_for
#include <memory>    // for std::allocator_traits
// #include <stdexcept>
#include <type_traits>  // for std::common_reference
#include <utility>      // for std::move

#if __has_include(<memory_resource>)
#include <memory_resource>
// Apple clang provides the header but not the compiled library prior to version 16
#if defined(__cpp_lib_memory_resource) && (__cpp_lib_memory_resource >= 201603) && !(defined(__APPLE__) && defined(__clang_major__) && __clang_major__ <= 15) && (!defined(_LIBCPP_VERSION) || !(_LIBCPP_VERSION <= 160001))
#define BOOST_MULTI_HAS_MEMORY_RESOURCE
#endif
#endif

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<concepts>) && __has_include(<ranges>)
#include <concepts>  // for constructible_from  // NOLINT(misc-include-cleaner)  // IWYU pragma: keep
#include <ranges>    // IWYU pragma: keep
#endif

// TODO(correaa) or should be (__CUDA__) or CUDA__ || HIP__
#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4626)  // assignment operator was implicitly defined as deleted
#endif

/// Branding namespace for the library
namespace boost {}

/// Main namespace for the library
namespace boost::multi {

namespace detail {

template<class Allocator>
struct array_allocator {
	using allocator_type = Allocator;
	array_allocator()    = default;

 private:
	BOOST_MULTI_NO_UNIQUE_ADDRESS allocator_type alloc_;

	using allocator_traits = multi::allocator_traits<allocator_type>;

	using size_type_ = typename allocator_traits::size_type;  // NOLINT(readability-redundant-typename) typename needed in C++17
	using pointer_   = typename allocator_traits::pointer;    // NOLINT(readability-redundant-typename) typename needed in C++17

 protected:
	constexpr auto alloc() & -> auto& { return alloc_; }
	constexpr auto alloc() const& -> allocator_type const& { return alloc_; }

	constexpr explicit array_allocator(allocator_type const& alloc) : alloc_{alloc} {}  // NOLINT(modernize-pass-by-value)

	constexpr auto allocate(size_type_ n) -> pointer_ {
		return n ? allocator_traits::allocate(alloc_, n) : pointer_{nullptr};
	}
	constexpr auto allocate(size_type_ n, typename allocator_traits::const_void_pointer hint) -> pointer_ {  // NOLINT(readability-redundant-typename) typename needed in C++17
		return n ? allocator_traits::allocate(alloc_, n, hint) : pointer_{nullptr};
	}

	constexpr auto uninitialized_fill_n(pointer_ first, size_type_ count, typename allocator_traits::value_type const& value) {  // NOLINT(readability-redundant-typename) typename needed in C++17
		return adl_alloc_uninitialized_fill_n(alloc_, first, count, value);
	}

	template<typename It, typename Size>
	constexpr auto uninitialized_copy_n(It first, Size n, pointer_ d_first) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename std::pointer_traits<pointer_>::element_type> && !multi::force_element_trivial_default_construction<typename std::pointer_traits<pointer_>::element_type>) {
			adl_alloc_uninitialized_default_construct_n(alloc_, d_first, n);
		}
		return adl_copy_n(first, n, d_first);
#else
		return adl_alloc_uninitialized_copy_n(alloc_, first, n, d_first);
#endif
	}

	template<typename It>
	auto uninitialized_move_n(It first, size_type_ count, pointer_ d_first) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename std::pointer_traits<pointer_>::element_type> && !multi::force_element_trivial_default_construction<typename std::pointer_traits<pointer_>::element_type>) {
			adl_alloc_uninitialized_default_construct_n(alloc_, d_first, count);
		}
		return adl_copy_n(std::make_move_iterator(first), count, d_first);
#else
		return adl_alloc_uninitialized_move_n(alloc_, first, count, d_first);
#endif
	}

	template<class ExecutionPolicy, typename It, typename Size>
	auto uninitialized_copy_n(ExecutionPolicy&& policy, It first, Size count, pointer_ d_first) {
		return adl_uninitialized_copy_n(std::forward<ExecutionPolicy>(policy), first, count, d_first);
	}

	template<typename It, typename Size>
	auto destroy_n(It first, Size n) { return adl_alloc_destroy_n(this->alloc(), first, n); }  // cppcheck-suppress functionStatic ; bug in cppcheck 2.19.0

 public:
#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	BOOST_MULTI_HD constexpr auto get_allocator() const noexcept -> allocator_type { return alloc_; }
};

}  // end namespace detail

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<class T, dimensionality_type D, class DummyAlloc = std::allocator<T>>  // DummyAlloc mechanism allows using the convention array<T, an_allocator<>>, is an_allocator supports void template argument
struct                                                                          // NOLINT(misc-multiple-inheritance) : used for composition
	dynamic_array
: protected detail::array_allocator<
	  typename allocator_traits<DummyAlloc>::template rebind_alloc<T>>
, public array_ref<T, D, typename multi::allocator_traits<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::pointer>
, private boost::multi::random_iterable<dynamic_array<T, D, typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>> {
	static_assert(
		std::is_same_v<
			std::remove_const_t<typename multi::allocator_traits<DummyAlloc>::value_type>,
			typename dynamic_array::element> ||
			std::is_same_v<
				std::remove_const_t<typename multi::allocator_traits<DummyAlloc>::value_type>,
				void>,  // allocator template can be redundant or void (which can be a default for the allocator)
		"allocator value type must match array value type"
	);

 protected:
	using array_alloc = detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>;

 public:
	using detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::get_allocator;

	/// Alloc
	using allocator_type = typename detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::allocator_type;  // NOLINT(readability-redundant-typename) needed for C++17
	/// Layout type (generally a strided layout `multi::layout_t<D>`)
	using layout_type    = typename array_ref<T, D, typename multi::allocator_traits<allocator_type>::pointer>::layout_type;  // NOLINT(readability-redundant-typename) needed for C++17
	/// Associalted array value type (generally itself, `multi::array<element, dimensionality, allocator_type>`)
	using decay_type     = array<T, D, allocator_type>;

	[[deprecated]] auto operator new(std::size_t count) -> void* { return ::operator new(count); }
	[[deprecated]] auto operator new(std::size_t count, void* ptr) -> void* { return ::operator new(count, ptr); }

	[[deprecated]] void operator delete(void* ptr) noexcept { ::operator delete(ptr); }  // this overrides the deleted delete operator in reference (base) class subarray

 protected:  // TODO(correaa) make private
	/// Associated array reference type, also its base class  (generally `multi::array_ref<element, dimensionality, allocator_type>`)
	using ref_ = array_ref<
		T, D,
		typename multi::allocator_traits<typename multi::allocator_traits<allocator_type>::template rebind_alloc<T>>::pointer>;

	using alloc_traits = typename multi::allocator_traits<allocator_type>;

	auto uninitialized_value_construct() {
		return adl_alloc_uninitialized_value_construct_n(dynamic_array::alloc(), this->base_, this->num_elements());
	}

	constexpr void uninitialized_default_construct() {
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), this->base_, this->num_elements());
		}
	}

	template<typename It> constexpr auto uninitialized_copy_elements(It first) {
		return array_alloc::uninitialized_copy_n(first, this->num_elements(), this->data_elements());
	}

	template<class ExecutionPolicy, typename It> auto uninitialized_copy_elements(ExecutionPolicy&& policy, It first) {
		return array_alloc::uninitialized_copy_n(std::forward<ExecutionPolicy>(policy), first, this->num_elements(), this->data_elements());
	}

	constexpr void destroy() {
		if constexpr(!(std::is_trivially_destructible_v<typename dynamic_array::element> || multi::force_element_trivial_destruction<typename dynamic_array::element>)) {
			array_alloc::destroy_n(this->data_elements(), this->num_elements());
		}
	}

	void allocate() {
		this->base_ = array_alloc::allocate(static_cast<typename multi::allocator_traits<typename dynamic_array::allocator_type>::size_type>(this->dynamic_array::num_elements()));  // NOLINT(readability-redundant-typename) needed for C++17
	}

 public:
	/// The value type associated after binding the first index (generally `multi::array<element, D - 1>`)
	using value_type = std::conditional_t<
		(D > 1),  // this parenthesis is needed
		array<typename dynamic_array::element, D - 1, allocator_type>,
		typename dynamic_array::element>;

	/// Signed integer type to represent difference between indices (usually `std::ptrdiff_t`)
	using typename ref_::difference_type;
	/// Integer type to represent sizes (usually `std::ptrdiff_t`)
	using typename ref_::size_type;

	explicit dynamic_array(allocator_type const& alloc) : array_alloc{alloc}, ref_(nullptr, {}) {}

	using ref_::operator();

	BOOST_MULTI_HD constexpr auto operator()() && -> decltype(auto) { return ref_::element_moved(); }

	using ref_::taked;

	constexpr auto taked(difference_type n) && -> decltype(auto) { return ref_::taked(n).element_moved(); }

	using ref_::dropped;

	constexpr auto dropped(difference_type n) && -> decltype(auto) { return ref_::dropped(n).element_moved(); }

	constexpr dynamic_array(dynamic_array&& other) /*noexcept(false)*/  // NOLINT(cppcoreguidelines-noexcept-move-operations,hicpp-noexcept-move,performance-noexcept-move-constructor,bugprone-exception-escape)
	: array_alloc{other.alloc()},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),  // NOLINT(readability-redundant-typename) needed for C++17
		  other.extents()
	  ) {
		adl_alloc_uninitialized_move_n(
			this->alloc(),
			other.data_elements(),
			other.num_elements(),
			this->data_elements()
		);
		(void)std::move(other);  //-V::530 this function moves the elements, not the whole object
	}

	constexpr dynamic_array(decay_type&& other, allocator_type const& alloc) noexcept
	: array_alloc{alloc}, ref_(std::exchange(other.base_, nullptr), other.extents()) {
		std::move(other).layout_mutable() = typename dynamic_array::layout_type(typename dynamic_array::extents_type{});  // = {};  careful! this is the place where layout can become invalid
	}

	explicit constexpr dynamic_array(decay_type&& other) noexcept
	: array_alloc{std::move(other.alloc())}, ref_(std::exchange(other.base_, nullptr), other.extents()) {
		std::move(other).layout_mutable() = typename dynamic_array::layout_type(typename dynamic_array::extents_type{});  // = {};  careful! this is the place where layout can become invalid
	}

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && (!defined(__clang_major__) || (__clang_major__ != 10))
	template<class It, std::sentinel_for<It> Sentinel = It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // NOLINT(readability-redundant-typename) needed for C++17
	constexpr explicit dynamic_array(It const& first, Sentinel const& last, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(layout_type{index_extension(adl_distance(first, last)) * multi::extents(*first)}.num_elements())),  // NOLINT(readability-redundant-typename) needed for C++17
		  index_extension(adl_distance(first, last)) * multi::extents(*first)
	  ) {
#if defined(__clang__) && defined(__CUDACC__)
		// TODO(correaa) add workaround for non-default constructible type and use adl_alloc_uninitialized_default_construct_n
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), ref_::data_elements(), ref_::num_elements());
		}
		adl_copy_n(first, last - first, ref_::begin());
#else
		adl_alloc_uninitialized_copy(dynamic_array::alloc(), first, last, ref_::begin());
#endif
	}
#else
	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>
	constexpr explicit dynamic_array(It const& first, It const& last, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(layout_type{index_extension(adl_distance(first, last)) * multi::extents(*first)}.num_elements())),
		  index_extension(adl_distance(first, last)) * multi::extents(*first)
	  ) {
#if defined(__clang__) && defined(__CUDACC__)
		// TODO(correaa) add workaround for non-default constructible type and use adl_alloc_uninitialized_default_construct_n
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), ref_::data_elements(), ref_::num_elements());
		}
		adl_copy_n(first, last - first, ref_::begin());
#else
		adl_alloc_uninitialized_copy(dynamic_array::alloc(), first, last, ref_::begin());
#endif
	}
#endif

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && (!defined(__clang_major__) || (__clang_major__ != 10))
	template<class It, std::sentinel_for<It> Sentinel, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // NOLINT(readability-redundant-typename) needed for C++17
	constexpr explicit dynamic_array(It const& first, Sentinel const& last)
	: dynamic_array(first, last, allocator_type{}) {}
#else
	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>
	constexpr explicit dynamic_array(It const& first, It const& last) : dynamic_array(first, last, allocator_type{}) {}
#endif

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)  //  && !defined(_MSC_VER)
 private:
	void extent_(typename dynamic_array::extents_type const& extensions) {  // NOLINT(readability-redundant-typename) needed for C++17
		auto new_layout = typename dynamic_array::layout_type{extensions};
		if(new_layout.num_elements() == 0) {
			return;
		}
		this->layout_mutable() = new_layout;  // typename array::layout_t{extensions};
		this->base_            = this->dynamic_array::array_alloc::allocate(
            static_cast<typename multi::allocator_traits<typename dynamic_array::allocator_type>::size_type>(  // NOLINT(readability-redundant-typename) needed for C++17
                new_layout.num_elements()
            ),
            this->data_elements()  // used as hint
        );
	}

 public:
	template<
		class Range, class = std::enable_if_t<!std::is_base_of_v<dynamic_array, std::decay_t<Range>>>,  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		class = decltype(std::declval<Range const&>().begin()),
		class = decltype(std::declval<Range const&>().end()),
		class = std::enable_if_t<!detail::is_subarray<Range const&>::value>>                                                // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
	requires std::is_convertible_v<std::ranges::range_reference_t<std::decay_t<std::ranges::range_reference_t<Range>>>, T>  //
		explicit dynamic_array(Range const& rng)                                                                            // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax  // NOSONAR
	: dynamic_array() {
		if(rng.size() == 0) {
			return;
		}
		auto outer_it = std::ranges::begin(rng);
		static_assert(D == 2);

		auto const& outer_ref   = *outer_it;
		auto        common_size = static_cast<size_type>(outer_ref.size());
		extent_({static_cast<size_type>(rng.size()), common_size});

		auto [is, js] = this->extents();
		{
			index const i        = 0;
			auto        inner_it = std::ranges::begin(outer_ref);
			for(auto j : js) {              // NOLINT(altera-unroll-loops) TODO(correa) change to algorithm applied on elements
				(*this)[i][j] = *inner_it;  // rng[i][j];
				++inner_it;
			}
			++outer_it;
		}

		for(index i = 1; i != is.last(); ++i) {
			auto const& outer_ref2 = *outer_it;
			assert(static_cast<multi::ssize_t>(outer_ref2.size()) == common_size);

			auto inner_it = std::ranges::begin(outer_ref2);
			for(auto j : js) {              // NOLINT(altera-unroll-loops) TODO(correa) change to algorithm applied on elements
				(*this)[i][j] = *inner_it;  // rng[i][j];
				++inner_it;
			}
			++outer_it;
		}
	}
#endif

	template<
		class Range, class = std::enable_if_t<!std::is_base_of_v<dynamic_array, std::decay_t<Range>>>,  // NOLINT(modernize-type-traits) bug in clang-tidy 19.1
		class = decltype(std::declval<Range const&>().begin()),
		class = decltype(std::declval<Range const&>().end()),
		class = std::enable_if_t<!detail::is_subarray<Range const&>::value>,  // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
		class = std::enable_if_t<(D == 1) || !std::is_convertible_v<decltype(*std::declval<Range const&>().begin()), T>>
#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)
		// for D==2 a genuine std::ranges range is owned by the constrained ranges-ctor above; step aside to avoid ambiguity
		// (input_range<Range> keeps the whole condition Range-dependent, so this is a soft SFINAE exclusion, not a hard error)
		,
		class = std::enable_if_t<!(D == 2 && std::ranges::input_range<Range>)>  // NOLINT(modernize-use-constraints) for C++20
#endif
		>
	// cppcheck-suppress noExplicitConstructor ; because I want to use assignment for lazy assigments form range-expressions // NOLINTNEXTLINE(runtime/explicit)
	/*explicit*/ dynamic_array(Range const& rng)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) : to allow terse syntax  // NOSONAR
	: dynamic_array(std::begin(rng), std::end(rng)) {
	}  // Sonar: Prefer free functions over member functions when handling objects of generic type "Range".

	template<class TT>
	auto uninitialized_fill_elements(TT const& value) {
		return array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), value);
	}

	template<class TT, class... As>
	dynamic_array(array_ref<TT, D, As...> const& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),  // NOLINT(readability-redundant-typename)
		  other.extents()
	  ) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
#else
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
#endif
	}

	explicit dynamic_array(typename dynamic_array::extents_type extensions, typename dynamic_array::element const& elem, allocator_type const& alloc)  // NOLINT(readability-redundant-typename)
	: array_alloc{alloc},
	  ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{extensions}.num_elements()), nullptr), extensions} {  // NOLINT(readability-redundant-typename)
		array_alloc::uninitialized_fill_n(this->data_elements(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);                             // NOLINT(readability-redundant-typename)
	}

	template<class Element>
	// for classic sfinae, needed by MSVC?
	explicit dynamic_array(Element const& elem, allocator_type const& alloc, std::enable_if_t<std::is_convertible_v<Element, typename dynamic_array::element> && (D == 0), int> /*dummy*/ = 0)  // if you get a compilation error here, you might be trying to initialize an array with a list of incorrect dimensionality
	: dynamic_array(typename dynamic_array::extents_type{}, elem, alloc) {}                                                                                                                  // NOLINT(readability-redundant-typename) for C++23

	template<
		class It,
		std::enable_if_t<std::is_convertible_v<typename std::iterator_traits<It>::value_type, T>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	// NOLINT(readability-redundant-typename)
	explicit constexpr dynamic_array(  // if you get a compilation error here, you might be trying to initialize an array with a list of incorrect dimensionality
		typename dynamic_array::extents_type exts, It elements_first
	)
	: array_alloc{},
	  array_ref<T, D, typename multi::allocator_traits<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::pointer>(  // NOLINT(readability-redundant-typename) for C++23
		  exts,
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type(exts).num_elements()),  // NOLINT(readability-redundant-typename) for C++23
								nullptr)
	  ) {
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), elements_first, this->num_elements(), this->elements().begin());
	}

	// NOLINT(readability-redundant-typename)
	explicit constexpr dynamic_array(  // if you get a compilation error here, you might be trying to initialize an array with a list of incorrect dimensionality
		typename dynamic_array::extents_type exts, typename dynamic_array::element const& elem
	)
	: array_alloc{},
	  array_ref<T, D, typename multi::allocator_traits<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::pointer>(  // NOLINT(readability-redundant-typename)
		  exts,
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type(exts).num_elements()),  // NOLINT(readability-redundant-typename)
								nullptr)
	  ) {
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element>) {
			array_alloc::uninitialized_fill_n(this->base(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);  // NOLINT(readability-redundant-typename) for C++23
		} else {                                                                                                                                                    // this workaround allows constexpr arrays for simple types
			adl_fill_n(this->base(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);                         // NOLINT(readability-redundant-typename) for C++23
		}
	}

	template<class ValueType, class = decltype(std::declval<ValueType>().extents()), std::enable_if_t<std::is_convertible_v<ValueType, typename dynamic_array::value_type>, int> = 0>                                                            // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(typename dynamic_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)                                                                                                        // fill constructor
	: array_alloc{alloc}, ref_(array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type(extension * value.extents()).num_elements())), extension * value.extents()) {  // NOLINT(readability-redundant-typename)
		static_assert(std::is_trivially_default_constructible_v<typename dynamic_array::element> || multi::force_element_trivial_default_construction<typename dynamic_array::element>);                                                         // TODO(correaa) not implemented for non-trivial types,
		adl_fill_n(this->begin(), this->size(), value);                                                                                                                                                                                          // TODO(correaa) implement via .elements()? substitute with uninitialized version of fill, uninitialized_fill_n?
	}

	template<class ValueType, class = decltype(std::declval<ValueType>().extents()), std::enable_if_t<std::is_convertible_v<ValueType, typename dynamic_array::value_type>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(typename dynamic_array::index_extension const& extension, ValueType const& value)
	: dynamic_array(extension, value, allocator_type{}) {}

	explicit dynamic_array(::boost::multi::extents_t<D> const& extensions, allocator_type const& alloc)
	: array_alloc{alloc}, ref_(array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{extensions}.num_elements())), extensions) {
		uninitialized_default_construct();
	}

	template<class... Args>
	static auto from_extensions(Args... exts) {
		return dynamic_array(::boost::multi::extents_t<D>(exts...));
	}

	explicit dynamic_array(::boost::multi::extents_t<D> const& exts)
	: dynamic_array(exts, allocator_type{}) {}

	// to make cling cppyy overload resolution easier
	template<class = void>  // gives low priority
	explicit dynamic_array(std::array<typename dynamic_array::size_type, static_cast<typename dynamic_array::dimensionality_type>(D)> const& exts)
	: dynamic_array(std::apply([](auto... sizes) -> auto { return typename dynamic_array::extents_type{sizes...}; }, exts)) {}

	template<class UninitilazedTag, std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_same_v<UninitilazedTag, ::boost::multi::uninitialized_elements_t>), int> = 0,                                                                  // NOLINT(modernize-use-constraints) for C++20
			 std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_trivially_default_constructible_v<typename dynamic_array::element> || multi::force_element_trivial_default_construction<typename dynamic_array::element>), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	explicit constexpr dynamic_array(::boost::multi::extents_t<D> const& extensions, UninitilazedTag /*unused*/, allocator_type const& alloc)
	: dynamic_array(extensions, alloc) {}

	template<class UninitilazedTag, std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_same_v<UninitilazedTag, ::boost::multi::uninitialized_elements_t>), int> = 0,                                                                    // NOLINT(modernize-use-constraints) for C++20
			 std::enable_if_t<sizeof(UninitilazedTag*) && (!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	[[deprecated("****element type cannot be partially formed (uninitialized), if you insists that this type should be treated as trivially constructible, consider opting-in to multi::force_trivial_default_construction at your own risk****")]]
	explicit constexpr dynamic_array(::boost::multi::extents_t<D> const& extensions, UninitilazedTag /*unusued*/) = delete /*[["****element type cannot be partially formed (uninitialized), if you insists that this type should be treated as trivially constructible, consider opting-in to multi::force_trivial_default_construction at your own risk****")]]*/;

	template<class UninitilazedTag, std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_same_v<UninitilazedTag, ::boost::multi::uninitialized_elements_t>), int> = 0,                                                                    // NOLINT(modernize-use-constraints) for C++20
			 std::enable_if_t<sizeof(UninitilazedTag*) && (!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	[[deprecated("****element type cannot be partially formed (uninitialized), if you insists that this type should be treated as trivially constructible, consider opting-in to multi::force_trivial_default_construction at your own risk****")]]
	explicit constexpr dynamic_array(::boost::multi::extents_t<D> const& extensions, UninitilazedTag /*unused*/, allocator_type const& /*alloc*/) = delete /*[["****element type cannot be partially formed (uninitialized), if you insists that this type should be treated as trivially constructible, consider opting-in to multi::force_trivial_default_construction at your own risk****"]]*/;

	template<class UninitilazedTag, std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_same_v<UninitilazedTag, ::boost::multi::uninitialized_elements_t>), int> = 0,                                                                  // NOLINT(modernize-use-constraints) for C++20
			 std::enable_if_t<sizeof(UninitilazedTag*) && (std::is_trivially_default_constructible_v<typename dynamic_array::element> || multi::force_element_trivial_default_construction<typename dynamic_array::element>), int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	explicit constexpr dynamic_array(::boost::multi::extents_t<D> const& extensions, UninitilazedTag /*unusued*/) : dynamic_array(extensions) {}

	template<class OtherT, class OtherEP, class OtherLayout, class = std::enable_if_t<std::is_assignable<typename ref_::element_ref, typename multi::subarray<OtherT, D, OtherEP, OtherLayout>::element>{}>, class = decltype(adl_copy(std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().begin(), std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().end(), std::declval<typename dynamic_array::iterator>()))>
	constexpr dynamic_array(multi::const_subarray<OtherT, D, OtherEP, OtherLayout> const& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{other.extents()}.num_elements())),
		  other.extents()
	  ) {
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), other.elements().begin(), this->num_elements(), this->data_elements());
	}

	template<class F>  // TODO(correaa) make more generic, e.g.: take ArrayWithElementsLike
	constexpr dynamic_array(multi::restriction<D, F> const& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{other.extents()}.num_elements())),
		  other.extents()
	  ) {
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), other.elements().begin(), this->num_elements(), this->data_elements());
	}

	template<class F>  // ArrayElementsLike, class = typename ArrayElementsLike::elements_t>
	// cppcheck-suppress noExplicitConstructor  // NOLINTNEXTLINE(runtime/explicit)
	constexpr dynamic_array(multi::restriction<D, F> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) to allow terse syntax
	: dynamic_array(other, allocator_type{}) {}

	template<class OtherT, class OtherEP, class OtherLayout, class = std::enable_if_t<std::is_assignable<typename ref_::element_ref, typename multi::subarray<OtherT, D, OtherEP, OtherLayout>::element>{}>, class = decltype(adl_copy(std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().begin(), std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().end(), std::declval<typename dynamic_array::iterator>()))>
	constexpr dynamic_array(multi::subarray<OtherT, D, OtherEP, OtherLayout>&& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{other.extents()}.num_elements())),
		  other.extents()
	  ) {
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), std::move(other).elements().begin(), this->num_elements(), this->data_elements());
	}

	template<
		class TT, class EElementPtr, class LLayout, std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout>&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename dynamic_array::iterator>()))>
	explicit dynamic_array(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other)
	: dynamic_array(other, allocator_type{}) {}

	template<
		class TT, class EElementPtr, class LLayout, std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename dynamic_array::iterator>()))>
	// cppcheck-suppress noExplicitConstructor  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ dynamic_array(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: dynamic_array(other, allocator_type{}) {}

	template<
		class TT, class EElementPtr, class LLayout, std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename dynamic_array::iterator>()))>
	// cppcheck-suppress noExplicitConstructor  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ dynamic_array(multi::subarray<TT, D, EElementPtr, LLayout>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: dynamic_array(std::move(other), allocator_type{}) {}

	// cppcheck-suppress noExplicitConstructor ; see below
	constexpr dynamic_array(multi::subarray<T, D, typename dynamic_array::element_ptr, typename dynamic_array::layout_type> const&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: dynamic_array(other, allocator_type{}) {}

	// cppcheck-suppress noExplicitConstructor ; see below
	constexpr dynamic_array(multi::const_subarray<T, D, typename dynamic_array::element_ptr, typename dynamic_array::layout_type> const&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: dynamic_array(other, allocator_type{}) {}

	// cppcheck-suppress noExplicitConstructor ; see below
	constexpr dynamic_array(multi::subarray<T, D, typename dynamic_array::element_ptr, typename dynamic_array::layout_type>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: dynamic_array(std::move(other), allocator_type{}) {}

	// NOLINTNEXTLINE(modernize-use-constraints) TODO(correaa) for C++20
	template<class TT, class... Args, std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&>().base()), T>, int> = 0>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ dynamic_array(array_ref<TT, D, Args...>& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: array_alloc{}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extents()} {
		dynamic_array::uninitialized_copy_elements(other.data_elements());
	}

	template<class TT, class... Args, std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&>().base()), T>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(array_ref<TT, D, Args...>& other)
	: array_alloc{}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extents()} {
		dynamic_array::uninitialized_copy_elements(other.data_elements());
	}

	// NOLINTNEXTLINE(modernize-use-constraints) TODO(correaa) for C++20
	template<class TT, class... Args, std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&&>().base()), T>, int> = 0>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ dynamic_array(array_ref<TT, D, Args...>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: array_alloc{}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extents()} {
		assert(this->stride() != 0);
		dynamic_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	template<class TT, class... Args, std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&&>().base()), T>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(array_ref<TT, D, Args...>&& other)
	: array_alloc{}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extents()} {
		assert(this->stride() != 0);
		dynamic_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	// NOLINTNEXTLINE(modernize-use-constraints) TODO(correaa) for C++20
	template<class TT, class... Args, std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...> const&>().base()), T>, int> = 0>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ dynamic_array(array_ref<TT, D, Args...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)  // NOSONAR
	: array_alloc{}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extents()} {
		assert(this->stride() != 0);
		dynamic_array::uninitialized_copy_elements(other.data_elements());
	}

	template<class TT, class... Args, std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...> const&>().base()), T>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(array_ref<TT, D, Args...> const& other)
	: array_alloc{},
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),
		  other.extents()
	  ) {
		assert(this->stride() != 0);
		dynamic_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	/// Copy constructor. Allocates new storage and copies all elements from @p other.
	/// The allocator is obtained respecting allocator propagation traits.
	/// @note O(n) where n is `other.num_elements()`.
	constexpr dynamic_array(dynamic_array const& other)
	: array_alloc(
		  multi::allocator_traits<allocator_type>::select_on_container_copy_construction(other.alloc())
	  ),
	  ref_(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),
		  other.extents()
	  ) {
		assert(this->stride() != 0);
		uninitialized_copy_elements(other.data_elements());
	}

	template<class ExecutionPolicy, std::enable_if_t<!std::is_convertible_v<ExecutionPolicy, typename dynamic_array::extents_type>, int> = 0>  // NOLINT(modernize-use-constraints,modernize-type-traits) TODO(correaa) for C++20
	explicit dynamic_array(ExecutionPolicy&& policy, dynamic_array const& other)
	: array_alloc{multi::allocator_traits<allocator_type>::select_on_container_copy_construction(other.alloc())}, ref_{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), extents(other)} {
		assert(this->stride() != 0);
		uninitialized_copy_elements(std::forward<ExecutionPolicy>(policy), other.data_elements());
	}

 private:
	using dynamic_value_type =
		std::conditional_t<
			(D != 1),
			dynamic_array<T, D - 1, allocator_type>,
			T>;

	template<typename, dimensionality_type, class> friend struct dynamic_array;

 public:
	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	constexpr dynamic_array(std::initializer_list<typename dynamic_array<T, D>::dynamic_value_type> values)
	: dynamic_array{(values.size() == 0) ? array<T, D>() : array<T, D>(values.begin(), values.end())} {}
	// ^^^ construct all with default constructor and copy to special memory at the end

 private:
	template<class Fun, class Tup>
	static auto std_apply_(Fun&& fun, Tup&& tup) {
		using std::apply;
		return apply(std::forward<Fun>(fun), std::forward<Tup>(tup));
	}
	struct make_from_tuple {
		template<typename... Elems>
		BOOST_MULTI_HD constexpr auto operator()(Elems... elems) const {
			return dynamic_array({static_cast<T const&>(elems)...});
		}
	};

 public:
	template<
		class Tuple,
		std::size_t                                                                                      = std::tuple_size<Tuple>::value,  // NOLINT(modernize-type-traits) TODO(correaa) remove or use tuple_size_v
		std::enable_if_t<                                                                                                                  // NOLINT(modernize-use-constraints)
			detail::all_elements_convertible_to<T, Tuple>::value && !multi::has_size<Tuple>::value, int> = 0>
	explicit constexpr dynamic_array(Tuple const& tup)
	: dynamic_array(std_apply_(make_from_tuple{}, tup)) {}

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	template<class TT = T, class = decltype(const_subarray<TT, D>(std::declval<std::initializer_list<std::initializer_list<TT>>>())), std::enable_if_t<multi::detail::is_implicitly_convertible_v<TT, T> && D == 2, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr dynamic_array(std::initializer_list<std::initializer_list<TT>> values)
	: dynamic_array{const_subarray<TT, D>(values)} {}  // construct all with default constructor and copy to special memory at the end

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	template<class TT = T, class = decltype(const_subarray<TT, D>(std::declval<std::initializer_list<std::initializer_list<std::initializer_list<TT>>>>())), std::enable_if_t<multi::detail::is_implicitly_convertible_v<TT, T> && D == 3, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	constexpr dynamic_array(std::initializer_list<std::initializer_list<std::initializer_list<TT>>> values)
	: dynamic_array{const_subarray<TT, D>(values)} {}  // construct all with default constructor and copy to special memory at the end

	explicit dynamic_array(
		std::initializer_list<typename dynamic_array<T, D>::value_type> values,
		allocator_type const&                                           alloc
	)
	: dynamic_array{(values.size() == 0) ? dynamic_array<T, D>() : dynamic_array<T, D>(values.begin(), values.end()), alloc} {}

	template<class TT, std::size_t N>
	constexpr explicit dynamic_array(TT (&array)[N])  // @SuppressWarnings(cpp:S5945) NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backward compatibility // NOSONAR
	: dynamic_array(std::begin(array), std::end(array)) {
		assert(this->stride() != 0);
	}

	/// returns an iterator to the beginning 
	constexpr auto begin() const& noexcept -> typename dynamic_array::const_iterator { return ref_::begin(); }

	/// returns an iterator to the end 
	constexpr auto end() const& noexcept -> typename dynamic_array::const_iterator { return ref_::end(); }

	constexpr auto begin() && noexcept -> typename dynamic_array::move_iterator { return ref_::begin(); }
	constexpr auto end() && noexcept -> typename dynamic_array::move_iterator { return ref_::end(); }

	constexpr auto begin() & noexcept -> typename dynamic_array::iterator { return ref_::begin(); }
	constexpr auto end() & noexcept -> typename dynamic_array::iterator { return ref_::end(); }

	/// Inherited indexed access
	using ref_::operator[];

	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> decltype(auto) {
		return multi::move(ref_::operator[](idx));
	}

	constexpr auto max_size() const noexcept { return static_cast<typename dynamic_array::size_type>(multi::allocator_traits<allocator_type>::max_size(this->alloc())); }  // TODO(correaa)  divide by nelements in under-dimensions?

 protected:
#ifdef __NVCC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress = 20011  // implicit __host__ __device__ ~dynamic_array [subobject] calls __host__ ~dynamic_array(); error attributed to deallocate() body
#endif
	constexpr void deallocate() {
		assert(this->stride() != 0);
		if(this->num_elements()) {
			multi::allocator_traits<allocator_type>::deallocate(this->alloc(), this->base_, static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()));
		}
	}

	void clear() noexcept {
		this->destroy();
		deallocate();
		this->layout_mutable() = typename dynamic_array::layout_type(typename dynamic_array::extents_type{});
		assert(this->stride() != 0);
	}

 public:
	constexpr dynamic_array() noexcept
	: array_alloc{}, ref_(nullptr, typename dynamic_array::extents_type{}) {
		assert(this->stride() != 0);
		assert(this->size() == 0);
	}

#if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
	constexpr ~dynamic_array()
#else
	~dynamic_array()
#endif
	{
		assert(this->stride() != 0);
		destroy();
		assert(this->stride() != 0);
		deallocate();
	}
#ifdef __NVCC__
#pragma nv_diagnostic pop
#endif

	/// The pointer type to a constant element (`std::allocator_traits<allocator_type>::const_pointer`, usually `T const*`)
	using element_const_ptr = typename std::allocator_traits<allocator_type>::const_pointer;

 private:
	/// A pointer wrapper that moves on dereferences
	using element_move_ptr = multi::move_ptr<typename dynamic_array::element_ptr>;

 public:
	/// Subarray reference after binding first index, `multi::subarray<T, D - 1, P>` or, for `D == 1`, `std::pointer_traits<pointer>::reference` (usually `T&`)
	using reference = std::conditional_t<
		(D > 1),
		subarray<typename dynamic_array::element, D - 1, typename dynamic_array::element_ptr>,
		std::conditional_t<
			D == 1,
			typename std::iterator_traits<typename dynamic_array::element_ptr>::reference,
			void>>;

	/// Subarray immutable reference after binding first index, `multi::const_subarray<T, D - 1, P>` or, for `D == 1`, `std::pointer_traits<const_pointer>::reference` (usually `T const&`)
	using const_reference = std::conditional_t<
		(D > 1),
		const_subarray<typename dynamic_array::element, D - 1, typename dynamic_array::element_ptr>,  // TODO(correaa) should be const_reference, but doesn't work witn rangev3?
		std::conditional_t<
			D == 1,
			decltype(*std::declval<typename dynamic_array::element_const_ptr>()),
			void>>;

	/// Random-access iterator in the leading dimension, in general they dereference to a subarray of lower dimension (`multi::subarray<...>`) or, for `D == 1` to an element reference (`T&`)
	using iterator       = multi::detail::array_iterator<T, D, typename dynamic_array::element_ptr>;
	/// Random-access iterator in the leading dimension, in general they dereference to an immutable subarrays of lower dimension (`multi::const_subarray<...>`) or, for `D == 1`, to an element immutable reference (`T const&`)
	using const_iterator = multi::detail::array_iterator<T, D, typename dynamic_array::element_ptr, true>;

	friend auto get_allocator(dynamic_array const& self) -> allocator_type { return self.get_allocator(); }

	// cppcheck-suppress duplInheritedMember ; to override
	BOOST_MULTI_HD constexpr auto data_elements() const& -> element_const_ptr { return this->base_; }

	// cppcheck-suppress duplInheritedMember ; to override
	BOOST_MULTI_HD constexpr auto data_elements() & -> typename dynamic_array::element_ptr { return this->base_; }

	// cppcheck-suppress duplInheritedMember ; to override
	BOOST_MULTI_HD constexpr auto data_elements() && -> typename dynamic_array::element_move_ptr { return std::make_move_iterator(this->base_); }

	// BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(dynamic_array const& self) { return self.data_elements(); }
	// BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(dynamic_array& self) { return self.data_elements(); }
	// BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(dynamic_array&& self) { return std::move(self).data_elements(); }

	/// Returns the base pointer of the array
	constexpr auto base() & -> typename dynamic_array::element_ptr { return ref_::base(); }
	constexpr auto base() const& -> typename dynamic_array::element_const_ptr { return typename dynamic_array::element_const_ptr{ref_::base()}; }

	constexpr auto origin() & -> typename dynamic_array::element_ptr { return ref_::origin(); }
	constexpr auto origin() const& -> typename dynamic_array::element_const_ptr { return ref_::origin(); }

	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(dynamic_array& self) -> typename dynamic_array::element_ptr { return self.origin(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(dynamic_array const& self) -> typename dynamic_array::element_const_ptr { return self.origin(); }

	template<class TT, typename EElementPtr, class LLayout>
	auto operator=(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other) -> dynamic_array& {
		ref_::operator=(other);  // TODO(correaa) : protect for self assigment
		assert(this->stride() != 0);
		return *this;
	}

	// (Restricted) copy-assignment, copies each element from the @p other array. Source and destination extents should match
	// @note Linear complexity in the number of elements
	auto operator=(dynamic_array const& other) & -> dynamic_array& {
		if(std::addressof(other) == this) {
			return *this;
		}  // cert-oop54-cpp
		assert(other.extents() == this->extents());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		assert(this->stride() != 0);
		return *this;
	}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

	// Element-move (deep move) assignment, moves each element from the @p other array. Source and destination extents should match
	// @note Linear complexity in the number of elements (cheaper than copy assignment if elements are effectively movable)
	constexpr auto operator=(dynamic_array&& other) noexcept -> dynamic_array& {                               // lints  (cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		assert(other.extents() == dynamic_array::extents());                                                   // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) there is no std::move_n algorithm
		assert(this->stride() != 0);
		return *this;
	}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

	/// Copy assignment from @p other array of a related typed
	template<class TT, class... As>  // , std::enable_if_t<std::is_assignable_v<T, TT>, int> = 0>  // NOLINT(modernize-use-constraints,modernize-type-traits) for C++20
	auto operator=(dynamic_array<TT, D, As...> const& other) & -> dynamic_array& {
		assert(extents(other) == dynamic_array::extents());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	/// Serializes elements into @p arxiv. Delegates to the base `array_ref::serialize`; shape is not saved.
	template<class Archive>
	void serialize(Archive& arxiv, unsigned int const version) { ref_::serialize(arxiv, version); }  // cppcheck-suppress duplInheritedMember ; to override

 private:
	void swap_(dynamic_array& other) noexcept { operator()().swap(other()); }  // cppcheck-suppress functionStatic

 public:
	friend void swap(dynamic_array& lhs, dynamic_array& rhs) noexcept { lhs.swap_(rhs); }

 private:
#ifdef __NVCOMPILER
#pragma diag_push
#pragma diag_suppress invalid_error_tag           // older nvc++ (e.g. 22) doesn't know the 'malloc_returns_non_pointer' tag below and would otherwise error on it
#pragma diag_suppress malloc_returns_non_pointer  // not recognized in nvc++ 22
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wignored-attributes"
#endif
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#endif
	template<class Ptr, std::enable_if_t<std::is_pointer_v<Ptr>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
#ifndef _BOOST_MULTI_SUPPRESS_GNU_MALLOC
#ifndef _MSC_VER
	[[gnu::malloc]] [[gnu::noinline]]
#endif
#endif
	BOOST_MULTI_HD static auto mallocate_me_(Ptr ptr) -> Ptr {
		return std::move(ptr);
	}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef __NVCOMPILER
#pragma diag_pop
#endif
	template<class Ptr, std::enable_if_t<!std::is_pointer_v<Ptr>, int> = 0>  // NOLINT(modernize-use-constraints) for C++20
	BOOST_MULTI_HD static auto mallocate_me_(Ptr ptr) -> Ptr { return std::move(ptr); }

 public:
	// BOOST_MULTI_HD constexpr auto splitted() & -> decltype(auto) { return ref_::splitted(); }

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
	BOOST_MULTI_HD constexpr auto splitted() && {
		multi::layout_t<1> const lyt1({}, this->layout().stride(), 0, this->layout().nelems() / this->layout().stride() / 2 * this->layout().stride());
		multi::layout_t<1> const lyt2({}, this->layout().stride(), 0, (this->layout().nelems() / this->layout().stride() + 1) / 2 * this->layout().stride());

		auto ptr1 = this->base_;                  // mallocate_me_(this->base_);                // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		auto ptr2 = this->base_ + lyt1.nelems();  // mallocate_me_(this->base_ + l1.nelems());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)

		return std::pair<
			subarray<T, 1, typename dynamic_array::element_ptr>,
			subarray<T, 1, typename dynamic_array::element_ptr>>(
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt1, ptr1),
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt2, ptr2)
		);
	}

	BOOST_MULTI_HD constexpr auto splitted() & {
		multi::layout_t<1> const lyt1({}, this->layout().stride(), 0, this->layout().nelems() / this->layout().stride() / 2 * this->layout().stride());
		multi::layout_t<1> const lyt2({}, this->layout().stride(), 0, (this->layout().nelems() / this->layout().stride() + 1) / 2 * this->layout().stride());

		auto ptr1 = this->base_;                  // mallocate_me_(this->base_);                // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		auto ptr2 = this->base_ + lyt1.nelems();  // mallocate_me_(this->base_ + l1.nelems());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)

		return std::pair<
			subarray<T, 1, typename dynamic_array::element_ptr>,
			subarray<T, 1, typename dynamic_array::element_ptr>>(
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt1, ptr1),
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt2, ptr2)
		);
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif
#ifdef __GNUC__
	[[gnu::always_inline]]
#endif
	BOOST_MULTI_HD constexpr auto fancy_splitted() && {
		multi::layout_t<1> const lyt1({}, this->layout().stride(), 0, this->layout().nelems() / this->layout().stride() / 2 * this->layout().stride());
		multi::layout_t<1> const lyt2({}, this->layout().stride(), 0, (this->layout().nelems() / this->layout().stride() + 1) / 2 * this->layout().stride());

		auto ptr1 = mallocate_me_(this->base_);                  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		auto ptr2 = mallocate_me_(this->base_ + lyt1.nelems());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)

#ifndef _BOOST_MULTI_SUPPRESS_ASSUMPTIONS
#if defined(__cpp_attributes_assume) && __cpp_attributes_assume >= 202207L
		[[assume(
			std::less_equal<>{}(ptr1 + lyt1.nelems(), ptr2) ||
			std::less_equal<>{}(ptr2 + lyt2.nelems(), ptr1)
		)]];
#endif
#endif
		return std::pair<
			subarray<T, 1, typename dynamic_array::element_ptr>,
			subarray<T, 1, typename dynamic_array::element_ptr>>(
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt1, ptr1),
			subarray<T, 1, typename dynamic_array::element_ptr>(lyt2, ptr2)
		);
	}

	BOOST_MULTI_HD constexpr auto fancy_splitted() & {
		multi::layout_t<1> const l1({}, this->layout().stride(), 0, this->layout().nelems() / this->layout().stride() / 2 * this->layout().stride());
		multi::layout_t<1> const l2({}, this->layout().stride(), 0, (this->layout().nelems() / this->layout().stride() + 1) / 2 * this->layout().stride());

		// auto p1 = this->base_;                // mallocate_me_(this->base_);                // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		// auto p2 = this->base_ + l1.nelems();  // mallocate_me_(this->base_ + l1.nelems());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)

		auto p1 = mallocate_me_(this->base_);                // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)
		auto p2 = mallocate_me_(this->base_ + l1.nelems());  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,llvm-qualified-auto,readability-qualified-auto)

#ifndef _BOOST_MULTI_SUPPRESS_ASSUMPTIONS
#if defined(__cpp_attributes_assume) && __cpp_attributes_assume >= 202207L
		[[assume(
			std::less_equal<>{}(p1 + l1.nelems(), p2) ||
			std::less_equal<>{}(p2 + l2.nelems(), p1)
		)]];
#endif
#endif
		return std::pair<
			subarray<T, 1, typename dynamic_array::element_ptr>,
			subarray<T, 1, typename dynamic_array::element_ptr>>(
			subarray<T, 1, typename dynamic_array::element_ptr>(l1, p1),
			subarray<T, 1, typename dynamic_array::element_ptr>(l2, p2)
		);
	}

#if defined(__clang__) && (__clang_major__ >= 16) && !defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#endif

	BOOST_MULTI_HD constexpr auto split() & {
		return std::move(*this).splitted();
	}
};

template<typename T, dimensionality_type D, class Alloc = std::allocator<T>>
using static_array [[deprecated("static_array has been renamed to dynamics_array (uses dynamic memory)")]] = dynamic_array<T, D, Alloc>;

#ifdef __clang__
#pragma clang diagnostic pop
#endif

/// A specialization for zero dimensions
///
/// The array might or might not contain an element
///
/// @tparam T Element type
/// @tparam Alloc Allocator type
template<typename T, class Alloc>
struct dynamic_array<T, 0, Alloc>  // NOLINT(misc-multiple-inheritance) : design
: protected detail::array_allocator<Alloc>
, public array_ref<T, 0, typename multi::allocator_traits<typename detail::array_allocator<Alloc>::allocator_type>::pointer> {
	static_assert(std::is_same_v<typename multi::allocator_traits<Alloc>::value_type, typename dynamic_array::element>, "allocator value type must match array value type");

 private:
	using array_alloc = detail::array_allocator<Alloc>;

 public:
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	/// (r-value address-of operator is deleted)
	/// @internal
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() && -> dynamic_array* = delete;  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : delete to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() & -> dynamic_array* { return this; }  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : override from base
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() const& -> dynamic_array const* { return this; }  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : override from base
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	using array_alloc::get_allocator;
	using allocator_type = typename dynamic_array::allocator_type;
	using decay_type     = array<T, 0, Alloc>;

	template<class Ptr>
	void assign(Ptr data) & {
		if(data) {
			assert(this->num_elements() == 1);
			adl_copy_n(data, this->num_elements(), this->base());
		}
	}

	template<
		class Singleton, std::enable_if_t<!std::is_base_of_v<dynamic_array, Singleton> && !std::is_same_v<Singleton, typename dynamic_array::element>, int> = 0,
		class = decltype(adl_copy_n(&std::declval<Singleton>(), 1, typename dynamic_array::element_ptr{}))>
	auto operator=(Singleton const& single) -> dynamic_array& {
		assign(&single);
		return *this;
	}

 protected:
	using alloc_traits = typename multi::allocator_traits<allocator_type>;
	using ref_         = array_ref<T, 0, typename multi::allocator_traits<typename multi::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;  // TODO(correaa) make private

	auto uninitialized_value_construct() {
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			return adl_alloc_uninitialized_value_construct_n(dynamic_array::alloc(), this->base_, this->num_elements());
		}
	}

	template<typename It> auto uninitialized_copy(It first) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(this->alloc(), this->data_elements(), this->num_elements());
		}
		return adl_copy(first, this->num_elements(), this->data_elements());
#else
		return adl_alloc_uninitialized_copy_n(this->alloc(), first, this->num_elements(), this->data_elements());
#endif
	}

	template<typename It>
	auto uninitialized_move(It first) {
		return adl_alloc_uninitialized_move_n(this->alloc(), first, this->num_elements(), this->data_elements());
	}

	constexpr void destroy() {
		if constexpr(!(std::is_trivially_destructible_v<typename dynamic_array::element> || multi::force_element_trivial_destruction<typename dynamic_array::element>)) {
			array_alloc::destroy_n(this->data_elements(), this->num_elements());
		}
	}

 public:
	using typename ref_::difference_type;
	using typename ref_::size_type;
	using typename ref_::value_type;

	constexpr explicit dynamic_array(allocator_type const& alloc) : array_alloc{alloc} {}

	BOOST_MULTI_HD constexpr dynamic_array(decay_type&& other, allocator_type const& alloc) noexcept  // 6b
	: array_alloc{alloc}, ref_{other.base_, other.extents()} {
		std::move(other).ref_::layout_t::operator=({});
	}

	using ref_::operator==;
	using ref_::operator!=;

	explicit dynamic_array(
		typename dynamic_array::extents_type const& extensions,
		typename dynamic_array::element const& elem, allocator_type const& alloc
	)
	: array_alloc{alloc},
	  ref_(dynamic_array::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{extensions}.num_elements())), extensions) {
		uninitialized_fill(elem);
	}

	explicit dynamic_array(typename dynamic_array::element const& elem, allocator_type const& alloc)
	: dynamic_array(typename dynamic_array::extents_type{}, elem, alloc) {}

	template<typename OtherT, typename OtherEPtr, class OtherLayout>
	explicit dynamic_array(multi::const_subarray<OtherT, 0, OtherEPtr, OtherLayout> const& other, allocator_type const& alloc)
	: array_alloc{alloc}, ref_(dynamic_array::allocate(other.num_elements()), extents(other)) {
		assert(other.num_elements() <= 1);
		if(other.num_elements()) {
#if defined(__clang__) && defined(__CUDACC__)
			if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
				adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), this->data_elements(), this->num_elements());
			}
			adl_copy(other.base(), other.base() + other.num_elements(), this->base());
#else
			adl_alloc_uninitialized_copy(dynamic_array::alloc(), other.base(), other.base() + other.num_elements(), this->base());
#endif
		}
	}

	template<class TT, class... Args>
	explicit dynamic_array(multi::dynamic_array<TT, 0, Args...> const& other, allocator_type const& alloc)  // TODO(correaa) : call other constructor (above)
	: array_alloc{alloc}, ref_(dynamic_array::allocate(static_cast<typename std::allocator_traits<Alloc>::size_type>(other.num_elements())), extents(other)) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
#else
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
#endif
	}

	template<class TT, class... Args>
	explicit dynamic_array(multi::dynamic_array<TT, 0, Args...> const& other)
	: dynamic_array(other, allocator_type{}) {}

	auto uninitialized_fill(typename dynamic_array::element const& elem) {
		array_alloc::uninitialized_fill_n(
			this->base_,
			static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()),
			elem
		);
	}

	template<class TT, class... Args>
	auto operator=(multi::const_subarray<TT, 0, Args...> const& other) -> dynamic_array& {
		adl_copy_n(other.base(), 1, this->base());
		return *this;
	}

	explicit dynamic_array(
		typename dynamic_array::extents_type const& exts,
		typename dynamic_array::element const&      elem
	)
	: array_alloc{}, ref_(dynamic_array::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename dynamic_array::layout_type{exts}.num_elements()), nullptr), exts) {
		uninitialized_fill(elem);
	}

	dynamic_array() : dynamic_array(multi::iextensions<0>{}) {}  // TODO(correaa) a noexcept will force a partially formed state for zero dimensional arrays

	explicit dynamic_array(typename dynamic_array::element const& elem)
	: dynamic_array(multi::iextensions<0>{}, elem) {}

	template<
		class Singleton, std::enable_if_t<!std::is_base_of_v<dynamic_array, Singleton> && !std::is_same_v<Singleton, typename dynamic_array::element>, int> = 0,  // NOLINT(modernize-type-traits) for C++20
		class = decltype(adl_copy_n(&std::declval<Singleton>(), 1, typename dynamic_array::element_ptr{}))>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax  // NOLINTNEXTLINE(runtime/explicit)
	/*implict*/ dynamic_array(Singleton const& single)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) this is used by the
	: ref_(dynamic_array::allocate(1), typename dynamic_array::extents_type{}) {
#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(!std::is_trivially_default_constructible_v<typename dynamic_array::element> && !multi::force_element_trivial_default_construction<typename dynamic_array::element>) {
			adl_alloc_uninitialized_default_construct_n(dynamic_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n(&single, typename multi::allocator_traits<Alloc>::size_type{1}, this->data_elements());
#else
		adl_alloc_uninitialized_copy_n(dynamic_array::alloc(), &single, typename multi::allocator_traits<Alloc>::size_type{1}, this->data_elements());
#endif
	}

	template<class ValueType, typename = std::enable_if_t<std::is_same_v<ValueType, value_type>>>                                          // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(typename dynamic_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)  // 3
	: dynamic_array(extension * extents(value), alloc) {
		assert(this->stride() != 0);
		using std::fill;
		fill(this->begin(), this->end(), value);
	}

	template<class ValueType, typename = std::enable_if_t<std::is_same_v<ValueType, value_type>>>             // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit dynamic_array(typename dynamic_array::index_extension const& extension, ValueType const& value)  // 3
	: dynamic_array(extension, value, allocator_type{}) {}

	explicit dynamic_array(typename dynamic_array::extents_type const& extensions, allocator_type const& alloc)  // 3
	: array_alloc{alloc}, ref_(dynamic_array::allocate(typename dynamic_array::layout_t{extensions}.num_elements()), extensions) {
		// assert(this->stride() != 0);
		uninitialized_value_construct();
	}

	explicit dynamic_array(typename dynamic_array::extents_type const& extensions)  // 3
	: dynamic_array(extensions, allocator_type{}) {}

	dynamic_array(dynamic_array const& other, allocator_type const& alloc)  // 5b
	: array_alloc{alloc}, ref_(dynamic_array::allocate(other.num_elements()), extents(other)) {
		assert(this->stride() != 0);
		uninitialized_copy_(other.data_elements());
	}

	dynamic_array(dynamic_array const& other)  // 5b
	: array_alloc{other.get_allocator()}, ref_{dynamic_array::allocate(other.num_elements(), other.data_elements()), {}} {
		assert(this->stride() != 0);
		uninitialized_copy(other.data_elements());
	}

	dynamic_array(dynamic_array&& other) noexcept
	: array_alloc{other.get_allocator()}, ref_(std::exchange(other.base_, nullptr), other.extents()) {  // should this move the elements? or move the object? or should be deleted?
		other.layout_mutable() = {};                                                                    // TODO(correaa) eliminate use of mutable member
		// other.layout_t<0>::operator=({});
		// , ref_(dynamic_array::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), other.extents()) {
		//  adl_alloc_uninitialized_move_n(
		//      this->alloc(),
		//      other.data_elements(),
		//      other.num_elements(),
		//      this->data_elements()
		//  );
		(void)std::move(other);  //-V::530 this function moves the elements, not the whole object
	}

 protected:
	void deallocate() {  // TODO(correaa) : move this to detail::array_allocator
		if(this->num_elements() && this->base_) {
			multi::allocator_traits<allocator_type>::deallocate(this->alloc(), this->base_, static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()));
		}
	}
	void clear() noexcept {
		this->destroy();
		deallocate();
		layout_t<0>::operator=({});
	}

 public:
	~dynamic_array() noexcept {
		this->destroy();
		deallocate();
	}
	using element_const_ptr = typename std::pointer_traits<typename dynamic_array::element_ptr>::template rebind<typename dynamic_array::element const>;

	BOOST_MULTI_FRIEND_CONSTEXPR auto get_allocator(dynamic_array const& self) -> allocator_type { return self.get_allocator(); }

	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	constexpr auto base() & -> typename dynamic_array::element_ptr { return ref_::base(); }
	constexpr auto base() const& -> typename dynamic_array::element_const_ptr { return ref_::base(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	BOOST_MULTI_FRIEND_CONSTEXPR auto base(dynamic_array& self) -> typename dynamic_array::element_ptr { return self.base(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto base(dynamic_array const& self) -> typename dynamic_array::element_const_ptr { return self.base(); }

	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	constexpr auto origin() & -> typename dynamic_array::element_ptr { return ref_::origin(); }
	constexpr auto origin() const& -> typename dynamic_array::element_const_ptr { return ref_::origin(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(dynamic_array& self) -> typename dynamic_array::element_ptr { return self.origin(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(dynamic_array const& self) -> typename dynamic_array::element_const_ptr { return self.origin(); }

	// NOSONAR
	constexpr operator typename std::iterator_traits<typename dynamic_array::element_const_ptr>::reference() const& {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
		return *(this->base_);
	}

	// NOSONAR
	constexpr operator std::add_rvalue_reference_t<typename std::iterator_traits<typename dynamic_array::element_ptr>::reference>() && {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
		return std::move(*(this->base_));
	}

	// NOSONAR
	constexpr operator typename std::iterator_traits<typename dynamic_array::element_ptr>::reference() & {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
		return *(this->base_);
	}

	template<class OtherElement, std::enable_if_t<std::is_convertible_v<typename dynamic_array::element, OtherElement>, int> = 0>  // NOLINT(modernize-use-constraints)
	// cppcheck-suppress duplInheritedMember ; to overwrite
	constexpr explicit operator OtherElement() const {
		return static_cast<OtherElement>(*(this->base_));
	}

	constexpr auto rotated() const& {  // cppcheck-suppress duplInheritedMember ; to overwrite
		typename dynamic_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename dynamic_array::element_const_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated() & {  // cppcheck-suppress duplInheritedMember ; to overwrite
		typename dynamic_array::layout_type new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename dynamic_array::element_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated() && {  // cppcheck-suppress duplInheritedMember ; to overwrite
		typename dynamic_array::layout_type new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename dynamic_array::element_ptr>{new_layout, this->base_};
	}

	friend constexpr auto rotated(dynamic_array& self) -> decltype(auto) { return self.rotated(); }
	friend constexpr auto rotated(dynamic_array const& self) -> decltype(auto) { return self.rotated(); }

 private:
	constexpr auto unrotated_aux_() {
		typename dynamic_array::layout_t new_layout = *this;
		new_layout.unrotate();
		return subarray<T, 0, typename dynamic_array::element_const_ptr>{new_layout, this->base_};
	}

 public:
	// cppcheck-suppress-begin duplInheritedMember ; to overwrite
	constexpr auto unrotated() & { return unrotated_aux_(); }
	constexpr auto unrotated() const& { return unrotated_aux_().as_const(); }
	// cppcheck-suppress-end duplInheritedMember ; to overwrite

	friend constexpr auto unrotated(dynamic_array& self) -> decltype(auto) { return self.unrotated(); }
	friend constexpr auto unrotated(dynamic_array const& self) -> decltype(auto) { return self.unrotated(); }

	constexpr auto operator=(dynamic_array const& other) -> dynamic_array& {
		assert(extents(other) == dynamic_array::extents());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		if(this == &other) {
			return *this;
		}  // lints (cert-oop54-cpp) : handle self-assignment properly
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

 private:
	constexpr auto equal_extensions_if_(std::true_type /*true */, dynamic_array const& other) { return this->extents() == extents(other); }
	constexpr auto equal_extensions_if_(std::false_type /*false*/, dynamic_array const& /*other*/) { return true; }

 public:
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

	constexpr auto operator=(dynamic_array&& other) noexcept -> dynamic_array& {
		assert(equal_extensions_if_(std::integral_constant<bool, (dynamic_array::rank_v != 0)>{}, other));     // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // there is no std::move_n algorithm  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return *this;
	}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

	template<class TT, class... As, class = std::enable_if_t<std::is_assignable<typename dynamic_array::element_ref, TT>{}>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto operator=(dynamic_array<TT, 0, As...> const& other) & -> dynamic_array& {
		assert(extents(other) == dynamic_array::extents());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	constexpr explicit operator subarray<value_type, 0, typename dynamic_array::element_const_ptr, typename dynamic_array::layout_type>() & {  // cppcheck-suppress duplInheritedMember ; to overwrite
		// cppcheck-suppress duplInheritedMember ; to overwrite
		return this->template dynamic_array_cast<value_type, typename dynamic_array::element_const_ptr>();  // cppcheck-suppress duplInheritedMember ; to overwrite
	}

	template<class Archive>
	void serialize(Archive& arxiv, unsigned int const version) {  // cppcheck-suppress duplInheritedMember ; to overwrite
		ref_::serialize(arxiv, version);
	}
};

namespace detail {

template<class T>
struct inplace_array_impl {
	using type = multi::dynamic_array<T, 0, multi::detail::static_allocator<T, 1>>;
};

template<class T, std::size_t N1>
struct inplace_array_impl<T[N1]> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
	using type = multi::dynamic_array<T, 1, multi::detail::static_allocator<T, N1>>;
};

template<class T, std::size_t N1, std::size_t N2>
struct inplace_array_impl<T[N1][N2]> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
	using type = multi::dynamic_array<T, 2, multi::detail::static_allocator<T, (N1 + 1) * N2>>;
};

template<class T, std::size_t N1, std::size_t N2, std::size_t N3>
struct inplace_array_impl<T[N1][N2][N3]> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
	using type = multi::dynamic_array<T, 3, multi::detail::static_allocator<T, (N1 + 1) * N2 * N3>>;
};

template<class T>
struct inplace_array_impl<T[]> {                            // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
	using type = typename inplace_array_impl<T[16]>::type;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
};

template<class T>
struct inplace_array_impl<T*> {
	using type = typename inplace_array_impl<T[16]>::type;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
};

template<class T>
struct inplace_array_impl<T**> {
	using type = typename inplace_array_impl<T[4][4]>::type;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
};

template<class T>
struct inplace_array_impl<T***> {
	using type = typename inplace_array_impl<T[2][2][2]>::type;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for notation
};

}  // namespace detail

/// Stack-allocated multidimensional array whose maximum size is encoded in the element type `T`
///
/// @tparam T C array type encoding the shape, e.g. `double[4][4]`
template<class T>
using inplace_array = typename detail::inplace_array_impl<T>::type;

/// A specialized multidimensional array value with zero dimensions
///
/// @tparam T Element type
/// @tparam Alloc Allocator type
template<typename T, class Alloc>
struct
	array<T, 0, Alloc> : dynamic_array<T, 0, Alloc> {
	using dynamic_array<T, 0, Alloc>::dynamic_array;

	using dynamic_array<T, 0, Alloc>::operator=;

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	BOOST_MULTI_HD constexpr auto operator[]() const& { return *(this->base()); }
#endif

#if !defined(__NVCOMPILER) || (__NVCOMPILER_MAJOR__ > 22 || (__NVCOMPILER_MAJOR__ == 22 && __NVCOMPILER_MINOR__ > 5))  // bug in nvcc 22.5: error: "operator=" has already been declared in the current scope
	template<class TT, class... Args>
	auto operator=(multi::array<TT, 0, Args...> const& other) & -> array& {
		if(other.base()) {
			adl_copy_n(other.base(), other.num_elements(), this->base());
		}
		return *this;
	}

	template<class TT, class... Args>
	auto operator=(multi::array<TT, 0, Args...> const& other) && -> array&& {  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) should assigment return auto& ?
		if(other.base()) {
			adl_copy_n(other.base(), other.num_elements(), this->base());
		}
		return std::move(*this);
	}
#endif

	auto reextent(typename array::extents_type const& /*empty_extensions*/) -> array& {  // NOLINT(readability-redundant-typename)
		return *this;
	}

	// cppcheck-suppress duplInheritedMember ; to overwrite  // NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() && -> array* = delete;  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

template<typename T, ::boost::multi::dimensionality_type D, class Alloc>
class unique_array : public dynamic_array<T, D, Alloc> {
	using dynamic_ = dynamic_array<T, D, Alloc>;

 public:
	using dynamic_array<T, D, Alloc>::dynamic_array;  // NOLINT(cppcoreguidelines-rvalue-reference-param-not-moved,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) false positive on inherited &&-taking constructor

	unique_array(unique_array const&) = default;
	unique_array(unique_array&&)      = default;

	auto operator=(unique_array const&) -> unique_array& = default;
	auto operator=(unique_array&&) -> unique_array&      = default;

#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	~unique_array() noexcept = default;  // pins execution space to match dynamic_array (host-only)
};

template<typename T, ::boost::multi::dimensionality_type D, class Alloc>
struct array : unique_array<T, D, Alloc> {
 private:
	using dynamic_ = dynamic_array<T, D, Alloc>;
	using unique_  = unique_array<T, D, Alloc>;

	static_assert(
		std::is_same_v<typename multi::allocator_traits<Alloc>::value_type, T> || std::is_same_v<typename multi::allocator_traits<Alloc>::value_type, void>,
		"only exact type of array element or void (default) is allowed as allocator value type"
	);

 public:
	// cppcheck-suppress duplInheritedMember ; to override  // NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() && -> array* = delete;  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
	// cppcheck-suppress duplInheritedMember ; to override  // NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() & -> array* { return this; }  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
	// cppcheck-suppress duplInheritedMember ; to override  // NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() const& -> array const* { return this; }  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary

	/// Serializes the array shape and elements into @p arxiv.
	/// On save: writes extensions then elements. On load: reads extensions and resizes the array if needed before loading elements.
	/// Compatible with Boost.Serialization and Cereal.
	template<class Archive, class ArTraits = multi::archive_traits<Archive>>
	void serialize(Archive& arxiv, unsigned int const version) {  // cppcheck-suppress duplInheritedMember ; to override
		auto extents_ = this->extents();

		arxiv& ArTraits::make_nvp("extensions", extents_);  // don't try `using ArTraits::make_nvp`, make_nvp is a static member
		if(this->extents() != extents_) {
			clear();
			this->reextent(extents_);
		}
		dynamic_::serialize(arxiv, version);
	}

	// vvv workaround for MSVC 14.3 and ranges, TODO(correaa) good solution would be to inherit from const_subarray
	BOOST_MULTI_HD operator subarray<T, D, typename array::element_const_ptr, typename array::layout_type> const&() const {     // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
		return reinterpret_cast<subarray<T, D, typename array::element_const_ptr, typename array::layout_type> const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	// move this to dynamic_array
	template<
		class Range, std::enable_if_t<!has_extents<std::decay_t<Range>>::value, int> = 0,
		class = decltype(Range{std::declval<typename array::const_iterator>(), std::declval<typename array::const_iterator>()})>
	constexpr explicit operator Range() const {  // cppcheck-suppress duplInheritedMember ; to overwrite
		// vvv Range{...} needed by Windows GCC?
		return Range{this->begin(), this->end()};  // e.g. std::vector(it, it, alloc = {})
	}

	// TODO(correaa) move this to dynamic_array
	/// Obtain a reference to a C-array reference (`T(&)[N]`) by casting
	template<class CArray, std::enable_if_t<std::is_array_v<CArray>, int> = 0>               // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-use-constraints) for C++20
	constexpr explicit operator CArray&() & { return this->template to_carray_<CArray>(); }  // cppcheck-suppress duplInheritedMember ; to override

	/// Obtain a reference to a constant C-array reference (`T(const&)[N]`) by casting
	template<class CArray, std::enable_if_t<std::is_array_v<CArray>, int> = 0>                          // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-use-constraints) for C++20
	constexpr explicit operator CArray const&() const& { return this->template to_carray_<CArray>(); }  // cppcheck-suppress duplInheritedMember ; to override

	/// Obtain a reference to a C-array reference (`T(&)[N]`) by casting (from an r-value)
	template<class CArray, std::enable_if_t<std::is_array_v<CArray>, int> = 0>                // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,modernize-use-constraints) for C++20
	constexpr explicit operator CArray&() && { return this->template to_carray_<CArray>(); }  // cppcheck-suppress duplInheritedMember ; to override

	// NOLINTNEXTLINE(cppcoreguidelines-rvalue-reference-param-not-moved) false positive in clang-tidy 17-20 ?
	using unique_array<T, D, Alloc>::unique_array;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) passing c-arrays to base
	using typename unique_array<T, D, Alloc>::value_type;

	/// Initializer list constructor from a (nested) list of (subarray) element  @p values. (Nested list should not be ragged.) (allocates)
	template<
		class Sub,
		std::enable_if_t<                                                                                                                                                           // NOLINT(modernize-use-constraints) for C++20
			std::is_constructible_v<typename dynamic_array<T, D>::value_type, Sub> && !std::is_convertible_v<Sub, typename dynamic_array<T, D>::value_type> && (D == 1), int> = 0>  // NOLINT(modernize-use-constraints,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) TODO(correaa) for C++20
	constexpr explicit array(std::initializer_list<Sub> values)                                                                                                                     // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor) inherit explicitness of conversion from the elements
	: unique_(
		  (values.size() == 0) ? array<T, D>()()
							   : array<T, D>(values.begin(), values.end()).element_transformed([](auto const& elem) noexcept -> auto { return static_cast<T>(elem); })
	  ) {}

#ifdef __circle_build__
	constexpr array(std::initializer_list<  // cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
					std::conditional_t<(D != 1), dynamic_array<T, D - 1, typename array::allocator_type>, T>>
						values)
	: unique_(
		  (values.size() == 0) ? array<T, D>{}
							   : array<T, D>(values.begin(), values.end())
	  ) {
	}
#endif

	/// Default constructor of an empty array (doesn't allocate, doesn't throw)
	array() noexcept = default;

	/// Copy constructor from @p other array (generally allocates)
	array(array const&) = default;

	/// Destructor, deallocates memory and destroys elements
#ifdef __NVCC__
#pragma nv_exec_check_disable
#endif
	~array() noexcept = default;

	/// Clear the values of array, making it empty (doesn't throw)
	auto clear() noexcept -> array& {  // cppcheck-suppress duplInheritedMember ; to override
		dynamic_::clear();
		assert(this->stride() != 0);
		return *this;
	}

	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array const& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array&& self) { return std::move(self).data_elements(); }

	// friend BOOST_MULTI_HD constexpr auto move(array& self) -> decltype(auto) { return std::move(self); }
	friend BOOST_MULTI_HD constexpr auto move(array&& self) -> decltype(auto) { return std::move(self); }

	/// Move constructor from @p other array that also sets the allocator @p alloc (may allocate)
	BOOST_MULTI_HD constexpr array(array&& other, Alloc const& alloc) noexcept  ///< Same as the move constructor, except that alloc is used as the allocator.
	: unique_array<T, D, Alloc>{std::move(other), alloc} {}

	BOOST_MULTI_HD constexpr array(array&& other) noexcept : unique_array<T, D, Alloc>{std::move(other)} {
		assert(this->stride() != 0);
	}

	/// Swaps the contents of this array with @p other (doesn't throw)
	void swap(array& other) noexcept {
		using std::swap;
		if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_swap::value) {
			swap(this->alloc(), other.alloc());
		}
		swap(this->base_, other.base_);
		swap(
			this->layout_mutable(),
			other.layout_mutable()
		);
		assert(this->stride() != 0);
	}

#ifndef NOEXCEPT_ASSIGNMENT
	auto operator=(array&& other) noexcept -> array& {
		if(this == std::addressof(other)) {
			return *this;
		}
		clear();
		this->base_ = other.base_;
		if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_move_assignment::value) {
			this->alloc() = std::move(other.alloc());
		}
		this->layout_mutable() = std::exchange(other.layout_mutable(), typename array::layout_type(typename array::extents_type{}));
		return *this;
	}

	/// Copy assignment from @p other array (allocates unless extents are already equal)
	auto operator=(array const& other) -> array& {
		if(this == &other) {  // required by cert-oop54-cpp
			return *this;
		}
		if(array::extents() == other.extents()) {
			if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			dynamic_::operator=(other);
		} else {
			clear();
			if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			this->layout_mutable() = other.layout();
			array::allocate();
			array::uninitialized_copy_elements(other.data_elements());
		}
		return *this;
	}
#else
	auto operator=(array o) noexcept -> array& { return swap(o), *this; }
#endif

	template<typename OtherT, typename OtherEP, class OtherLayout>
	auto operator=(multi::const_subarray<OtherT, D, OtherEP, OtherLayout> const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			dynamic_::operator=(other);  // TODO(correaa) : protect for self assigment
		} else {
			operator=(array{other});
		}
		return *this;
	}

	/// Assignment from @p other array of a different element type (allocates, unless extents are already equal or the number of elements is the same)
	template<class TT, class AAlloc>
	auto operator=(multi::array<TT, D, AAlloc> const& other) -> array& {
		if(array::extents() == other.extents()) {
			dynamic_::operator=(other);
		} else if(this->num_elements() == other.extents().num_elements()) {
			this->layout_mutable() = typename array::layout_type(other.extents());
			dynamic_::operator=(other);
		} else {
			operator=(static_cast<array>(other));
		}
		return *this;
	}

	/// Assignment from @p other multdimensional range (allocates, unless extents are already equal or the number of elements is the same)
	template<
		class Range, class = decltype(std::declval<dynamic_&>().operator=(std::declval<Range&&>())),
		std::enable_if_t<!has_data_elements<std::decay_t<Range>>::value, int> = 0,
		std::enable_if_t<has_extents<std::decay_t<Range>>::value, int>        = 0,
		std::enable_if_t<!std::is_base_of_v<array, std::decay_t<Range>>, int> = 0>  // NOLINT(modernize-use-constraints,modernize-type-traits) for C++20
	auto operator=(Range&& other) -> array& {
		if(array::extents() == other.extents()) {
			this->operator()() = std::forward<Range>(other);
		} else if(this->num_elements() == other.extents().num_elements()) {
			typename array::layout_type const new_layout{other.extents()};
			assert(new_layout.num_elements() == this->num_elements());
			this->layout_mutable() = new_layout;
			assert(this->stride() != 0);
			// return *this;
			// reshape(other.extents());
			this->operator()() = std::forward<Range>(other);
		} else {
			operator=(static_cast<array>(std::forward<Range>(other)));
		}
		return *this;
	}

	template<
		class Range, class = decltype(std::declval<dynamic_&>().operator=(std::declval<Range&&>())),
		std::enable_if_t<!std::is_base_of_v<array, std::decay_t<Range>>, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto from(Range&& other) -> array& {                                            // TODO(correaa) : check that LHS is not read-only?
		if(array::extents() == other.extents()) {
			this->operator()() = other;
		} else if(this->num_elements() == other.extents().num_elements()) {
			reshape(other.extents());
			this->operator()() = other;
		} else {
			operator=(static_cast<array>(std::forward<Range>(other)));
		}
		return *this;
	}

	friend void swap(array& self, array& other) noexcept(true /*noexcept(self.swap(other))*/) { self.swap(other); }

	void assign(typename array::extents_type extensions, typename array::element const& elem) {
		if(array::extents() == extensions) {
			adl_fill_n(this->base_, this->num_elements(), elem);
		} else {
			this->clear();
			(*this).array::layout_t::operator=(layout_t<D>{extensions});
			this->base_ = this->dynamic_::array_alloc::allocate(this->num_elements(), nullptr);
			adl_alloc_uninitialized_fill_n(this->alloc(), this->base_, this->num_elements(), elem);
		}
	}

	/// Assigns elements from an iterator range [`first`, `last`), resizing if necessary. complexity: O(n)
	template<class It>
	void assign(It first, It last) {  // cppcheck-suppress duplInheritedMember ; to overwrite
		using std::all_of;
		using std::next;
		if(adl_distance(first, last) == this->size()) {
			dynamic_::ref_::assign(first);
		} else {
			this->operator=(array(first, last));
		}
	}

	void assign(std::initializer_list<value_type> values) {
		if(values.size() != 0) {
			assign(values.begin(), values.end());
		}
	}

	template<class Range> auto assign(Range&& other) & -> decltype(assign(adl_begin(std::forward<Range>(other)), adl_end(std::forward<Range>(other)))) {
		return assign(adl_begin(std::forward<Range>(other)), adl_end(std::forward<Range>(other)));  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
	}

	// Assignment from a (nested) list of (subarray) element  @p values. (Nested list should not be ragged.) (Allocates unless extents match)
	auto operator=(std::initializer_list<value_type> values) -> array& {
		if(values.size() == 0) {
			this->clear();
		} else {
			assign(values.begin(), values.end());
		}
		return *this;
	}

	/// Change the extents of the array to @p exts, preserving elements when possible. (generally allocates, elements are discarded unless extents do not change).
	auto reextent(typename array::extents_type const& exts) && -> array&& {  // NOLINT(readability-redundant-typename)
		if(exts == this->extents()) {
			return std::move(*this);
		}

		auto new_layout = typename array::layout_type{exts};

		if(new_layout.num_elements() != this->layout().num_elements()) {
			this->destroy();
			this->deallocate();

			this->layout_mutable() = new_layout;  // typename array::layout_t{extensions};

			this->base_ = this->dynamic_::array_alloc::allocate(
				static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(
					new_layout.num_elements()
				),
				this->data_elements()  // used as hint
			);

			if constexpr(!(std::is_trivially_default_constructible_v<typename array::element> || multi::force_element_trivial_default_construction<typename array::element>)) {
				adl_alloc_uninitialized_value_construct_n(this->alloc(), this->base_, this->num_elements());
			}
		} else {
			this->layout_mutable() = new_layout;
		}

		return std::move(*this);
	}

	auto reextent(typename array::extents_type const& extensions) & -> array& {  // NOLINT(readability-redundant-typename)
		if(extensions == this->extents()) {
			return *this;
		}
		auto&& tmp = typename array::ref_(
			this->dynamic_::array_alloc::allocate(
				static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(
					typename array::layout_type{extensions}.num_elements()
				),
				this->data_elements()  // used as hint
			),
			extensions
		);
		if constexpr(!(std::is_trivially_default_constructible_v<typename array::element> || multi::force_element_trivial_default_construction<typename array::element>)) {
			adl_alloc_uninitialized_value_construct_n(this->alloc(), tmp.data_elements(), tmp.num_elements());
		}

		auto const intersect = intersection(this->extents(), extensions);

		tmp.apply(intersect) = this->apply(intersect);  // TODO(correaa) : use `.moved_elements()`? or move_n?

		this->destroy();
		this->deallocate();

		this->base_            = tmp.base();
		this->layout_mutable() = tmp.layout();

		return *this;
	}

	[[nodiscard]] constexpr auto operator+() const& { return array{*this}; }         // cppcheck-suppress duplInheritedMember ; to overwrite
	[[nodiscard]] constexpr auto operator+() && { return array{std::move(*this)}; }  // cppcheck-suppress duplInheritedMember ; to overwrite  //-V::659 move version

	auto reextent(typename array::extents_type const& exs, typename array::element const& elem) & -> array& {  // NOLINT(readability-redundant-typename)
		if(exs == this->extents()) {
			return *this;
		}

		// implementation with hint
		auto&& tmp = typename array::ref_(
			this->dynamic_::array_alloc::allocate(
				static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(typename array::layout_type{exs}.num_elements()),
				this->data_elements()  // use as hint
			),
			exs
		);
		this->uninitialized_fill_n(tmp.data_elements(), static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(tmp.num_elements()), elem);
		auto const intersect = intersection(this->extents(), exs);
		tmp.apply(intersect) = this->apply(intersect);  // TODO(correaa) use moved_elements?
		this->destroy();
		this->deallocate();
		this->base_            = tmp.base();
		this->layout_mutable() = tmp.layout();

		return *this;
	}
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __cpp_deduction_guides

#define BOOST_MULTI_IL std::initializer_list  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing, TODO(correaa) use template typedef instead of macro

// vvv MSVC 14.3 in c++17 mode needs paranthesis in dimensionality_type(d)
// template<class T> dynamic_array(BOOST_MULTI_IL<T>) -> dynamic_array<T, static_cast<dimensionality_type>(1U), std::allocator<T>>;  // MSVC needs the allocator argument error C2955: 'boost::multi::dynamic_array': use of class template requires template argument list
// template<class T> dynamic_array(BOOST_MULTI_IL<BOOST_MULTI_IL<T>>) -> dynamic_array<T, static_cast<dimensionality_type>(2U), std::allocator<T>>;
// template<class T> dynamic_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>) -> dynamic_array<T, static_cast<dimensionality_type>(3U), std::allocator<T>>;
// template<class T> dynamic_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>) -> dynamic_array<T, static_cast<dimensionality_type>(4U), std::allocator<T>>;
// template<class T> dynamic_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>>) -> dynamic_array<T, static_cast<dimensionality_type>(5U), std::allocator<T>>;

// TODO(correaa) add zero dimensional case?
// template<class T> array(BOOST_MULTI_IL<T>) -> array<T, static_cast<dimensionality_type>(1U)>;
// template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<T>>) -> array<T, static_cast<dimensionality_type>(2U)>;
// template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>) -> array<T, static_cast<dimensionality_type>(3U)>;
// template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>) -> array<T, static_cast<dimensionality_type>(4U)>;
// template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>>) -> array<T, static_cast<dimensionality_type>(5U)>;

#undef BOOST_MULTI_IL

template<class T> array(T[]) -> array<T, static_cast<dimensionality_type>(1U)>;  // NOSONAR(cpp:S5945) NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

//  vvv these are necessary to catch {n, m, ...} notation (or single integer notation)
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<0>, T) -> array<T, static_cast<dimensionality_type>(0U)>;  // TODO(correaa) use some std::allocator_traits instead of is_allocator
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<1>, T) -> array<T, static_cast<dimensionality_type>(1U)>;
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<2>, T) -> array<T, static_cast<dimensionality_type>(2U)>;
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<3>, T) -> array<T, static_cast<dimensionality_type>(3U)>;
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<4>, T) -> array<T, static_cast<dimensionality_type>(4U)>;
template<class T, class = std::enable_if_t<!multi::is_allocator_v<T>>> array(iextensions<5>, T) -> array<T, static_cast<dimensionality_type>(5U)>;

// generalization, will not work with naked {n, m, ...} notation (or single integer notation)
template<dimensionality_type D, class T, class = std::enable_if_t<!boost::multi::is_allocator_v<T>>>
array(iextensions<D>, T) -> array<T, D>;

template<class MatrixRef, class DT = typename MatrixRef::decay_type, class T = typename DT::element, dimensionality_type D = DT::rank_v, class Alloc = typename DT::allocator_type>
array(MatrixRef) -> array<T, D, Alloc>;

template<class MatValues, class T = typename MatValues::element, dimensionality_type D = MatValues::rank_v>
array(MatValues) -> array<T, D>;

template<class T, dimensionality_type D>
array(const_subarray<T, D>) -> array<T, D>;

template<class MatValues, class T = typename MatValues::element, dimensionality_type D = MatValues::rank_v, class Alloc = std::allocator<T>, class = std::enable_if_t<multi::is_allocator_v<Alloc>>>  /// , class Alloc = typename DT::allocator_type>
array(MatValues, Alloc) -> array<T, D, Alloc>;

template<typename T, dimensionality_type D, typename P> array(subarray<T, D, P>) -> array<T, D>;

template<
	class Range, std::enable_if_t<!has_extents<Range>::value, int> = 0,
	typename V = decltype(*::std::begin(std::declval<Range const&>()))
	// typename V = typename std::iterator_traits<decltype(::std::begin(std::declval<Range const&>()))>::value_type
	>
array(Range) -> array<V, 1>;

template<class Reference>
auto operator+(Reference const& ref) -> decltype(array<typename Reference::element, Reference::dimensionality>(ref)) {
	return array<typename Reference::element, Reference::dimensionality>(ref);
}

#endif  // ends defined(__cpp_deduction_guides)

// template<class T, std::size_t N>
// // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
// auto decay(T const (&arr)[N]) noexcept -> multi::array<std::remove_all_extents_t<T[N]>, std::rank_v<T[N]>> {
// 	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
// 	return multi::array_cref<std::remove_all_extents_t<T[N]>, std::rank_v<T[N]>>(data_elements(arr), extents(arr));
// }

template<class T, std::size_t N>
struct detail::array_traits<T[N], void, void> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using reference  = T&;
	using element    = std::remove_all_extents_t<T[N]>;  // NOSONAR(cpp:S5945) NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using decay_type = multi::array<T, 1>;
};

}  // end namespace boost::multi

/// Convenience aliases using `std::pmr::polymorphic_allocator` as the allocator.
namespace boost::multi::pmr {

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
template<class T, boost::multi::dimensionality_type D> using array         = ::boost::multi::array<T, D, std::pmr::polymorphic_allocator<T>>;
template<class T, boost::multi::dimensionality_type D> using dynamic_array = ::boost::multi::dynamic_array<T, D, std::pmr::polymorphic_allocator<T>>;
#else
template<class T, boost::multi::dimensionality_type D> struct [[deprecated("no PMR allocator")]] array;          // your version of C++ doesn't provide polymorphic_allocators
template<class T, boost::multi::dimensionality_type D> struct [[deprecated("no PMR allocator")]] dynamic_array;  // your version of C++ doesn't provide polymorphic_allocators
#endif

}  // end namespace boost::multi::pmr

namespace boost::serialization {

template<typename T> struct version;  // in case serialization was not included before

template<typename T, boost::multi::dimensionality_type D, class A>
struct version<boost::multi::array<T, D, A>> {
	using type = std::integral_constant<int, BOOST_MULTI_SERIALIZATION_ARRAY_VERSION>;  // TODO(correaa) use constexpr variable here, not a macro
	// NOLINTNEXTLINE(cppcoreguidelines-use-enum-class) for backward compatibility with Boost Serialization
	enum /*class value_t*/ { value = type::value };  // NOSONAR(cpp:S3642)  // https://community.sonarsource.com/t/suppress-issue-in-c-source-file/43154/24
};

}  // end namespace boost::serialization

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_ARRAY_HPP
