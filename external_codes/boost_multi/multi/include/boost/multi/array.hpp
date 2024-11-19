// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ARRAY_HPP_
#define BOOST_MULTI_ARRAY_HPP_

#include <boost/multi/array_ref.hpp>  // IWYU pragma: export

#include <boost/multi/detail/config/NO_UNIQUE_ADDRESS.hpp>

#include <boost/multi/detail/adl.hpp>
#include <boost/multi/detail/is_trivial.hpp>
#include <boost/multi/detail/memory.hpp>

#include <memory>  // for std::allocator_traits
#include <tuple>  // needed by a deprecated function
#include <type_traits>  // for std::common_reference
#include <utility>  // for std::move

#if __has_include(<memory_resource>)
#  include <memory_resource>
// Apple clang provides the header but not the compiled library prior to version 16
#  if (defined(__cpp_lib_memory_resource) && (__cpp_lib_memory_resource >= 201603)) && !(defined(__APPLE__) && defined(__clang_major__) && __clang_major__ <= 15) && (!defined(_LIBCPP_VERSION) || !(_LIBCPP_VERSION <= 160001) )
#    define BOOST_MULTI_HAS_MEMORY_RESOURCE
#  endif
#endif

// TODO(correaa) or should be (__CUDA__) or CUDA__ || HIP__
#if defined(__NVCC__)
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif

namespace boost::multi {

namespace detail {

template<class Allocator>
struct array_allocator {
	using allocator_type = Allocator;
	array_allocator()    = default;

 private:
	BOOST_MULTI_NO_UNIQUE_ADDRESS allocator_type alloc_;

	using allocator_traits = typename multi::allocator_traits<allocator_type>;
	using size_type_       = typename allocator_traits::size_type;
	using pointer_         = typename allocator_traits::pointer;

 protected:
	constexpr auto alloc() & -> auto& { return alloc_; }
	constexpr auto alloc() const& -> allocator_type const& { return alloc_; }

	constexpr explicit array_allocator(allocator_type const& alloc) : alloc_{alloc} {}  // NOLINT(modernize-pass-by-value)

	constexpr auto allocate(size_type_ n) -> pointer_ {
		return n ? allocator_traits::allocate(alloc_, n) : pointer_{nullptr};
	}
	constexpr auto allocate(size_type_ n, typename allocator_traits::const_void_pointer hint) -> pointer_ {
		return n ? allocator_traits::allocate(alloc_, n, hint) : pointer_{nullptr};
	}

	constexpr auto uninitialized_fill_n(pointer_ first, size_type_ count, typename allocator_traits::value_type const& value) {
		return adl_alloc_uninitialized_fill_n(alloc_, first, count, value);
	}

	template<typename It>
	auto uninitialized_copy_n(It first, size_type count, pointer_ d_first) {
		#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename std::pointer_traits<pointer_>::element_type> && ! multi::force_element_trivial_default_construction<typename std::pointer_traits<pointer_>::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(alloc_, d_first, count);
		}
		return adl_copy_n                    (        first, count, d_first);
		#else
		return adl_alloc_uninitialized_copy_n(alloc_, first, count, d_first);
		#endif
	}

	template<typename It>
	auto uninitialized_move_n(It first, size_type count, pointer_ d_first) {
		#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename std::pointer_traits<pointer_>::element_type> && ! multi::force_element_trivial_default_construction<typename std::pointer_traits<pointer_>::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(alloc_, d_first, count);
		}
		return                     adl_copy_n(        std::make_move_iterator(first), count, d_first);
		#else
		return adl_alloc_uninitialized_move_n(alloc_,                         first , count, d_first);
		#endif
	}

	template<class EP, typename It>
	auto uninitialized_copy_n(EP&& ep, It first, size_type count, pointer_ d_first) {
		// #if defined(__clang__) && defined(__CUDACC__)
		// if constexpr(! std::is_trivially_default_constructible_v<typename std::pointer_traits<pointer_>::element_type> && ! multi::force_element_trivial_default_construction<typename std::pointer_traits<pointer_>::element_type> ) {
		//  adl_alloc_uninitialized_default_construct_n(alloc_, d_first, count);
		// }
		// return adl_copy_n                    (        first, count, d_first);
		// #else
		return adl_uninitialized_copy_n(std::forward<EP>(ep), first, count, d_first);
		// return adl_alloc_uninitialized_copy_n(std::forward<EP>(ep), alloc_, first, count, d_first);
		// #endif
	}

	template<typename It>
	auto destroy_n(It first, size_type n) { return adl_alloc_destroy_n(this->alloc(), first, n); }

 public:
	constexpr auto get_allocator() const -> allocator_type { return alloc_; }
};
}  // end namespace detail

template<class T, dimensionality_type D, class DummyAlloc = std::allocator<T>>  // DummyAlloc mechanism allows using the convention array<T, an_allocator<>>, is an_allocator supports void template argument
struct static_array  // NOLINT(fuchsia-multiple-inheritance) : multiple inheritance used for composition
: protected detail::array_allocator<
	  // Alloc
	  typename allocator_traits<DummyAlloc>::template rebind_alloc<T>>
, public array_ref<T, D, typename multi::allocator_traits<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::pointer>
, boost::multi::random_iterable<static_array<T, D, typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>> {
	static_assert(
		std::is_same_v<
			std::remove_const_t<typename multi::allocator_traits<DummyAlloc>::value_type>,
			typename static_array::element_type
		> 
		||
		std::is_same_v<
			std::remove_const_t<typename multi::allocator_traits<DummyAlloc>::value_type>,
			void
		>,  // allocator template can be redundant or void (which can be a default for the allocator)
		"allocator value type must match array value type"
	);

 private:
	// using Alloc = typename allocator_traits<DummyAlloc>::template rebind_alloc<T>;

 protected:
	using array_alloc = detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T> >;

 public:
	// constexpr auto get_allocator() const -> allocator_type { return alloc_; }
	using detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T> >::get_allocator;

	using allocator_type = typename detail::array_allocator<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::allocator_type;
	using decay_type     = array<T, D, allocator_type>;
	using layout_type    = typename array_ref<T, D, typename multi::allocator_traits<allocator_type>::pointer>::layout_type;

	using ref = array_ref<
		T, D,
		typename multi::allocator_traits<typename multi::allocator_traits<allocator_type>::template rebind_alloc<T>>::pointer
	>;

	auto operator new(std::size_t count) -> void* { return ::operator new(count); }
	auto operator new(std::size_t count, void* ptr) -> void* { return ::operator new(count, ptr); }
	void operator delete(void* ptr) noexcept { ::operator delete(ptr); }  // this overrides the deleted delete operator in reference (base) class subarray

 protected:
	using alloc_traits = typename multi::allocator_traits<allocator_type>;

	auto uninitialized_value_construct() {
		return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}

	constexpr auto uninitialized_default_construct() {
		if constexpr(!std::is_trivially_default_constructible_v<typename static_array::element_type> &&  ! multi::force_element_trivial_default_construction<typename static_array::element_type>) {
			return adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->base_, this->num_elements());
		}
	}

	template<typename It> auto uninitialized_copy_elements(It first) {
		return array_alloc::uninitialized_copy_n(first, this->num_elements(), this->data_elements());
	}

	template<typename It> auto uninitialized_move_elements(It first) {
		return array_alloc::uninitialized_move_n(first, this->num_elements(), this->data_elements());
	}

	template<class EP, typename It> auto uninitialized_copy_elements(EP&& ep, It first) {
		return array_alloc::uninitialized_copy_n(std::forward<EP>(ep), first, this->num_elements(), this->data_elements());
	}

	constexpr void destroy() {
		if constexpr(!(std::is_trivially_destructible_v<typename static_array::element_type> || multi::force_element_trivial_destruction<typename static_array::element_type>)) {
			array_alloc::destroy_n(this->data_elements(), this->num_elements());
		}
	}

	void allocate() {
		this->base_ = array_alloc::allocate(static_cast<typename multi::allocator_traits<typename static_array::allocator_type>::size_type>(this->static_array::num_elements()));
	}

 public:
	using value_type = typename std::conditional_t<
		(D > 1),  // this parenthesis is needed
		array<typename static_array::element_type, D - 1, allocator_type>,
		typename static_array::element_type>;

	using typename ref::difference_type;
	using typename ref::size_type;
	explicit static_array(allocator_type const& alloc) : array_alloc{alloc}, ref(nullptr, {}) {}

	using ref::       operator();
	BOOST_MULTI_HD constexpr auto operator()() && -> decltype(auto) { return ref::element_moved(); }

	using ref::taked;
	constexpr auto taked(difference_type n) && -> decltype(auto) { return ref::taked(n).element_moved(); }

	using ref::dropped;
	constexpr auto dropped(difference_type n) && -> decltype(auto) { return ref::dropped(n).element_moved(); }

	static_array(static_array&& other) noexcept :
		array_alloc{other.alloc()},
		ref{
			array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),
			other.extensions()
		}
	{
		uninitialized_move_elements(other.data_elements());
	}

	constexpr static_array(decay_type&& other, allocator_type const& alloc) noexcept
	: array_alloc{alloc}, ref(std::exchange(other.base_, nullptr), other.extensions()) {
		std::move(other).layout_mutable() = typename static_array::layout_type(typename static_array::extensions_type{});  // = {};  careful! this is the place where layout can become invalid
	}

	constexpr explicit static_array(decay_type&& other) noexcept
	: static_array(std::move(other), allocator_type{}) {}  // 6b

	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>
	constexpr explicit static_array(It first, It last, allocator_type const& alloc)
	: array_alloc{alloc}, ref{
		                      array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(layout_type{index_extension(adl_distance(first, last)) * multi::extensions(*first)}.num_elements())),
		                      index_extension(adl_distance(first, last)) * multi::extensions(*first)} {
	#if defined(__clang__) && defined(__CUDACC__)
		// TODO(correaa) add workaround for non-default constructible type and use adl_alloc_uninitialized_default_construct_n
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(static_array::alloc(), ref::data_elements(), ref::num_elements());
		}
		adl_copy_n(first, last - first, ref::begin());
	#else
		adl_alloc_uninitialized_copy(static_array::alloc(), first, last, ref::begin());
	#endif
	}

	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>
	constexpr explicit static_array(It first, It last) : static_array(first, last, allocator_type{}) {}

	template<
		class Range, class = std::enable_if_t<!std::is_base_of<static_array, std::decay_t<Range>>{}>,
		class = decltype(/*static_array*/ (std::declval<Range const&>().begin() - std::declval<Range const&>().end())),  // instantiation of static_array here gives a compiler error in 11.0, partially defined type?
		class = std::enable_if_t<! is_subarray<Range const&>::value>  // NOLINT(modernize-use-constraints) TODO(correaa) in C++20
	>
	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions // NOLINTNEXTLINE(runtime/explicit)
	static_array(Range const& rng)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array{std::begin(rng), std::end(rng)} {}  // Sonar: Prefer free functions over member functions when handling objects of generic type "Range".
	// : static_array{rng.begin(), rng.end()} {}   // Sonar: Prefer free functions over member functions when handling objects of generic type "Range".

	template<class TT>
	auto uninitialized_fill_elements(TT const& value) {
		return array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), value);
	}

	template<class TT, class... As>
	static_array(array_ref<TT, D, As...> const& other, allocator_type const& alloc)
	: array_alloc{alloc}, ref{
		                      array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),
		                      other.extensions()} {
	#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
	#else
		adl_alloc_uninitialized_copy_n(static_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
	#endif
	}

	static_array(typename static_array::extensions_type extensions, typename static_array::element_type const& elem, allocator_type const& alloc)  // 2
	: array_alloc{alloc}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), nullptr), extensions} {
		array_alloc::uninitialized_fill_n(this->data_elements(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);
	}

	// template<class... Exts, class... Ts>
	// explicit static_array(std::tuple<Exts...> extensions, Ts&&... args)  // this is important to pass arguments to boost::interprocess::construct
	// : static_array{
	//  std::apply([](auto... exts) {return typename static_array::extensions_type{exts...};}, extensions),
	//  std::forward<Ts>(args)...
	// } {}

	template<class Element>
	explicit static_array(
		Element const& elem, allocator_type const& alloc,
		std::enable_if_t<std::is_convertible_v<Element, typename static_array::element_type> && (D == 0), int> /*dummy*/ = 0  // NOLINT(fuchsia-default-arguments-declarations) for classic sfinae, needed by MSVC?
	)
	: static_array(typename static_array::extensions_type{}, elem, alloc) {}

	constexpr static_array(typename static_array::extensions_type exts, typename static_array::element_type const& elem)
	:
	array_alloc{},
	array_ref<T, D, typename multi::allocator_traits<typename multi::allocator_traits<DummyAlloc>::template rebind_alloc<T>>::pointer>(
		exts,
		array_alloc::allocate(
			static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t(exts).num_elements()) ,
			nullptr
		)
	) {
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type>) {
			array_alloc::uninitialized_fill_n(this->base(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);
		} else {  // this workaround allows constexpr arrays for simple types
		                           adl_fill_n(this->base(), static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);
		}
	}

	template<class ValueType, class = decltype(std::declval<ValueType>().extensions()),
		std::enable_if_t<std::is_convertible_v<ValueType, typename static_array::value_type>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)  // fill constructor
	: array_alloc{alloc}, ref(array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t(extension*value.extensions()).num_elements())), extension*value.extensions())
	{
		static_assert(std::is_trivially_default_constructible_v<typename static_array::element_type> || multi::force_element_trivial_default_construction<typename static_array::element_type> );  // TODO(correaa) not implemented for non-trivial types,
		adl_fill_n(this->begin(), this->size(), value);  // TODO(correaa) implement via .elements()? substitute with uninitialized version of fill, uninitialized_fill_n?
	}

	template<class ValueType, class = decltype(std::declval<ValueType>().extensions()),
		std::enable_if_t<std::is_convertible_v<ValueType, typename static_array::value_type>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value)  // fill constructor
	: static_array(extension, value, allocator_type{}) {}

	constexpr explicit static_array(typename static_array::extensions_type extensions, allocator_type const& alloc)
	: array_alloc{alloc}, ref(array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements())), extensions) {
		uninitialized_default_construct();
		assert(this->stride() != 0);
	}

	constexpr explicit static_array(typename static_array::extensions_type extensions)
	: static_array(extensions, allocator_type{}) {}

	template<class OtherT, class OtherEP, class OtherLayout,  // class... Args,
	         class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::subarray<OtherT, D, OtherEP, OtherLayout>::element_type>{}>,
	         class = decltype(adl_copy(std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().begin(), std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().end(), std::declval<typename static_array::iterator>()))>
	constexpr static_array(multi::const_subarray<OtherT, D, OtherEP, OtherLayout> const& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{other.extensions()}.num_elements())),
		  other.extensions()
	  ) {
		// #if (defined(__clang__) && defined(__CUDACC__))
		// if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
		//  adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->data_elements(), this->num_elements());
		// }
		// adl_copy_n                    (                       other.elements().begin(), this->num_elements(), this->data_elements());
		// #else
		adl_alloc_uninitialized_copy_n(static_array::alloc(), other.elements().begin(), this->num_elements(), this->data_elements());
		// #endif
	}

	template<class OtherT, class OtherEP, class OtherLayout,  // class... Args,
	         class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::subarray<OtherT, D, OtherEP, OtherLayout>::element_type>{}>,
	         class = decltype(adl_copy(std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().begin(), std::declval<multi::subarray<OtherT, D, OtherEP, OtherLayout> const&>().end(), std::declval<typename static_array::iterator>()))>
	constexpr static_array(multi::subarray<OtherT, D, OtherEP, OtherLayout>&& other, allocator_type const& alloc)
	: array_alloc{alloc},
	  ref(
		  array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{other.extensions()}.num_elements())),
		  other.extensions()
	  ) {
		adl_alloc_uninitialized_copy_n(static_array::alloc(), std::move(other).elements().begin(), this->num_elements(), this->data_elements());
	}

	template<
		class TT, class EElementPtr, class LLayout,
		std::enable_if_t< ! multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout>&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	explicit static_array(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other)
	: static_array(other, allocator_type{}) {}

	template<
		class TT, class EElementPtr, class LLayout,
		std::enable_if_t<   multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	// cppcheck-suppress noExplicitConstructor  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ static_array(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: static_array(other, allocator_type{}) {}

	template<
		class TT, class EElementPtr, class LLayout,
		std::enable_if_t<   multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().base()), T>, int> = 0,
		class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::const_subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	// cppcheck-suppress noExplicitConstructor  // NOLINTNEXTLINE(runtime/explicit)
	constexpr /*implicit*/ static_array(multi::subarray<TT, D, EElementPtr, LLayout>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: static_array(std::move(other), allocator_type{}) {}

	constexpr static_array(multi::subarray<T, D, typename static_array::element_ptr, typename static_array::layout_type> const&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: static_array(other, allocator_type{}) {}

	constexpr static_array(multi::const_subarray<T, D, typename static_array::element_ptr, typename static_array::layout_type> const&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: static_array(other, allocator_type{}) {}

	constexpr static_array(multi::subarray<T, D, typename static_array::element_ptr, typename static_array::layout_type>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: static_array(std::move(other), allocator_type{}) {}

	// template<class TT, class EElementPtr, class LLayout,  // class... Args,
	//          std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().base()), T>, int> = 0,
	//          class                                                                                                                                           = decltype(adl_copy(std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().begin(), std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename static_array::iterator>()))>
	// explicit static_array(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other)
	// : static_array(other, allocator_type{}) {}



	// template<typename TT, typename EElementPtr, class LLayout,
	//          std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::subarray<TT, D, EElementPtr, LLayout>&&>().base()), T>, int> = 0,
	//          class                                                                                                                                     = decltype(adl_copy(std::declval<multi::subarray<TT, D, EElementPtr, LLayout>&&>().begin(), std::declval<multi::subarray<TT, D, EElementPtr, LLayout> const&>().end(), std::declval<typename static_array::iterator>()))>
	// // cppcheck-suppress noExplicitConstructor ; // NOLINTNEXTLINE(runtime/explicit)
	// /*mplct*/ static_array(multi::subarray<TT, D, EElementPtr, LLayout>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	// : static_array(std::move(other), allocator_type{}) {}

	// template<typename TT, typename... Args,
	//          std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<multi::const_subarray<TT, D, Args...>&&>().base()), T>, int> = 0,
	//          class = decltype(adl_copy(std::declval<multi::const_subarray<TT, D, Args...>&&>().begin(), std::declval<multi::const_subarray<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))>
	// explicit static_array(multi::const_subarray<TT, D, Args...>&& other)
	// : static_array(std::move(other), allocator_type{}) {}

	template<class TT, class... Args,
		std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ static_array(array_ref<TT, D, Args...>& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: array_alloc{}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extensions()} {
		static_array::uninitialized_copy_elements(other.data_elements());
	}

	template<class TT, class... Args,
		std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(array_ref<TT, D, Args...>& other)  // NOLINT(fuchsia-default-arguments-declarations)
	: array_alloc{}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extensions()} {
		assert(this->stride() != 0);
		static_array::uninitialized_copy_elements(other.data_elements());
	}

	template<class TT, class... Args,
		std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ static_array(array_ref<TT, D, Args...>&& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: array_alloc{}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extensions()} {
		assert(this->stride() != 0);
		static_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	template<class TT, class... Args,
		std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...>&&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(array_ref<TT, D, Args...>&& other)  // NOLINT(fuchsia-default-arguments-declarations)
	: array_alloc{}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extensions()} {
		assert(this->stride() != 0);
		static_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	template<class TT, class... Args,
		std::enable_if_t<multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...> const&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax
	/*mplct*/ static_array(array_ref<TT, D, Args...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: array_alloc{}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())), other.extensions()} {
		assert(this->stride() != 0);
		static_array::uninitialized_copy_elements(other.data_elements());
	}

	template<class TT, class... Args,
		std::enable_if_t<!multi::detail::is_implicitly_convertible_v<decltype(*std::declval<array_ref<TT, D, Args...> const&>().base()), T>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(array_ref<TT, D, Args...> const& other)  // NOLINT(fuchsia-default-arguments-declarations)
	:
	array_alloc{},
	ref(
		array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements())),
		other.extensions()
	) {
		assert(this->stride() != 0);
		static_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	static_array(static_array const& other)  // 5b
	: 
		array_alloc{
			multi::allocator_traits<allocator_type>::select_on_container_copy_construction(other.alloc())
		},
		ref{
			array_alloc::allocate(
				static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements()) //,
				// other.data_elements()
			),
			other.extensions()
		}
	{
		assert(this->stride() != 0);
		uninitialized_copy_elements(other.data_elements());
	}

	template<class ExecutionPolicy, std::enable_if_t<!std::is_convertible_v<ExecutionPolicy, typename static_array::extensions_type>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	static_array(ExecutionPolicy&& policy, static_array const& other)
	: array_alloc{multi::allocator_traits<allocator_type>::select_on_container_copy_construction(other.alloc())}, ref{array_alloc::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), extensions(other)} {
		assert(this->stride() != 0);
		uninitialized_copy_elements(std::forward<ExecutionPolicy>(policy), other.data_elements());
	}

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	constexpr static_array(std::initializer_list<typename static_array<T, D>::value_type> values)
	: static_array{array<T, D>(values.begin(), values.end())} {}  // construct all with default constructor and copy to special memory at the end

	static_array(
		std::initializer_list<typename static_array<T, D>::value_type> values,
		allocator_type const&                                          alloc
	)
	: static_array{static_array<T, D>(values.begin(), values.end()), alloc} {
		assert(this->stride() != 0);
	}

	template<class TT, std::size_t N>
	constexpr explicit static_array(TT (&array)[N])  // @SuppressWarnings(cpp:S5945) NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backward compatibility // NOSONAR
	: static_array(std::begin(array), std::end(array)) {
		assert(this->stride() != 0);
	}

	constexpr auto begin() const& -> typename static_array::const_iterator { return ref::begin(); }
	constexpr auto end() const& -> typename static_array::const_iterator { return ref::end(); }

	constexpr auto begin() && -> typename static_array::move_iterator { return ref::begin(); }
	constexpr auto end() && -> typename static_array::move_iterator { return ref::end(); }

	constexpr auto begin() & -> typename static_array::iterator { return ref::begin(); }
	constexpr auto end() & -> typename static_array::iterator { return ref::end(); }

	using ref::operator[];
	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> decltype(auto) {
		return multi::move(ref::operator[](idx));
		// if constexpr(D == 1) {
		//  return std::move(ref::operator[](idx));
		// } else {
		//  return ref::operator[](idx).element_moved();
		// }  // NOLINT(readability/braces)
	}

	constexpr auto max_size() const noexcept { return static_cast<typename static_array::size_type>(multi::allocator_traits<allocator_type>::max_size(this->alloc())); }  // TODO(correaa)  divide by nelements in under-dimensions?

 protected:
	constexpr void deallocate() {
		assert(this->stride() != 0);
		if(this->num_elements()) {
			multi::allocator_traits<allocator_type>::deallocate(this->alloc(), this->base_, static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()));
		}
	}
	void clear() noexcept {
		this->destroy();
		deallocate();
		this->layout_mutable() = typename static_array::layout_type(typename static_array::extensions_type{});
		assert(this->stride() != 0);
	}
	template<class... Indices>
	constexpr auto reindex(Indices... idxs) & -> static_array& {
		static_array::layout_t::reindex(idxs...);
		return *this;
	}
	template<class... Indices>
	constexpr auto reindex(Indices... idxs) && -> static_array&& {
		reindex(idxs...);
		return std::move(*this);
	}

 public:
	constexpr static_array() noexcept
	: array_alloc{}, ref(nullptr, typename static_array::extensions_type{}) {
		assert(this->stride() != 0);
		assert(this->size() == 0);
	}

#if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
	constexpr
#endif
	~static_array() /*noexcept*/ {
		assert(this->stride() != 0);
		destroy();
		assert(this->stride() != 0);
		deallocate();
	}

	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element_type const>;
	using element_move_ptr  = multi::move_ptr<typename static_array::element_ptr>;

	using reference = std::conditional_t<
		(D > 1),
		subarray<typename static_array::element_type, D - 1, typename static_array::element_ptr>,
		std::conditional_t<
			D == 1,
			typename std::iterator_traits<typename static_array::element_ptr>::reference,
			void
		>
	>;
	using const_reference = std::conditional_t<
		(D > 1),
		const_subarray<typename static_array::element_type, D - 1, typename static_array::element_ptr>,  // TODO(correaa) should be const_reference, but doesn't work witn rangev3?
		std::conditional_t<
			D == 1,
			decltype(*std::declval<typename static_array::element_const_ptr>()),
			void
		>
	>;

	using iterator       = multi::array_iterator<T, D, typename static_array::element_ptr>;
	using const_iterator = multi::array_iterator<T, D, typename static_array::element_ptr, true>;

	friend auto get_allocator(static_array const& self) -> allocator_type { return self.get_allocator(); }

	BOOST_MULTI_HD constexpr auto           data_elements() const& -> element_const_ptr { return this->base_; }
	BOOST_MULTI_HD constexpr auto           data_elements() & -> typename static_array::element_ptr { return this->base_; }
	BOOST_MULTI_HD constexpr auto           data_elements() && -> typename static_array::element_move_ptr { return std::make_move_iterator(this->base_); }

	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(static_array const& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(static_array& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(static_array&& self) { return std::move(self).data_elements(); }

	constexpr auto base() & -> typename static_array::element_ptr { return ref::base(); }
	constexpr auto base() const& -> typename static_array::element_const_ptr { return typename static_array::element_const_ptr{ref::base()}; }

	// BOOST_MULTI_FRIEND_CONSTEXPR auto base(static_array& self) -> typename static_array::element_ptr { return self.base(); }
	// BOOST_MULTI_FRIEND_CONSTEXPR auto base(static_array const& self) -> typename static_array::element_const_ptr { return self.base(); }

	constexpr auto              origin() & -> typename static_array::element_ptr { return ref::origin(); }
	constexpr auto              origin() const& -> typename static_array::element_const_ptr { return ref::origin(); }

	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(static_array& self) -> typename static_array::element_ptr { return self.origin(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(static_array const& self) -> typename static_array::element_const_ptr { return self.origin(); }

	//  private:
	//  constexpr auto rotated_aux() const {
	//      typename static_array::layout_t new_layout = this->layout();
	//      new_layout.rotate();
	//      return subarray<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	//  }

	// constexpr auto rotated() const& {return std::move(*this).rotated_aux();}
	// constexpr auto rotated()      & {return std::move(*this).rotated_aux();}
	// constexpr auto rotated()     && {return std::move(*this).rotated_aux();}

	// friend constexpr auto rotated(static_array&       self) -> decltype(auto) {return self.rotated();}
	// friend constexpr auto rotated(static_array const& self) -> decltype(auto) {return self.rotated();}

	// constexpr auto unrotated() const& -> subarray<T, D, typename static_array::element_ptr> const {
	//  typename static_array::layout_t new_layout = this->layout();
	//  new_layout.unrotate();
	//  return subarray<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	// }
	// constexpr auto unrotated()      & {
	//  typename static_array::layout_t new_layout = this->layout();
	//  new_layout.unrotate();
	//  return subarray<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	// }

	// friend constexpr auto unrotated(static_array      & self) -> decltype(auto) {return self.unrotated();}
	// friend constexpr auto unrotated(static_array const& self) -> decltype(auto) {return self.unrotated();}

	template<class TT, typename EElementPtr, class LLayout>
	auto operator=(multi::const_subarray<TT, D, EElementPtr, LLayout> const& other) -> static_array& {
		ref::operator=(other);  // TODO(correaa) : protect for self assigment
		assert(this->stride() != 0);
		return *this;
	}
	auto operator=(static_array const& other) & -> static_array& {
		if(std::addressof(other) == this) {
			return *this;
		}  // cert-oop54-cpp
		assert(other.extensions() == this->extensions());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		assert(this->stride() != 0);
		return *this;
	}

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	constexpr auto operator=(static_array&& other) noexcept -> static_array& {  // lints  (cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		assert(extensions(other) == static_array::extensions());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // there is no std::move_n algorithm
		assert(this->stride() != 0);
		return *this;
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	template<class TT, class... As>
	auto operator=(static_array<TT, D, As...> const& other) & -> static_array& {
		assert(extensions(other) == static_array::extensions());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	// constexpr explicit operator subarray<typename static_array::value_type, D, typename static_array::element_const_ptr, typename static_array::layout_t>()& {
	// return this->template static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	// }

	template<class Archive>
	void serialize(Archive& arxiv, unsigned int const version) {
		ref::serialize(arxiv, version);
	}

 private:
	void swap_(static_array& other) noexcept { operator()().swap(other()); assert(this->stride() != 0);}

 public:
	friend void swap(static_array& lhs, static_array& rhs) noexcept {
		lhs.swap_(rhs);
	}
};

template<typename T, class Alloc>
struct static_array<T, ::boost::multi::dimensionality_type{0}, Alloc>  // NOLINT(fuchsia-multiple-inheritance) : design
: protected detail::array_allocator<Alloc>
, public array_ref<T, 0, typename multi::allocator_traits<typename detail::array_allocator<Alloc>::allocator_type>::pointer> {
	static_assert(std::is_same_v<typename multi::allocator_traits<Alloc>::value_type, typename static_array::element_type>,
	              "allocator value type must match array value type");

 private:
	using array_alloc = detail::array_allocator<Alloc>;

 public:
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() && -> static_array* = delete;  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : delete to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() & -> static_array* { return this; }  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : override from base
	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() const& -> static_array const* { return this; }  // NOSONAR(cpp:S877) NOLINT(google-runtime-operator) : override from base

	using array_alloc::get_allocator;
	using allocator_type = typename static_array::allocator_type;
	using decay_type     = array<T, 0, Alloc>;

	template<class Ptr>
	void assign(Ptr data) & {
		if(data) {
			assert(this->num_elements() == 1);
			adl_copy_n(data, this->num_elements(), this->base());
		}
	}

	template<class Singleton,
	         std::enable_if_t<!std::is_base_of_v<static_array, Singleton> && !std::is_same_v<Singleton, typename static_array::element_type>, int> = 0,
	         class                                                                                                                                 = decltype(adl_copy_n(&std::declval<Singleton>(), 1, typename static_array::element_ptr{}))>
	auto operator=(Singleton const& single) -> static_array& {
		assign(&single);
		return *this;
	}

 protected:
	using alloc_traits = typename multi::allocator_traits<allocator_type>;
	using ref          = array_ref<T, 0, typename multi::allocator_traits<typename multi::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;

	auto uninitialized_value_construct() {
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type>) {
			return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());
		}
	}

	template<typename It> auto uninitialized_copy(It first) {
	#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(this->alloc(), this->data_elements(), this->num_elements());
		}
		return adl_copy                      (               first, this->num_elements(), this->data_elements());
	#else
		return adl_alloc_uninitialized_copy_n(this->alloc(), first, this->num_elements(), this->data_elements());
	#endif
	}
	template<typename It>
	auto uninitialized_move(It first) {
		return adl_alloc_uninitialized_move_n(this->alloc(), first, this->num_elements(), this->data_elements());
	}

	constexpr void destroy() {
		if constexpr(!(std::is_trivially_destructible_v<typename static_array::element_type> || multi::force_element_trivial_destruction<typename static_array::element_type>)) {
			array_alloc::destroy_n(this->data_elements(), this->num_elements());
		}
	}
	// auto destroy() {
	//  return adl_alloc_destroy_n(this->alloc(), this->data_elements(), this->num_elements());
	//  // array_alloc::destroy_n(this->data_elements(), this->num_elements());
	// }

 public:
	using typename ref::difference_type;
	using typename ref::size_type;
	using typename ref::value_type;
	constexpr explicit static_array(allocator_type const& alloc) : array_alloc{alloc} {}

	constexpr static_array(decay_type&& other, allocator_type const& alloc)  // 6b
	: array_alloc{alloc}, ref{other.base_, other.extensions()} {
		std::move(other).ref::layout_t::operator=({});
	}

	using ref::operator==;
	using ref::operator!=;

	static_array(
		typename static_array::extensions_type const& extensions,
		typename static_array::element const& elem, allocator_type const& alloc
	)
	: array_alloc{ alloc },
	  ref(
		  static_array::allocate(
			  static_cast<typename multi::allocator_traits<allocator_type>::size_type>(
				typename static_array::layout_t{ extensions }.num_elements()
			)
		  ),
		  extensions
	  ) {
		uninitialized_fill(elem);
	}

	static_array(typename static_array::element_type const& elem, allocator_type const& alloc)
	: static_array(typename static_array::extensions_type{}, elem, alloc) {}

	template<typename OtherT, typename OtherEPtr, class OtherLayout>
	explicit static_array(multi::const_subarray<OtherT, 0, OtherEPtr, OtherLayout> const& other, allocator_type const& alloc)
	: array_alloc{alloc}, ref(static_array::allocate(other.num_elements()), extensions(other)) {
		assert(other.num_elements() <= 1);
		if(other.num_elements()) {
			#if defined(__clang__) && defined(__CUDACC__)
			if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
				adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->data_elements(), this->num_elements());
			}
			adl_copy                    (                       other.base(), other.base() + other.num_elements(), this->base());
			#else
			adl_alloc_uninitialized_copy(static_array::alloc(), other.base(), other.base() + other.num_elements(), this->base());
			#endif
		}
	}

	template<class TT, class... Args>
	explicit static_array(multi::static_array<TT, 0, Args...> const& other, allocator_type const& alloc)  // TODO(correaa) : call other constructor (above)
	: array_alloc{alloc}, ref(static_array::allocate(other.num_elements()), extensions(other)) {
		#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n                    (                       other.data_elements(), other.num_elements(), this->data_elements());
		#else
		adl_alloc_uninitialized_copy_n(static_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
		#endif
	}

	template<class TT, class... Args>
	explicit static_array(multi::static_array<TT, 0, Args...> const& other)
	: static_array(other, allocator_type{}) {}

	auto uninitialized_fill(typename static_array::element_type const& elem) {
		array_alloc::uninitialized_fill_n(
			this->base_,
			static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()),
			elem
		);
	}

	template<class TT, class... Args>
	auto operator=(multi::const_subarray<TT, 0, Args...> const& other) -> static_array& {
		adl_copy_n(other.base(), 1, this->base());
		return *this;
	}

	static_array(
		typename static_array::extensions_type const& extensions,
		typename static_array::element_type const&         elem
	)  // 2
	: array_alloc{}, ref(static_array::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), nullptr), extensions) {
		uninitialized_fill(elem);
	}

	static_array() : static_array(multi::iextensions<0>{}) {}  // TODO(correaa) a noexcept will force a partially formed state for zero dimensional arrays

	explicit static_array(typename static_array::element_type const& elem)
	: static_array(multi::iextensions<0>{}, elem) {}

	template<class Singleton,
	         std::enable_if_t<!std::is_base_of_v<static_array, Singleton> && !std::is_same_v<Singleton, typename static_array::element_type>, int> = 0,
	         class                                                                                                                                 = decltype(adl_copy_n(&std::declval<Singleton>(), 1, typename static_array::element_ptr{}))>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax  // NOLINTNEXTLINE(runtime/explicit)
	static_array(Singleton const& single)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
	: ref(static_array::allocate(1), typename static_array::extensions_type{}) {
		#if defined(__clang__) && defined(__CUDACC__)
		if constexpr(! std::is_trivially_default_constructible_v<typename static_array::element_type> && ! multi::force_element_trivial_default_construction<typename static_array::element_type> ) {
			adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->data_elements(), this->num_elements());
		}
		adl_copy_n                    (                       &single, 1, this->data_elements());
		#else
		adl_alloc_uninitialized_copy_n(static_array::alloc(), &single, 1, this->data_elements());
		#endif
	}

	template<class ValueType,
		typename = std::enable_if_t<std::is_same_v<ValueType, value_type>> >  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)  // 3
	: static_array(extension * extensions(value), alloc) {
		assert(this->stride() != 0);
		using std::fill;
		fill(this->begin(), this->end(), value);
	}
	template<class ValueType,
		typename = std::enable_if_t<std::is_same_v<ValueType, value_type> >>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value)  // 3  // TODO(correaa) : call other constructor (above)
	: static_array(extension * extensions(value)) {
		assert(this->stride() != 0);
		using std::fill;
		fill(this->begin(), this->end(), value);
	}

	explicit static_array(typename static_array::extensions_type const& extensions, allocator_type const& alloc)  // 3
	: array_alloc{alloc}, ref(static_array::allocate(typename static_array::layout_t{extensions}.num_elements()), extensions) {
		assert(this->stride() != 0);
		uninitialized_value_construct();
	}
	explicit static_array(typename static_array::extensions_type const& extensions)  // 3
	: static_array(extensions, allocator_type{}) {
		assert(this->stride() != 0);
	}

	static_array(static_array const& other, allocator_type const& alloc)  // 5b
	: array_alloc{alloc}, ref(static_array::allocate(other.num_elements()), extensions(other)) {
		assert(this->stride() != 0);
		uninitialized_copy_(other.data_elements());
	}

	static_array(static_array const& other)  // 5b
	: array_alloc{other.get_allocator()}, ref{static_array::allocate(other.num_elements(), other.data_elements()), {}} {
		assert(this->stride() != 0);
		uninitialized_copy(other.data_elements());
	}

	static_array(static_array&& other) noexcept  // it is private because it is a valid operation for derived classes //5b
	: array_alloc{other.get_allocator()}, ref{static_array::allocate(static_cast<typename multi::allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), other.extensions()} {
		assert(this->stride() != 0);
		uninitialized_move(std::move(other).data_elements());
	}
	//  template<class It> static auto distance(It a, It b) {using std::distance; return distance(a, b);}

 protected:
	void deallocate() {  // TODO(correaa) : move this to detail::array_allocator
		if(this->num_elements()) {
			multi::allocator_traits<allocator_type>::deallocate(this->alloc(), this->base_, static_cast<typename multi::allocator_traits<allocator_type>::size_type>(this->num_elements()));
		}
	}
	void clear() noexcept {
		this->destroy();
		deallocate();
		layout_t<0>::operator=({});
	}

 public:
	~static_array() noexcept {
		this->destroy();
		deallocate();
	}
	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element_type const>;

	BOOST_MULTI_FRIEND_CONSTEXPR auto get_allocator(static_array const& self) -> allocator_type { return self.get_allocator(); }

	constexpr auto              base() & -> typename static_array::element_ptr { return ref::base(); }
	constexpr auto              base() const& -> typename static_array::element_const_ptr { return ref::base(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto base(static_array& self) -> typename static_array::element_ptr { return self.base(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto base(static_array const& self) -> typename static_array::element_const_ptr { return self.base(); }

	constexpr auto              origin() & -> typename static_array::element_ptr { return ref::origin(); }
	constexpr auto              origin() const& -> typename static_array::element_const_ptr { return ref::origin(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(static_array& self) -> typename static_array::element_ptr { return self.origin(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto origin(static_array const& self) -> typename static_array::element_const_ptr { return self.origin(); }

	constexpr operator typename std::iterator_traits<typename static_array::element_const_ptr>::reference() const& {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
		return *(this->base_);
	}
	constexpr operator std::add_rvalue_reference_t<typename std::iterator_traits<typename static_array::element_ptr>::reference>() && {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
		return std::move(*(this->base_));
	}
	constexpr operator typename std::iterator_traits<typename static_array::element_ptr>::reference() & {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
		return *(this->base_);
	}

	constexpr explicit operator typename static_array::element_type() const {
		return *(this->base_);
	}

	constexpr auto rotated() const& {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated() & {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated() && {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return subarray<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

	friend constexpr auto rotated(static_array& self) -> decltype(auto) { return self.rotated(); }
	friend constexpr auto rotated(static_array const& self) -> decltype(auto) { return self.rotated(); }

 private:
	constexpr auto unrotated_aux_() {
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate();
		return subarray<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}

 public:
	constexpr auto unrotated() & { return unrotated_aux_(); }
	constexpr auto unrotated() const& { return unrotated_aux_().as_const(); }

	friend constexpr auto unrotated(static_array& self) -> decltype(auto) { return self.unrotated(); }
	friend constexpr auto unrotated(static_array const& self) -> decltype(auto) { return self.unrotated(); }

	//  TODO(correaa) find a symbolic way to express rotations, A << 1, A >> 1, A <<o; A >>o; ~A; !A; ++A; A++; --A; A--; -A; +A; e<<A; A>>e; e>>A; <<A; ~A;
	//  constexpr auto operator<<(dimensionality_type d)       -> decltype(auto) {return   rotated(d);}
	//  constexpr auto operator>>(dimensionality_type d)       -> decltype(auto) {return unrotated(d);}
	//  constexpr auto operator<<(dimensionality_type d) const -> decltype(auto) {return   rotated(d);}
	//  constexpr auto operator>>(dimensionality_type d) const -> decltype(auto) {return unrotated(d);}

	constexpr auto operator=(static_array const& other) -> static_array& {
		assert(extensions(other) == static_array::extensions());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		if(this == &other) {
			return *this;
		}  // lints (cert-oop54-cpp) : handle self-assignment properly
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

 private:
	constexpr auto equal_extensions_if_(std::true_type /*true */, static_array const& other) { return this->extensions() == extensions(other); }
	constexpr auto equal_extensions_if_(std::false_type /*false*/, static_array const& /*other*/) { return true; }

 public:
	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	constexpr auto operator=(static_array&& other) noexcept -> static_array& {
		assert(equal_extensions_if_(std::integral_constant<bool, (static_array::rank_v != 0)>{}, other));  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // there is no std::move_n algorithm
		return *this;
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif

	template<class TT, class... As,
		class = std::enable_if_t<std::is_assignable<typename static_array::element_ref, TT>{}>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto operator=(static_array<TT, 0, As...> const& other) & -> static_array& {
		assert(extensions(other) == static_array::extensions());
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	constexpr explicit operator subarray<value_type, 0, typename static_array::element_const_ptr, typename static_array::layout_t>() & {
		return this->template static_array_cast<value_type, typename static_array::element_const_ptr>();
		// return static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}

	template<class Archive>
	void serialize(Archive& arxiv, unsigned int const version) {
		ref::serialize(arxiv, version);
	}
};

template<typename T, class Alloc>
struct array<T, 0, Alloc> : static_array<T, 0, Alloc> {
	// using static_ = static_array<T, 0, Alloc>;
	using static_array<T, 0, Alloc>::static_array;


	using static_array<T, 0, Alloc>::operator=;

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
	
	template<class Other, 
		std::enable_if_t<!std::is_base_of<array, std::decay_t<Other>>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto operator=(Other const& other) -> array& {
		this->assign(&other);
		return *this;
	}

	auto reextent(typename array::extensions_type const& /*empty_extensions*/) -> array& {
		return *this;
	}

	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() && -> array* = delete;  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
	// auto operator&()      & -> array      *{return this;}
	// auto operator&() const& -> array const*{return this;}
};

template<class T, ::boost::multi::dimensionality_type D, class Alloc>
struct array : static_array<T, D, Alloc> {
	~array() = default;
	using static_ = static_array<T, D, Alloc>;
	static_assert(
		std::is_same_v<
			typename multi::allocator_traits<Alloc>::value_type, std::remove_const_t<T>
			// typename array::alloc_traits::value_type, std::remove_const_t<T>
		>
		||
		std::is_same_v<
			typename multi::allocator_traits<Alloc>::value_type, void
			// typename array::alloc_traits::value_type, void
		>,
		"only exact type of array element or void (default?) is allowed as allocator value type"
	);

	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() && -> array* = delete;  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() & -> array* { return this; }  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	BOOST_MULTI_HD constexpr auto operator&() const& -> array const* { return this; }  // NOLINT(google-runtime-operator) //NOSONAR delete operator&& defined in base class to avoid taking address of temporary

	// friend auto sizes(array const& self) -> typename array::sizes_type { return self.sizes(); }

	template<class Archive, class ArTraits = multi::archive_traits<Archive>>
	void serialize(Archive& arxiv, unsigned int const version) {
		auto extensions_ = this->extensions();
		arxiv& ArTraits::make_nvp("extensions", extensions_); // don't try `using ArTraits::make_nvp`, make_nvp is a static member
		if(this->extensions() != extensions_) {
			clear();
			this->reextent(extensions_);
		}
		static_::serialize(arxiv, version);
	}

	// vvv workaround for MSVC 14.3 and ranges, TODO(correaa) good solution would be to inherit from const_subarray
	BOOST_MULTI_HD operator subarray<T, D, typename array::element_const_ptr, typename array::layout_type> const&() const {  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)
		return reinterpret_cast<subarray<T, D, typename array::element_const_ptr, typename array::layout_type> const&>(*this);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	}

	// move this to static_array
	template<
	 class Range,
	 std::enable_if_t<! has_extensions<std::decay_t<Range>>::value, int> =0,
	 class = decltype(Range{std::declval<typename array::const_iterator>(), std::declval<typename array::const_iterator>()})
	>
	constexpr explicit operator Range() const {
		// vvv Range{...} needed by Windows GCC?
		return Range{this->begin(), this->end()};  // NOLINT(fuchsia-default-arguments-calls) e.g. std::vector(it, it, alloc = {})
	}

	// move this to static_array
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	constexpr explicit operator TTN const&() const& { return this->template to_carray_<TTN>(); }
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	constexpr explicit operator TTN&() && { return this->template to_carray_<TTN>(); }
	template<class TTN, std::enable_if_t<std::is_array_v<TTN>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	constexpr explicit operator TTN&() & { return this->template to_carray_<TTN>(); }

	// NOLINTNEXTLINE(cppcoreguidelines-rvalue-reference-param-not-moved) false positive in clang-tidy 17 ?
	using static_array<T, D, Alloc>::static_array;  // MSVC wants fullname here? // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) passing c-arrays to base
	using typename static_array<T, D, Alloc>::value_type;  // MSVC wants fullname here?

#ifdef _MSC_VER
	array(typename array::extensions_type exts, typename array::allocator_type const& alloc)
	: static_array<T, D, Alloc>(exts, alloc) {assert(this->stride() != 0);}

	array(typename array::extensions_type exts)
	: static_array<T, D, Alloc>(exts) {assert(this->stride() != 0);}
#endif

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	constexpr array(std::initializer_list<typename static_array<T, D>::value_type> ilv)
	: static_{array<T, D>(ilv.begin(), ilv.end())} {
		assert(this->stride() != 0);
	}

	template<class OtherT,
		class = std::enable_if_t<std::is_constructible_v<typename static_array<T, D>::value_type, OtherT> && !std::is_convertible_v<OtherT, typename static_array<T, D>::value_type> && (D == 1)>>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	constexpr explicit array(std::initializer_list<OtherT> ilv)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) inherit explicitness of conversion from the elements
	: static_{array<T, D>(ilv.begin(), ilv.end()).element_transformed([](auto const& elem) noexcept { return static_cast<T>(elem); })} {
		assert(this->stride() != 0);
	}  // TODO(correaa) investigate why noexcept is necessary

	array()             = default;
	array(array const&) = default;

	auto reshape(typename array::extensions_type extensions) & -> array& {
		typename array::layout_t const new_layout{extensions};  // TODO(correaa) implement move-reextent in terms of reshape
		assert(new_layout.num_elements() == this->num_elements());
		this->layout_mutable() = new_layout;
		assert(this->stride() != 0);
		return *this;
	}

	auto clear() noexcept -> array& {
		static_::clear();
		assert(this->stride() != 0);
		return *this;
	}
	friend auto clear(array& self) noexcept -> array& { return self.clear(); }

	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array const& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array& self) { return self.data_elements(); }
	BOOST_MULTI_FRIEND_CONSTEXPR auto data_elements(array&& self) { return std::move(self).data_elements(); }

	// auto move() & -> subarray<typename array::element_type, D, multi::move_ptr<typename array::element_type>> {
	//  subarray<typename array::element_type, D, multi::move_ptr<typename array::element_type>>
	//      ret = multi::static_array_cast<typename array::element_type, multi::move_ptr<typename array::element_type>>(*this);

	//  layout_t<D>::operator=({});

	//  assert(this->stride() != 0);
	//  return ret;
	// }
	// friend auto move(array& self) -> subarray<typename array::element_type, D, multi::move_ptr<typename array::element_type>> {
	//  return self.move();
	// }

	friend BOOST_MULTI_HD constexpr auto move(array& self) -> decltype(auto) { return std::move(self); }
	friend BOOST_MULTI_HD constexpr auto move(array&& self) -> decltype(auto) { return std::move(self); }

	array(array&& other, typename array::allocator_type const& alloc) noexcept : static_array<T, D, Alloc>{std::move(other), alloc} {
		assert(this->stride() != 0);
	}
	array(array&& other) noexcept : array{std::move(other), other.get_allocator()} {
		assert(this->stride() != 0);
	}

	friend auto get_allocator(array const& self) -> typename array::allocator_type { return self.get_allocator(); }

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
		this->layout_mutable() = std::exchange(other.layout_mutable(), typename array::layout_type(typename array::extensions_type{}));
		assert(this->stride() != 0);
		assert(other.stride() != 0);
		return *this;
	}

	auto operator=(array const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			if(this == &other) {
				return *this;
			}  // required by cert-oop54-cpp
			if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			static_::operator=(other);
		} else {
			clear();
			if constexpr(multi::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			this->layout_mutable() = other.layout();
			array::allocate();
			array::uninitialized_copy_elements(other.data_elements());
		}
		assert(this->stride() != 0);
		return *this;
	}
#else
	auto operator=(array o) noexcept -> array& { return swap(o), *this; }
#endif

	template<typename OtherT, typename OtherEP, class OtherLayout>
	auto operator=(multi::const_subarray<OtherT, D, OtherEP, OtherLayout> const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			static_::operator=(other);  // TODO(correaa) : protect for self assigment
		} else {
			operator=(array{other});
		}
		return *this;
	}

	template<class TT, class AAlloc>  // , std::enable_if_t<not std::is_base_of_v<array, multi::array<TT, D, AAlloc>> , int> =0>
	auto operator=(multi::array<TT, D, AAlloc> const& other) -> array& {  // TODO(correaa) : check that LHS is not read-only?
		if(array::extensions() == other.extensions()) {
			//  this->operator()() = other;
			static_::operator=(other);
		} else if(this->num_elements() == other.extensions().num_elements()) {
			reshape(other.extensions());
			static_::operator=(other);
			//  this->operator()() = other;
		} else {
			operator=(static_cast<array>(other));
		}
		assert(this->stride() != 0);
		return *this;
	}

	template<
		class Range,
		class                                                                 = decltype(std::declval<static_&>().operator=(std::declval<Range&&>())),
		std::enable_if_t<!has_data_elements<std::decay_t<Range>>::value, int> = 0,
		std::enable_if_t<has_extensions<std::decay_t<Range>>::value, int> = 0,
		std::enable_if_t<!std::is_base_of<array, std::decay_t<Range>>{}, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
	auto operator=(Range&& other) -> array& {  // TODO(correaa) : check that LHS is not read-only?
		if(array::extensions() == other.extensions()) {
			this->operator()() = std::forward<Range>(other);
			//  static_::operator=(other);
		} else if(this->num_elements() == other.extensions().num_elements()) {
			reshape(other.extensions());
			//  static_::operator=(other);
			this->operator()() = std::forward<Range>(other);
		} else {
			operator=(static_cast<array>(std::forward<Range>(other)));
		}
		return *this;
	}

	template<
		class Range,
		class                                                                 = decltype(std::declval<static_&>().operator=(std::declval<Range&&>())),
		std::enable_if_t<!std::is_base_of<array, std::decay_t<Range>>{}, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
	auto from(Range&& other) -> array& {  // TODO(correaa) : check that LHS is not read-only?
		if(array::extensions() == other.extensions()) {
			this->operator()() = other;
			//  static_::operator=(other);
		} else if(this->num_elements() == other.extensions().num_elements()) {
			reshape(other.extensions());
			this->operator()() = other;
			//  static_::operator=(other);
		} else {
			operator=(static_cast<array>(std::forward<Range>(other)));
		}
		return *this;
	}

	friend void swap(array& self, array& other) noexcept(true /*noexcept(self.swap(other))*/) { self.swap(other); }

	void assign(typename array::extensions_type extensions, typename array::element_type const& elem) {
		if(array::extensions() == extensions) {
			adl_fill_n(this->base_, this->num_elements(), elem);
		} else {
			this->clear();
			(*this).array::layout_t::operator=(layout_t<D>{extensions});
			this->base_ = this->static_::array_alloc::allocate(this->num_elements(), nullptr);
			adl_alloc_uninitialized_fill_n(this->alloc(), this->base_, this->num_elements(), elem);
		}
	}

	template<class It>
	auto assign(It first, It last) -> array& {
		using std::all_of;
		using std::next;
		if(adl_distance(first, last) == this->size()) {
			static_::ref::assign(first);
		} else {
			this->operator=(array(first, last));
		}
		return *this;
	}
	void assign(std::initializer_list<value_type> values) { assign(values.begin(), values.end()); }

	template<class Range> auto assign(Range&& other) & -> decltype(assign(adl_begin(std::forward<Range>(other)), adl_end(std::forward<Range>(other)))) {  // TODO(correaa) use forward
		return assign(adl_begin(std::forward<Range>(other)), adl_end(std::forward<Range>(other)));
	}

	auto operator=(std::initializer_list<value_type> values) -> array& {
		assign(values.begin(), values.end());
		return *this;
	}

	// template<class... TTs>
	// [[deprecated("use extensions for reextents, not tuples")]]
	// auto reextent(std::tuple<TTs...> const& other) -> array& {
	//  return reextent(
	//    std::apply([](auto const&... extensions) {return typename array::extensions_type(extensions...);}, other)
	//  );  // paren is important here ext_type(...) for allow narrowing casts ^^^
	// }

	auto reextent(typename array::extensions_type const& extensions) && -> array&& {
		if(extensions == this->extensions()) {
			return std::move(*this);
		}
		this->destroy();
		this->deallocate();
		this->layout_mutable() = typename array::layout_t{extensions};
		this->base_            = this->static_::array_alloc::allocate(
            static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(
                typename array::layout_t{extensions}.num_elements()
            ),
            this->data_elements()  // used as hint
        );
		if constexpr(!(std::is_trivially_default_constructible_v<typename array::element_type> || multi::force_element_trivial_default_construction<typename array::element_type>)) {
			adl_alloc_uninitialized_value_construct_n(this->alloc(), this->base_, this->num_elements());
		}
		return std::move(*this);
	}

	auto reextent(typename array::extensions_type const& extensions) & -> array& {
		if(extensions == this->extensions()) {
			return *this;
		}
		auto&& tmp = typename array::ref(
			this->static_::array_alloc::allocate(
				static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(
					typename array::layout_t{extensions}.num_elements()
				),
				this->data_elements()  // used as hint
			),
			extensions
		);
		if constexpr(!(std::is_trivially_default_constructible_v<typename array::element_type> || multi::force_element_trivial_default_construction<typename array::element_type>)) {
			adl_alloc_uninitialized_value_construct_n(this->alloc(), tmp.data_elements(), tmp.num_elements());
		}
		auto const is = intersection(this->extensions(), extensions);
		tmp.apply(is) = this->apply(is);  // TODO(correaa) : use (and implement) `.move();`
		this->destroy();
		this->deallocate();
		this->base_            = tmp.base();
		this->layout_mutable() = tmp.layout();
		return *this;
	}

	[[nodiscard]] constexpr auto operator+() const& { return array{*this}; }
	[[nodiscard]] constexpr auto operator+() && { return array{*this}; }

	auto reextent(typename array::extensions_type const& exs, typename array::element_type const& elem) & -> array& {
		if(exs == this->extensions()) {
			return *this;
		}

		// array tmp(x, e, this->get_allocator());  // TODO(correaa) opportunity missed to use hint allocation
		// auto const is = intersection(this->extensions(), x);
		// tmp.apply(is) = this->apply(is);
		// swap(tmp);

		// implementation with hint
		auto&& tmp = typename array::ref(
			this->static_::array_alloc::allocate(
				static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(typename array::layout_t{exs}.num_elements()),
				this->data_elements()  // use as hint
			),
			exs
		);
		this->uninitialized_fill_n(tmp.data_elements(), static_cast<typename multi::allocator_traits<typename array::allocator_type>::size_type>(tmp.num_elements()), elem);
		auto const is = intersection(this->extensions(), exs);
		tmp.apply(is) = this->apply(is);
		this->destroy();
		this->deallocate();
		this->base_            = tmp.base();  // TODO(correaa) : use (and implement) `.move();`
		this->layout_mutable() = tmp.layout();
		//  (*this).array::layout_t::operator=(tmp.layout());

		return *this;
	}
	template<class... Indices> constexpr auto reindex(Indices... idxs) && -> array&& {
		this->layout_mutable().reindex(idxs...);
		return std::move(*this);
	}
	template<class... Indices> constexpr auto reindex(Indices... idxs) & -> array& {
		this->layout_mutable().reindex(idxs...);
		return *this;
	}

	// ~array() {
	//  assert(this->stride() != 0);
	// }
};

#if defined(__cpp_deduction_guides)

#define BOOST_MULTI_IL std::initializer_list  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing TODO(correaa) remove

// vvv MSVC 14.3 in c++17 mode needs paranthesis in dimensionality_type(d)
template<class T> static_array(BOOST_MULTI_IL<T>) -> static_array<T, static_cast<dimensionality_type>(1U), std::allocator<T>>;  // MSVC needs the allocator argument error C2955: 'boost::multi::static_array': use of class template requires template argument list
template<class T> static_array(BOOST_MULTI_IL<BOOST_MULTI_IL<T>>) -> static_array<T, static_cast<dimensionality_type>(2U), std::allocator<T>>;
template<class T> static_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>) -> static_array<T, static_cast<dimensionality_type>(3U), std::allocator<T>>;
template<class T> static_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>) -> static_array<T, static_cast<dimensionality_type>(4U), std::allocator<T>>;
template<class T> static_array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>>) -> static_array<T, static_cast<dimensionality_type>(5U), std::allocator<T>>;

// TODO(correaa) add zero dimensional case?
template<class T> array(BOOST_MULTI_IL<T>) -> array<T, static_cast<dimensionality_type>(1U)>;
template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<T>>) -> array<T, static_cast<dimensionality_type>(2U)>;
template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>) -> array<T, static_cast<dimensionality_type>(3U)>;
template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>) -> array<T, static_cast<dimensionality_type>(4U)>;
template<class T> array(BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<BOOST_MULTI_IL<T>>>>>) -> array<T, static_cast<dimensionality_type>(5U)>;

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

template<class MatrixRef, class DT = typename MatrixRef::decay_type, class T = typename DT::element_type, dimensionality_type D = DT::rank_v, class Alloc = typename DT::allocator_type>
array(MatrixRef) -> array<T, D, Alloc>;

template<typename T, dimensionality_type D, typename P> array(subarray<T, D, P>) -> array<T, D>;

template<
	class Range, 
	std::enable_if_t<!has_extensions<Range>::value, int> = 0,
	typename V = decltype(*::std::begin(std::declval<Range const&>()))
	// typename V = typename std::iterator_traits<decltype(::std::begin(std::declval<Range const&>()))>::value_type
>
array(Range) -> array<V, 1>;

template<class Reference>
auto operator+(Reference&& ref)
->decltype(array(std::forward<Reference>(ref))) {
	return array(std::forward<Reference>(ref)); }

#endif  // ends defined(__cpp_deduction_guides)

template<class T, std::size_t N>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
auto decay(const T (&arr)[N]) noexcept -> multi::array<std::remove_all_extents_t<T[N]>, std::rank_v<T[N]>> {
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	return multi::array_cref<std::remove_all_extents_t<T[N]>, std::rank_v<T[N]>>(data_elements(arr), extensions(arr));
}

template<class T, std::size_t N>
struct array_traits<T[N], void, void> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using reference  = T&;
	using element    = std::remove_all_extents_t<T[N]>;  // NOSONAR(cpp:S5945) NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using decay_type = multi::array<T, 1>;
};

}  // end namespace boost::multi

namespace boost::multi::pmr {

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
template<class T, boost::multi::dimensionality_type D>
using array = boost::multi::array<T, D, std::pmr::polymorphic_allocator<T>>;
#else
template<class T, boost::multi::dimensionality_type D>
struct [[deprecated("no PMR allocator")]] array;  // your version of C++ doesn't provide polymorphic_allocators
#endif

}  // end namespace boost::multi::pmr

// common_reference for compatibility with ranges
#if defined(__cpp_lib_common_reference) || defined(__cpp_lib_ranges)
// TODO(correaa) achieve this by normal inheritance
// NOLINTBEGIN(cert-dcl58-cpp)
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<typename ::boost::multi::array<T, D, A...>::basic_const_array     &&,          ::boost::multi::array<T, D, A...>                         &> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array     &&; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<typename ::boost::multi::array<T, D, A...>::basic_const_array     &&,          ::boost::multi::array<T, D, A...>                    const&> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array     &&; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<         ::boost::multi::array<T, D, A...>                         &, typename ::boost::multi::array<T, D, A...>::basic_const_array     &&> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array     &&; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<         ::boost::multi::array<T, D, A...>                    const&, typename ::boost::multi::array<T, D, A...>::basic_const_array     &&> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array     &&; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<typename ::boost::multi::array<T, D, A...>::basic_const_array       ,          ::boost::multi::array<T, D, A...>                         &> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array       ; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<         ::boost::multi::array<T, D, A...>                    const&, typename ::boost::multi::array<T, D, A...>::basic_const_array const&> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array const&; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<typename ::boost::multi::array<T, D, A...>::basic_const_array const&,          ::boost::multi::array<T, D, A...>                    const&> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array const&; };

// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<typename ::boost::multi::array<T, D, A...>::basic_const_array      &,          ::boost::multi::array<T, D, A...>                         &> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array      &; };
// template<class T, ::boost::multi::dimensionality_type D, class... A> struct std::common_reference<         ::boost::multi::array<T, D, A...>                         &, typename ::boost::multi::array<T, D, A...>::basic_const_array      &> { using type = typename ::boost::multi::array<T, D, A...>::basic_const_array      &; };
// NOLINTEND(cert-dcl58-cpp)
#endif

namespace boost::serialization {

template<typename T> struct version;  // in case serialization was not included before

template<typename T, boost::multi::dimensionality_type D, class A>
struct version<boost::multi::array<T, D, A>> {
	using type = std::integral_constant<int, BOOST_MULTI_SERIALIZATION_ARRAY_VERSION>;  // TODO(correaa) use constexpr variable here, not macro
	enum /*class value_t*/ { value = type::value };  // NOSONAR(cpp:S3642)  // https://community.sonarsource.com/t/suppress-issue-in-c-source-file/43154/24
};

}  // end namespace boost::serialization

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_ARRAY_HPP_
