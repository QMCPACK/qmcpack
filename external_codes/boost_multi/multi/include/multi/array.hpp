// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MULTI_ARRAY_HPP_
#define BOOST_MULTI_ARRAY_HPP_

#include "./array_ref.hpp"
#include "./config/NO_UNIQUE_ADDRESS.hpp"

#include "./detail/adl.hpp"
#include "./detail/memory.hpp"
#include "./detail/type_traits.hpp"

#include "./memory/allocator.hpp"

#include<algorithm>  // for copy
#include<memory>     // for allocator_traits
#include<tuple>      // needed by a deprecated function
#include<utility>    // for move

namespace boost::multi {

template<class Allocator>
struct array_allocator {
	using allocator_type = Allocator;

 private:
	MULTI_NO_UNIQUE_ADDRESS allocator_type alloc_;

	using allocator_traits = typename std::allocator_traits<allocator_type>;
	using size_type_ = typename allocator_traits::size_type;
	using pointer_   = typename allocator_traits::pointer;

 protected:
	auto alloc()      & -> allocator_type      & {return alloc_;}
	auto alloc() const& -> allocator_type const& {return alloc_;}

	array_allocator() = default;
	explicit array_allocator(allocator_type const& a) : alloc_{a} {}

	auto allocate(size_type_ n) -> pointer_ {
		return n?allocator_traits::allocate(alloc_, n):pointer_{nullptr};
	}
	auto allocate(size_type_ n, typename allocator_traits::const_void_pointer hint) -> pointer_ {
		return n?allocator_traits::allocate(alloc_, n, hint):pointer_{nullptr};
	}

	auto uninitialized_fill_n(pointer_ base, size_type_ n, typename allocator_traits::value_type e) {
		return adl_alloc_uninitialized_fill_n(alloc_, base, n, e);
	}
	template<typename It>
	auto uninitialized_copy_n(It first, size_type n, pointer_ data) {
		return adl_alloc_uninitialized_copy_n(alloc_, first, n, data);
	}
	template<typename It>
	auto destroy_n(It first, size_type n) {return adl_alloc_destroy_n(this->alloc(), first, n);}

 public:
	constexpr auto get_allocator() const -> allocator_type {return alloc_;}
};

template<class T, dimensionality_type D, class Alloc = std::allocator<T>>
struct static_array  // NOLINT(fuchsia-multiple-inheritance) : multiple inheritance used for composition
: protected array_allocator<Alloc>
, public array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>
, boost::multi::random_iterable<static_array<T, D, Alloc>> {
 protected:
	using array_alloc = array_allocator<Alloc>;

 public:
	static_assert( std::is_same<typename std::allocator_traits<Alloc>::value_type, typename static_array::element>{},
		"allocator value type must match array value type");
//  static_assert( std::is_same<typename std::allocator_traits<Alloc>::pointer, typename static_array::element_ptr>{},
//  	"allocator pointer type must match array pointer type");
	using array_alloc::get_allocator;
	using typename array_allocator<Alloc>::allocator_type;
//  using allocator_type = typename static_array::allocator_type;
	using decay_type = array<T, D, Alloc>;
	using layout_type = typename array_ref<T, D, typename std::allocator_traits<Alloc>::pointer>::layout_type;

 protected:
	using alloc_traits = typename std::allocator_traits<typename static_array::allocator_type>;
	using ref = array_ref<T, D, typename std::allocator_traits<typename std::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;

	auto uninitialized_value_construct() {
		return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}

	auto uninitialized_default_construct_if(std::true_type /*true*/ ) {}
	auto uninitialized_default_construct_if(std::false_type/*false*/) {
		return adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}

	auto uninitialized_default_construct() {
		return uninitialized_default_construct_if(multi::is_trivially_default_constructible<typename static_array::element_type>{});
	}

	template<typename It> auto uninitialized_copy_elements(It first) {
		return array_alloc::uninitialized_copy_n(first, this->num_elements(), this->data_elements());
	}

	void destroy_if_not(std::true_type /*true */) {}
	void destroy_if_not(std::false_type/*false*/) {array_alloc::destroy_n(this->data_elements(), this->num_elements());}
	void destroy() {destroy_if_not(std::is_trivially_destructible<typename static_array::element>{});}

	void allocate() {this->base_ = array_alloc::allocate(static_array::num_elements());}

 public:
	using value_type = typename std::conditional<
		(D > 1),  // this parenthesis is needed
		array<typename static_array::element, D-1, allocator_type>,
		typename static_array::element
	>::type;

	using typename ref::size_type;
	using typename ref::difference_type;
	explicit static_array(typename static_array::allocator_type const& a) : array_alloc{a} {}

 protected:
	static_array(static_array&& other, typename static_array::allocator_type const& a) noexcept  // 6b
	: array_alloc{a}  // TODO(correaa) : handle allocation propagation here
	, ref{other.base_, other.extensions()} {
		other.ref::layout_t::operator=({});
		other.base_ = nullptr;
	}

	static_array(static_array&& other) noexcept
	: static_array(std::move(other), typename static_array::allocator_type{}) {}  // 6b

 public:
	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // decltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>
	// analogous to std::vector::vector (5) https://en.cppreference.com/w/cpp/container/vector/vector
	static_array(It first, It last, typename static_array::allocator_type const& a)
	: array_alloc{a}
	, ref {
		array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(layout_type {index_extension {adl_distance(first, last)}*multi::extensions(*first)}.num_elements())),
		index_extension {adl_distance(first, last)}*multi::extensions(*first)
	} {
		adl_alloc_uninitialized_copy(static_array::alloc(), first, last, ref::begin());
	}

	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // decltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>
	// analogous to std::vector::vector (5) https://en.cppreference.com/w/cpp/container/vector/vector
	static_array(It first, It last) : static_array(first, last, allocator_type{}) {}

	template<
		class Range, class = std::enable_if_t<not std::is_base_of<static_array, std::decay_t<Range>>{}>,
		class = decltype(/*static_array*/(std::declval<Range&&>().begin() - std::declval<Range&&>().end())),  // instantiation of static_array here gives a compiler error in 11.0, partially defined type?
		class = std::enable_if_t<not is_basic_array<Range&&>{}>
	>
	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions // NOLINTNEXTLINE(runtime/explicit)
	static_array(Range&& rng)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array{std::forward<Range>(rng).begin(), std::forward<Range>(rng).end()} {}

	template<class TT>
	auto uninitialized_fill_elements(TT const& value) {
		return array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), value);
	}

	// vvv TODO(correaa) : check if really necessary
	template<class TT, class... As>
	static_array(array_ref<TT, D, As...> const& other, typename static_array::allocator_type const& a)
	: array_alloc{a}
	, ref{
		array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(other.num_elements())),
		other.extensions()
	} {
		adl_alloc_uninitialized_copy_n(static_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
	}

	template<class TT, class... As>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented  // NOLINTNEXTLINE(runtime/explicit)
	static_array(array_ref<TT, D, As...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array(other, allocator_type{}) {}
	// ^^^ TODO(correaa) : check if really necessary

	static_array(typename static_array::extensions_type x, typename static_array::element const& e, typename static_array::allocator_type const& a)  // 2
	: array_alloc{a}
	, ref{array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{x}.num_elements())), x} {
		array_alloc::uninitialized_fill_n(this->data_elements(), static_cast<typename std::allocator_traits<allocator_type>::size_type>(this->num_elements()), e);
	}

	template<class Element, std::enable_if_t<std::is_convertible<Element, typename static_array::element>{} and (D == 0), int> = 0>
	explicit static_array(Element const& e, allocator_type const& a)
	: static_array(typename static_array::extensions_type{}, e, a) {}

	static_array(typename static_array::extensions_type x, typename static_array::element const& e)  // 2
	: array_alloc{}
	, ref{array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{x}.num_elements())), x} {
		array_alloc::uninitialized_fill_n(this->base(), static_cast<typename std::allocator_traits<allocator_type>::size_type>(this->num_elements()), e);
	}

	explicit static_array(typename static_array::extensions_type x, typename std::allocator_traits<Alloc>::const_void_pointer hint)
	: array_alloc{}
	, ref{array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{x}.num_elements()), hint), x} {}

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v, typename static_array::allocator_type const& a)  // 3
	= delete;

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v)  // 3
	= delete;

// analgous to std::vector::vector ((4)) https://en.cppreference.com/w/cpp/container/vector/vector
	explicit static_array(typename static_array::extensions_type x, typename static_array::allocator_type const& a)
	: array_alloc{a}
	, ref{array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{x}.num_elements())), x} {
		uninitialized_default_construct();
	}

	explicit static_array(typename static_array::extensions_type x)
	: static_array(x, typename static_array::allocator_type{}) {}

	template<class TT, class... Args,
		class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::basic_array<TT, D, Args...>::element>{}>,
		class = decltype(adl_copy(std::declval<multi::basic_array<TT, D, Args...> const&>().begin(), std::declval<multi::basic_array<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	static_array(multi::basic_array<TT, D, Args...> const& o, typename static_array::allocator_type const& a)
	: static_array(o.extensions(), a) {
		adl_uninitialized_copy(o.begin(), o.end(), this->begin());  // TODO(correaa): call this conditionally on T properties
	}

	template<class TT, class... Args,
	//  class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::basic_array<TT, D, Args...>::element>{}>,
		class = decltype(adl_copy(std::declval<multi::basic_array<TT, D, Args...> const&>().begin(), std::declval<multi::basic_array<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented  // NOLINTNEXTLINE(runtime/explicit)
	static_array(multi::basic_array<TT, D, Args...> const& o)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array(o, typename static_array::allocator_type{}) {}

	template<class TT, class... Args>
	explicit static_array(array_ref<TT, D, Args...>&& o)
	: array_alloc{}
	, ref{array_alloc::allocate(o.num_elements()), o.extensions()} {
		static_array::uninitialized_copy_elements(std::move(o).data_elements());
	}

	static_array(static_array const& o)                                    // 5b
	: array_alloc{std::allocator_traits<Alloc>::select_on_container_copy_construction(o.alloc())}
	, ref{array_alloc::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(o.num_elements()), o.data_elements()), extensions(o)} {
		uninitialized_copy_elements(o.data_elements());
	}

// TODO(correaa) static_array(static_array&& o)                           // 5b'
//  : array_alloc{o.get_allocator()}, ref{array_alloc::allocate(num_elements(o)), extensions(o)} {
//      array_alloc::uninitialized_move_elements(data_elements(o));
//  }

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	constexpr static_array(std::initializer_list<typename static_array<T, D>::value_type> mil)
	: static_array(static_array<T, D>(mil.begin(), mil.end())) {}

	static_array(
		std::initializer_list<typename static_array<T, D>::value_type> mil,
		typename static_array::allocator_type const& a
	) : static_array{static_array<T, D>(mil.begin(), mil.end()), a} {}

	template<class TT, std::size_t N>
	constexpr explicit static_array(TT(&array)[N]) : static_array(std::begin(array), std::end(array)) {}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backward compatibility
	template<class It> static auto distance(It a, It b) {
		using std::distance;
		return distance(a, b);
	}

 protected:
	void deallocate() {
		if(this->num_elements()) {
			alloc_traits::deallocate(this->alloc(), this->base_, static_cast<typename alloc_traits::size_type>(this->num_elements()));
		}
	}
	void clear() noexcept {
		this->destroy();
		deallocate();
		layout_t<D>::operator=({});
	}
	template<class... Ts>
	constexpr auto reindex(Ts... a)& -> static_array& {
		static_array::layout_t::reindex(a...);
		return *this;
	}
	template<class... Ts>
	constexpr auto reindex(Ts... a)&& -> static_array&& {
	//  return std::move(reindex(a...));
		reindex(a...);
		return std::move(*this);
	}

 public:
	static_array() = default;
	~static_array() noexcept {destroy(); deallocate();}

	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element const>;
	using element_move_ptr  = std::move_iterator<typename static_array::element_ptr>;

	using reference = typename std::conditional<
		(D > 1),
		basic_array<typename static_array::element, D-1, typename static_array::element_ptr>,
		typename std::conditional<
			D == 1,
			typename std::iterator_traits<typename static_array::element_ptr>::reference,
			void
		>::type
	>::type;
	using const_reference = typename std::conditional<
		(D > 1),
		basic_array<typename static_array::element, D-1, typename static_array::element_const_ptr>,  // TODO(correaa) should be const_reference, but doesn't work witn rangev3?
		typename std::conditional<
			D == 1,
			decltype(*std::declval<typename static_array::element_const_ptr>()),
			void
		>::type
	>::type;

	using       iterator = multi::array_iterator<T, D, typename static_array::element_ptr      >;
	using const_iterator = multi::array_iterator<T, D, typename static_array::element_const_ptr>;

	friend
	#if not defined(__NVCC__) and not defined(__INTEL_COMPILER)
	constexpr
	#endif
	auto get_allocator(static_array const& s) -> typename static_array::allocator_type {return s.get_allocator();}

	       constexpr auto data_elements()            const& ->                        element_const_ptr {return this->base_;}
	       constexpr auto data_elements()                 & -> typename static_array::element_ptr       {return this->base_;}
	       constexpr auto data_elements()                && -> typename static_array::element_move_ptr  {return std::make_move_iterator(this->base_);}
	friend constexpr auto data_elements(static_array const& s) {return           s .data_elements();}
	friend constexpr auto data_elements(static_array      & s) {return           s .data_elements();}
	friend constexpr auto data_elements(static_array     && s) {return std::move(s).data_elements();}

	       constexpr auto base()                 &    -> typename static_array::element_ptr       {return ref::base();}
	       constexpr auto base()            const&    -> typename static_array::element_const_ptr {return typename static_array::element_const_ptr{ref::base()};}
	friend constexpr auto base(static_array      & s) -> typename static_array::element_ptr       {return s.base();}
	friend constexpr auto base(static_array const& s) -> typename static_array::element_const_ptr {return s.base();}

	       constexpr auto origin()                 &    -> typename static_array::element_ptr       {return ref::origin();}
	       constexpr auto origin()            const&    -> typename static_array::element_const_ptr {return ref::origin();}
	friend constexpr auto origin(static_array      & s) -> typename static_array::element_ptr       {return    s.origin();}
	friend constexpr auto origin(static_array const& s) -> typename static_array::element_const_ptr {return    s.origin();}

 private:
	constexpr auto rotated_aux() const {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}

 public:
//  constexpr auto rotated(dimensionality_type d) const& {
//  	typename static_array::layout_t new_layout = this->layout();
//  	new_layout.rotate(d);
//  	return basic_array<T, D, typename static_array::element_const_ptr>{new_layout, this->base_};
//  }
//  constexpr auto rotated(dimensionality_type d)& {
//  	typename static_array::layout_t new_layout = this->layout();
//  	new_layout.rotate(d);
//  	return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
//  }
//  constexpr auto rotated(dimensionality_type d)&& {
//  	typename static_array::layout_t new_layout = this->layout();
//  	new_layout.rotate(d);
//  	return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
//  }

	constexpr auto rotated() const& {return std::move(*this).rotated_aux();}
	constexpr auto rotated()      & {return std::move(*this).rotated_aux();}
	constexpr auto rotated()     && {return std::move(*this).rotated_aux();}

	friend constexpr auto rotated(static_array&       s) -> decltype(auto) {return s.rotated();}
	friend constexpr auto rotated(static_array const& s) -> decltype(auto) {return s.rotated();}

	constexpr auto unrotated() const& {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return basic_array<T, D, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto unrotated()      & {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}

//  constexpr auto unrotated() const& {return unrotated(1);}
//  constexpr auto unrotated()      & {return unrotated(1);}

	friend constexpr auto unrotated(static_array      & self) -> decltype(auto) {return self.unrotated();}
	friend constexpr auto unrotated(static_array const& self) -> decltype(auto) {return self.unrotated();}

//  constexpr auto operator<<(dimensionality_type d)       -> decltype(auto) {return   rotated(d);}
//  constexpr auto operator>>(dimensionality_type d)       -> decltype(auto) {return unrotated(d);}
//  constexpr auto operator<<(dimensionality_type d) const -> decltype(auto) {return   rotated(d);}
//  constexpr auto operator>>(dimensionality_type d) const -> decltype(auto) {return unrotated(d);}

	template<class TT, class... Args>
	auto operator=(multi::basic_array<TT, D, Args...> const& other) -> static_array& {
		ref::operator=(other);  // TODO(correaa) : protect for self assigment
		return *this;
	}
	auto operator=(static_array const& other) & -> static_array& {
		assert( extensions(other) == static_array::extensions() );
		if(this == &other) {return *this;}  // lints (cert-oop54-cpp) : handle self-assignment properly
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}
	constexpr auto operator=(static_array&& other) noexcept -> static_array& {  // lints  (cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		assert( extensions(other) == static_array::extensions() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // there is no std::move_n algorithm
		return *this;
	}
	template<class TT, class... As>
	auto operator=(static_array<TT, D, As...> const& other) & -> static_array& {
		assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}
	constexpr explicit operator basic_array<typename static_array::value_type, D, typename static_array::element_const_ptr, typename static_array::layout_t>()& {
		return this->template static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		ref::serialize(ar, version);
	}
};

template<class T, class Alloc>
struct static_array<T, 0, Alloc>  // NOLINT(fuchsia-multiple-inheritance) : design
: protected array_allocator<Alloc>
, public array_ref<T, 0, typename std::allocator_traits<typename array_allocator<Alloc>::allocator_type>::pointer> {
	static_assert( std::is_same<typename std::allocator_traits<Alloc>::value_type, typename static_array::element>{},
		"allocator value type must match array value type");

 private:
	using array_alloc = array_allocator<Alloc>;

 public:
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&()     && -> static_array      * = delete;       // NOLINT(google-runtime-operator) : delete to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&()      & -> static_array      * {return this;}  // NOLINT(google-runtime-operator) : override from base
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&() const& -> static_array const* {return this;}  // NOLINT(google-runtime-operator) : override from base

	using array_alloc::get_allocator;
	using allocator_type = typename static_array::allocator_type;

 protected:
	using alloc_traits = typename std::allocator_traits<allocator_type>;
	using ref = array_ref<T, 0, typename std::allocator_traits<typename std::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;

	auto uninitialized_value_construct_if_not(std::true_type /*true */) {}
	auto uninitialized_value_construct_if_not(std::false_type/*false*/) {
		return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}
	auto uninitialized_value_construct() {
		uninitialized_value_construct_if_not(std::is_trivially_default_constructible<typename static_array::element>{});
	}

	template<typename It> auto uninitialized_copy(It first) {return adl_alloc_uninitialized_copy_n(this->alloc(), first, this->num_elements(), this->data_elements());}
	template<typename It>
	auto uninitialized_move(It first) {
		return adl_alloc_uninitialized_move_n(this->alloc(), first, this->num_elements(), this->data_elements());
	}
	void destroy() {array_alloc::destroy_n(this->data_elements(), this->num_elements());}

 public:
	using typename ref::value_type;
	using typename ref::size_type;
	using typename ref::difference_type;
	constexpr explicit static_array(allocator_type const& a) : array_alloc{a} {}

 protected:
	constexpr static_array(static_array&& other, allocator_type const& a)  // 6b
	: array_alloc{a}
	, ref{other.base_, other.extensions()} {
		other.ref::layout_t::operator=({});
	}

 public:
	using ref::operator==;
	using ref::operator!=;

	static_array(typename static_array::extensions_type x, typename static_array::element const& e, allocator_type const& a)  // 2
	: array_alloc{a}
	, ref(static_array::allocate(typename static_array::layout_t{x}.num_elements()), x) {
		uninitialized_fill(e);
	}

	static_array(typename static_array::element_type const& e, allocator_type const& a)
	: static_array(typename static_array::extensions_type{}, e, a) {}

	auto uninitialized_fill(typename static_array::element const& e) {
		array_alloc::uninitialized_fill_n(this->base_, static_cast<typename std::allocator_traits<allocator_type>::size_type>(this->num_elements()), e);
	}

	static_array(typename static_array::extensions_type const& x, typename static_array::element const& e)  // 2
	: array_alloc{}, ref(static_array::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{x}.num_elements())), x) {
		uninitialized_fill(e);
	}

	static_array() : static_array(multi::iextensions<0>{}) {}

	explicit static_array(typename static_array::element const& e)          // 2
	: static_array(multi::iextensions<0>{}, e) {}

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v, allocator_type const& a)  // 3
	: static_array(e*extensions(v), a) {
		using std::fill; fill(this->begin(), this->end(), v);
	}
	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v)  // 3  // TODO(correaa) : call other constructor (above)
	: static_array(e*extensions(v)) {
		using std::fill; fill(this->begin(), this->end(), v);
	}

	explicit static_array(typename static_array::extensions_type const& x, allocator_type const& a)  // 3
	: array_alloc{a}, ref{static_array::allocate(typename static_array::layout_t{x}.num_elements()), x} {
		uninitialized_value_construct();
	}
	explicit static_array(typename static_array::extensions_type const& x)  // 3
	: static_array(x, allocator_type{}) {}

	template<class TT, class... Args>
	explicit static_array(multi::basic_array<TT, 0, Args...> const& other, allocator_type const& a)
	: array_alloc{a}
	, ref(static_array::allocate(other.num_elements()), extensions(other)) {
		using std::copy; copy(other.begin(), other.end(), this->begin());
	}
	template<class TT, class... Args>
	explicit static_array(multi::basic_array<TT, 0, Args...> const& other)  // TODO(correaa) : call other constructor (above)
	: array_alloc{}, ref(static_array::allocate(other.num_elements())
	, extensions(other)) {
		using std::copy; copy(other.begin(), other.end(), this->begin());
	}

	template<class TT, class... Args>
	explicit static_array(array_ref<TT, 0, Args...> const& other)
	: array_alloc{}, ref{static_array::allocate(other.num_elements()), extensions(other)} {
		uninitialized_copy_(other.data_elements());
	}

	static_array(static_array const& other, allocator_type const& a)       // 5b
	: array_alloc{a}
	, ref{static_array::allocate(other.num_elements()), extensions(other)} {
		uninitialized_copy_(other.data_elements());
	}

	static_array(static_array const& o)                                    // 5b
	: array_alloc{o.get_allocator()}
	, ref{static_array::allocate(o.num_elements(), o.data_elements()), {}} {
		uninitialized_copy(o.data_elements());
	}

	static_array(static_array&& o) noexcept :  // it is private because it is a valid operation for derived classes //5b
		array_alloc{o.get_allocator()},
		ref{static_array::allocate(static_cast<typename std::allocator_traits<allocator_type>::size_type>(o.num_elements()), o.data_elements()), o.extensions()} {
		uninitialized_move(o.data_elements());
	}
	template<class It> static auto distance(It a, It b) {using std::distance; return distance(a, b);}

 protected:
	void deallocate() {  // TODO(correaa) : move this to array_allocator
		if(this->num_elements()) {
			alloc_traits::deallocate(this->alloc(), this->base_, static_cast<typename alloc_traits::size_type>(this->num_elements()));
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
	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element const>;

	friend
	#if not defined(__NVCC__) and not defined(__INTEL_COMPILER)
	constexpr
	#endif
	auto get_allocator(static_array const& s) -> allocator_type {return s.get_allocator();}

	       constexpr auto base()                 &    -> typename static_array::element_ptr       {return ref::base();}
	       constexpr auto base()            const&    -> typename static_array::element_const_ptr {return ref::base();}
	friend constexpr auto base(static_array      & s) -> typename static_array::element_ptr       {return    s.base();}
	friend constexpr auto base(static_array const& s) -> typename static_array::element_const_ptr {return    s.base();}

	       constexpr auto origin()                 &    -> typename static_array::element_ptr       {return ref::origin();}
	       constexpr auto origin()            const&    -> typename static_array::element_const_ptr {return ref::origin();}
	friend constexpr auto origin(static_array      & s) -> typename static_array::element_ptr       {return    s.origin();}
	friend constexpr auto origin(static_array const& s) -> typename static_array::element_const_ptr {return    s.origin();}

	constexpr explicit operator typename std::iterator_traits<typename static_array::element_const_ptr>::reference() const& {
		return *(this->base_);
	}
	constexpr explicit operator typename std::add_rvalue_reference<typename std::iterator_traits<typename static_array::element_ptr>::reference>::type()&& {
		return *(this->base_);
	}
	constexpr explicit operator typename std::iterator_traits<typename static_array::element_ptr>::reference()& {
		return *(this->base_);
	}

//  constexpr auto rotated(dimensionality_type d) const& {
//  	typename static_array::layout_t new_layout = *this;
//  	new_layout.rotate(d);
//  	return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
//  }
	constexpr auto rotated() const& {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}

//  constexpr auto rotated(dimensionality_type d)& {
//  	typename static_array::layout_t new_layout = *this;
//  	new_layout.rotate(d);
//  	return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
//  }
	constexpr auto rotated() & {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

//  constexpr auto rotated(dimensionality_type d) && {
//  	typename static_array::layout_t new_layout = *this;
//  	new_layout.rotate(d);
//  	return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
//  }
	constexpr auto rotated() && {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

	friend constexpr auto rotated(static_array&       s) -> decltype(auto) {return s.rotated();}
	friend constexpr auto rotated(static_array const& s) -> decltype(auto) {return s.rotated();}

	constexpr auto unrotated(dimensionality_type d) const& {
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto unrotated(dimensionality_type d)& {
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

 private:
	constexpr auto unrotated_aux() {
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate();
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}

 public:
	constexpr auto unrotated()      & {return unrotated_aux();}
	constexpr auto unrotated() const& {return unrotated_aux().as_const();}

	friend constexpr auto unrotated(static_array      & self) -> decltype(auto) {return self.unrotated();}
	friend constexpr auto unrotated(static_array const& self) -> decltype(auto) {return self.unrotated();}

//  constexpr auto operator<<(dimensionality_type d)       -> decltype(auto) {return   rotated(d);}
//  constexpr auto operator>>(dimensionality_type d)       -> decltype(auto) {return unrotated(d);}

//  constexpr auto operator<<(dimensionality_type d) const -> decltype(auto) {return   rotated(d);}
//  constexpr auto operator>>(dimensionality_type d) const -> decltype(auto) {return unrotated(d);}

	constexpr auto operator=(static_array const& other) -> static_array& {
		assert( extensions(other) == static_array::extensions() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		if(this == &other) {return *this;}  // lints (cert-oop54-cpp) : handle self-assignment properly
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

 private:
	constexpr auto equal_extensions_if(std::true_type  /*true */, static_array const& other   ) {return this->extensions() == extensions(other);}
	constexpr auto equal_extensions_if(std::false_type /*false*/, static_array const&/*other*/) {return true;}

 public:
	constexpr auto operator=(static_array&& other) noexcept -> static_array& {
		assert( equal_extensions_if(std::integral_constant<bool, (static_array::rank_v != 0)>{}, other) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : allow a constexpr-friendly assert
		adl_move(other.data_elements(), other.data_elements() + other.num_elements(), this->data_elements());  // there is no std::move_n algorithm
		return *this;
	}

	template<class TT, class... As, class = std::enable_if_t<std::is_assignable<typename static_array::element_ref, TT>{}> >
	auto operator=(static_array<TT, 0, As...> const& other)& -> static_array& {
		assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	constexpr explicit operator basic_array<typename static_array::value_type, 0, typename static_array::element_const_ptr, typename static_array::layout_t>()& {
		return this->template static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>();
		  // return static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {  // case D = 0
		ref::serialize(ar, version);
	}
};

template<typename T, class Alloc>
struct array<T, 0, Alloc> : static_array<T, 0, Alloc>{
	using static_ = static_array<T, 0, Alloc>;
	using static_::static_;

	auto reextent(typename array::extensions_type const& /*empty_extensions*/) -> array& {
		return *this;
	}

	// NOLINTNEXTLINE(runtime/operator)
	auto operator&()     && -> array      * = delete;  // NOLINT(google-runtime-operator) : delete operator&& defined in base class to avoid taking address of temporary
	// auto operator&()      & -> array      *{return this;}
	// auto operator&() const& -> array const*{return this;}
};

template<class T, dimensionality_type D, class Alloc>
struct array : static_array<T, D, Alloc> {
	using static_ = static_array<T, D, Alloc>;
	static_assert(
		   std::is_same<typename array::alloc_traits::value_type, T   >{}
		or std::is_same<typename array::alloc_traits::value_type, void>{}, "!"
	);

 public:
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&()     && -> array      * = delete;  // NOLINT(google-runtime-operator) : delete operator&& defined in base class to avoid taking address of temporary
	//  auto operator&()      & -> array      *{return this;}
	//  auto operator&() const& -> array const*{return this;}

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {  // NOLINT(fuchsia-default-arguments-declarations) version is used for threshold of big vs small data
		using AT = multi::archive_traits<Archive>;
		auto extensions_ = this->extensions();
		ar &                   AT::make_nvp("extensions", extensions_);
	//  ar & boost::serialization::make_nvp("extensions", extensions);
	//  ar &        cereal       ::make_nvp("extensions", extensions);
	//  ar &        BOOST_SERIALIZATION_NVP              (extensions);
	//  ar &                     CEREAL_NVP              (extensions);
	//  ar &                                              extensions ;
		if(extensions_ != this->extensions()) {clear(); this->reextent(extensions_);}
		static_::serialize(ar, version);
	}

	using static_::static_;
	using typename static_::value_type;

	array() = default;
	array(array const&) = default;
	auto reshape(typename array::extensions_type x) & -> array& {
		typename array::layout_t new_layout{x};
		assert( new_layout.num_elements() == this->num_elements() );
		static_cast<typename array::layout_t&>(*this) = new_layout;
		return *this;
	}
	auto clear() noexcept -> array& {
		static_::clear();
		return *this;
	}
	friend auto clear(array& self) noexcept -> array& {return self.clear();}

	friend auto data_elements(array const& self) {return self.data_elements();}
	friend auto data_elements(array      & self) {return self.data_elements();}
	friend auto data_elements(array     && self) {return std::move(self).data_elements();}

	auto move() & -> basic_array<typename array::element, D, multi::move_ptr<typename array::element> >{
		basic_array<typename array::element, D, multi::move_ptr<typename array::element> >
		ret = multi::static_array_cast<typename array::element, multi::move_ptr<typename array::element>>(*this);
		layout_t<D>::operator=({});
		return ret;
	}
	friend auto move(array& self) -> basic_array<typename array::element, D, multi::move_ptr<typename array::element> >{
		return self.move();
	}

	array(array&& o, typename array::allocator_type const& a) noexcept : static_{std::move(o), a} {}
	array(array&& o) noexcept : array{std::move(o), o.get_allocator()} {}

	friend
	#if not defined(__NVCC__) and not defined(__NVCOMPILER) and not defined(__INTEL_COMPILER)
	constexpr
	#endif
	auto get_allocator(array const& s) -> typename array::allocator_type {return s.get_allocator();}

 private:
	template<class TrueType>
	void swap_allocator_if(TrueType&&     /*T*/, typename array::allocator_type&  source  ) {using std::swap; swap(this->alloc(), source);}
	void swap_allocator_if(std::true_type /*T*/, typename array::allocator_type&  source  ) {using std::swap; swap(this->alloc(), source);}
	void swap_allocator_if(std::false_type/*F*/, typename array::allocator_type&/*source*/) {}

 public:
	void swap(array& other) noexcept {
		using std::swap;
		swap_allocator_if(typename std::allocator_traits<typename array::allocator_type>::propagate_on_container_swap{}, other.alloc());
		swap(this->base_, other.base_);
		swap(
			static_cast<typename array::layout_t&>(*this),
			static_cast<typename array::layout_t&>(other)
		);
	}

#ifndef NOEXCEPT_ASSIGNMENT

 private:
	template<class TrueType>
	void move_allocator_if(TrueType&&     /*T*/, typename array::allocator_type&&  source  ) {this->alloc() = std::move(source);}
	void move_allocator_if(std::true_type /*T*/, typename array::allocator_type&&  source  ) {this->alloc() = std::move(source);}
	void move_allocator_if(std::false_type/*F*/, typename array::allocator_type&&/*source*/) {}

 public:
	auto operator=(array&& other) noexcept -> array& {
		clear();
	//  this->base_ = std::exchange(other.base_, nullptr);  // final null assigment shouldn't be necessary?
		this->base_ = other.base_;
		move_allocator_if(typename std::allocator_traits<typename array::allocator_type>::propagate_on_container_move_assignment{}, std::move(other.alloc()));
		// this->alloc_ = std::move(other.alloc_);
		static_cast<typename array::layout_t&>(*this) = std::exchange(static_cast<typename array::layout_t&>(other), {});
		return *this;
	}

 private:
	template<class TrueType>
	void copy_allocator_if(TrueType&&     /*T*/, typename array::allocator_type const&  source  ) {this->alloc() = source;}
	void copy_allocator_if(std::true_type /*T*/, typename array::allocator_type const&  source  ) {this->alloc() = source;}
	void copy_allocator_if(std::false_type/*F*/, typename array::allocator_type const&/*source*/) {}

 public:
	auto operator=(array const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			if(this == &other) {return *this;}  // required by cert-oop54-cpp
			copy_allocator_if(typename std::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment{}, other.alloc());
			static_::operator=(other);
		} else {
			clear();
			copy_allocator_if(typename std::allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment{}, other.alloc());
			static_cast<typename array::layout_t&>(*this) = static_cast<typename array::layout_t const&>(other);
			array::allocate();
			array::uninitialized_copy_elements(other.data_elements());
			// operator=(array{other}); // calls operator=(array&&)
		}
		return *this;
	}
#else
	auto operator=(array o) noexcept -> array& {return swap(o), *this;}
#endif

	template<
		class Range,
		class = decltype(std::declval<static_&>().operator=(std::declval<Range&&>())),
		std::enable_if_t<not std::is_base_of<array, std::decay_t<Range>>{}, int> = 0
	>
	auto operator=(Range&& other) ->array& {  // TODO(correaa) : check that LHS is not read-only?
		if(array::extensions() == other.extensions()) {
			static_::operator=(other);
		} else {
			array tmp(other);
			operator=(std::move(tmp));  // operator=(array{other}); produces an error in nvcc 11.2
		}
		return *this;
	}

	template<class TT, class... Args>
	auto operator=(multi::basic_array<TT, D, Args...> const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			static_::operator=(other);  // TODO(correaa) : protect for self assigment
		} else {
			operator=(array{other});
		}
		return *this;
	}

	auto operator=(basic_array<T, D, multi::move_ptr<typename array::element, typename array::element_ptr>>& other) -> array& {
		if(other.layout() != this->layout()) {
			array::operator=(other.template static_array_cast<typename array::element, typename array::element_ptr>());
			return *this;
		}
		if(this->base_ != other.base_) {
			other.base_ = nullptr;
		}
		return *this;
	}
	friend void swap(array& a, array& b) {a.swap(b);}
	void assign(typename array::extensions_type x, typename array::element const& e) {
		if(array::extensions() == x) {
			adl_fill_n(this->base_, this->num_elements(), e);
		} else {
			this->clear();
			(*this).array::layout_t::operator=(layout_t<D>{x});
			this->base_ = this->static_::array_alloc::allocate(this->num_elements());
			adl_alloc_uninitialized_fill_n(this->alloc(), this->base_, this->num_elements(), e);  //  recursive_uninitialized_fill<dimensionality>(alloc(), begin(), end(), e);
		}
	}

	template<class It>
	auto assign(It first, It last) -> array& {
		using std::next; using std::all_of;
		if(adl_distance(first, last) == array::size()) {  // and multi::extensions(*first) == multi::extensions(*array::begin())){
			static_::ref::assign(first);
		} else {
			this->operator=(array(first, last));
		}
		return *this;
	}
	void assign(std::initializer_list<typename array::value_type> il) {assign(il.begin(), il.end());}

	template<class Range> auto assign(Range&& r) &
	->decltype(assign(adl_begin(r), adl_end(r))) {
		return assign(adl_begin(r), adl_end(r)); }

	auto operator=(std::initializer_list<typename array::value_type> il) -> array& {
		assign(il.begin(), il.end()); return *this;
	}

	template<class... TTs>
	[[deprecated]] auto reextent(std::tuple<TTs...> const& other) -> array& {
		return reextent(std::apply([](auto const&... es) {return typename array::extensions_type(es...);}, other));  // paren is important here ext_type(...) for allow narrowing casts
	}

	auto reextent(typename array::extensions_type const& x) -> array& {
		if(x == this->extensions()) {
			return *this;
		}
#if 0
		array tmp(x, this->get_allocator());  // TODO(correaa) opportunity missed to use hint allocation
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
#else
		auto&& tmp = typename array::ref{
			this->static_::array_alloc::allocate(
				static_cast<typename std::allocator_traits<typename array::allocator_type>::size_type>(
					typename array::layout_t{x}.num_elements()
				),
				this->data_elements()  // used as hint
			),
			x
		};
		if(not std::is_trivially_default_constructible<typename array::element>{}) {  // TODO(correaa) convert into constexpr if
			adl_alloc_uninitialized_value_construct_n(this->alloc(), tmp.data_elements(), tmp.num_elements());
		}
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);  // TODO(correaa) : use (and implement) `.move();`
		this->destroy();
		this->deallocate();
		this->base_ = tmp.base();
		(*this).array::layout_t::operator=(tmp.layout());
#endif
		return *this;
	}
	auto reextent(typename array::extensions_type const& x, typename array::element const& e) -> array& {
		if(x == this->extensions()) {
			return *this;
		}
#if 0
		array tmp(x, e, this->get_allocator());  // TODO(correaa) opportunity missed to use hint allocation
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
#else  // implementation with hint
		auto&& tmp = typename array::ref{this->static_::array_alloc::allocate(static_cast<typename std::allocator_traits<typename array::allocator_type>::size_type>(typename array::layout_t{x}.num_elements()), this->data_elements()), x};
		this->uninitialized_fill_n(tmp.data_elements(), static_cast<typename std::allocator_traits<typename array::allocator_type>::size_type>(tmp.num_elements()), e);
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		this->destroy();
		this->deallocate();
		this->base_ = tmp.base();  // TODO(correaa) : use (and implement) `.move();`
		(*this).array::layout_t::operator=(tmp.layout());
#endif
		return *this;
	}
	template<class... Ts> constexpr auto reindex(Ts... a)&& -> array&& {array::layout_t::reindex(a...); return std::move(*this);}
	template<class... Ts> constexpr auto reindex(Ts... a) & -> array & {array::layout_t::reindex(a...); return           *this ;}

	~array() = default;
};

#if defined(__cpp_deduction_guides)

#define IL std::initializer_list

template<class T> static_array(IL<T>                ) -> static_array<T, 1>;
template<class T> static_array(IL<IL<T>>            ) -> static_array<T, 2>;
template<class T> static_array(IL<IL<IL<T>>>        ) -> static_array<T, 3>;
template<class T> static_array(IL<IL<IL<IL<T>>>>    ) -> static_array<T, 4>;
template<class T> static_array(IL<IL<IL<IL<IL<T>>>>>) -> static_array<T, 5>;

template<class T>        array(IL<T>                ) ->        array<T, 1>;
template<class T>        array(IL<IL<T>>            ) ->        array<T, 2>;
template<class T>        array(IL<IL<IL<T>>>        ) ->        array<T, 3>;
template<class T>        array(IL<IL<IL<IL<T>>>>    ) ->        array<T, 4>;
template<class T>        array(IL<IL<IL<IL<IL<T>>>>>) ->        array<T, 5>;

#undef IL

template<class T>        array(T[]                  ) ->        array<T, 1>;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

//  vvv these are necessary to catch {n, m, ...} notation (or single integer notation)
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<0>, T) -> array<T, 0>;
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<1>, T) -> array<T, 1>;
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<2>, T) -> array<T, 2>;
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<3>, T) -> array<T, 3>;
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<4>, T) -> array<T, 4>;
template<class T, class = std::enable_if_t<not is_allocator<T>{}> > array(iextensions<5>, T) -> array<T, 5>;

// generalization, will not work with naked {n, m, ...} notation (or single integer notation)
template<dimensionality_type D, class T, class = std::enable_if_t<not is_allocator<T>{}> >
array(iextensions<D>, T) -> array<T, D>;

template<class T, class MR, class A = memory::allocator<T, MR>> array(extensions_t<1>, T, MR*) -> array<T, 1, A>;
template<class T, class MR, class A = memory::allocator<T, MR>> array(extensions_t<2>, T, MR*) -> array<T, 2, A>;
template<class T, class MR, class A = memory::allocator<T, MR>> array(extensions_t<3>, T, MR*) -> array<T, 3, A>;
template<class T, class MR, class A = memory::allocator<T, MR>> array(extensions_t<4>, T, MR*) -> array<T, 4, A>;
template<class T, class MR, class A = memory::allocator<T, MR>> array(extensions_t<5>, T, MR*) -> array<T, 5, A>;

template<class MatrixRef, class DT = typename MatrixRef::decay_type, class T = typename DT::element, dimensionality_type D = typename DT::rank{}, class Alloc = typename DT::allocator_type>
array(MatrixRef)->array<T, D, Alloc>;

template<typename T, dimensionality_type D, typename P> array(basic_array<T, D, P>)->array<T, D>;

#endif  // ends defined(__cpp_deduction_guides)

template <class T, std::size_t N>
auto decay(const T(&t)[N]) noexcept -> multi::array<typename std::remove_all_extents<T[N]>::type, std::rank<T[N]>{}>{return multi::array_cref<typename std::remove_all_extents<T[N]>::type, std::rank<T[N]>{}>(data_elements(t), extensions(t));}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N>
struct array_traits<T[N], void, void> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using reference = T&;
	using element = std::remove_all_extents_t<T[N]>;  //  NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using decay_type = multi::array<T, 1>;
};

}  // end namespace boost::multi

#if(__cpp_lib_memory_resource >= 201603)
namespace boost::multi::pmr {

template <class T, boost::multi::dimensionality_type D>
using array = boost::multi::array<T, D, std::pmr::polymorphic_allocator<T>>;

}  // end namespace boost::multi::pmr
#endif

namespace boost::serialization {

template<typename T, boost::multi::dimensionality_type D, class A>
struct version< boost::multi::array<T, D, A> > {
	using type = std::integral_constant<int, MULTI_SERIALIZATION_ARRAY_VERSION>;  // typedef mpl::int_<1> type;
	enum { value = type::value };
};

}  // end namespace boost::serialization

#endif
