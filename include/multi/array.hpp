// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_ARRAY_HPP
#define MULTI_ARRAY_HPP

#include "./array_ref.hpp"
#include "./config/NO_UNIQUE_ADDRESS.hpp"

#include "./detail/adl.hpp"
#include "./detail/memory.hpp"
#include "./detail/type_traits.hpp"

// #include "./memory/allocator.hpp"

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

	using allocator_traits = typename multi::allocator_traits<allocator_type>;
	using size_type_ = typename allocator_traits::size_type;
	using pointer_   = typename allocator_traits::pointer;

 protected:
	auto alloc()      & -> allocator_type      & {return alloc_;}
	auto alloc() const& -> allocator_type const& {return alloc_;}

	array_allocator() = default;
	explicit array_allocator(allocator_type const& alloc) : alloc_{alloc} {}

	auto allocate(size_type_ n) -> pointer_ {
		return n?allocator_traits::allocate(alloc_, n):pointer_{nullptr};
	}
	auto allocate(size_type_ n, typename allocator_traits::const_void_pointer hint) -> pointer_ {
		return n?allocator_traits::allocate(alloc_, n, hint):pointer_{nullptr};
	}

	auto uninitialized_fill_n(pointer_ first, size_type_ count, typename allocator_traits::value_type value) {
		return adl_alloc_uninitialized_fill_n(alloc_, first, count, value);
	}
	template<typename It>
	auto uninitialized_copy_n(It first, size_type count, pointer_ d_first) {
		if constexpr(std::is_trivial_v<typename std::iterator_traits<pointer_>::value_type>) {
			return                     adl_copy_n(        first, count, d_first);
		} else {
			return adl_alloc_uninitialized_copy_n(alloc_, first, count, d_first);
		}
	}
	template<typename It>
	auto destroy_n(It first, size_type n) {return adl_alloc_destroy_n(this->alloc(), first, n);}

 public:
	constexpr auto get_allocator() const -> allocator_type {return alloc_;}
};

template<class T, dimensionality_type D, class Alloc = std::allocator<T>>
struct static_array  // NOLINT(fuchsia-multiple-inheritance) : multiple inheritance used for composition
: protected array_allocator<Alloc>
, public array_ref<T, D, typename allocator_traits<Alloc>::pointer>
, boost::multi::random_iterable<static_array<T, D, Alloc>> {
 protected:
	using array_alloc = array_allocator<Alloc>;

 public:
	static_assert( std::is_same<typename allocator_traits<Alloc>::value_type, typename static_array::element>{},
		"allocator value type must match array value type");

	using array_alloc::get_allocator;
	using allocator_type = typename array_allocator<Alloc>::allocator_type;
	using decay_type = array<T, D, Alloc>;
	using layout_type = typename array_ref<T, D, typename allocator_traits<Alloc>::pointer>::layout_type;

	using ref = array_ref<T, D, typename allocator_traits<typename allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;

 protected:
	using alloc_traits = typename multi::allocator_traits<allocator_type>;

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

	void destroy() {
		if constexpr(not std::is_trivially_destructible_v<typename static_array::element>) {
			array_alloc::destroy_n(this->data_elements(), this->num_elements());
		}
	}

	void allocate() {this->base_ = array_alloc::allocate(static_cast<typename allocator_traits<typename static_array::allocator_type>::size_type>(static_array::num_elements()));}

 public:
	using value_type = typename std::conditional<
		(D > 1),  // this parenthesis is needed
		array<typename static_array::element, D-1, allocator_type>,
		typename static_array::element
	>::type;

	using typename ref::size_type;
	using typename ref::difference_type;
	explicit static_array(allocator_type const& alloc) : array_alloc{alloc} {}

	using ref::operator();
//  HD constexpr auto operator()()      & -> decltype(auto) {return ref::operator()();}
	HD constexpr auto operator()()     && -> decltype(auto) {return ref::element_moved();}
//  HD constexpr auto operator()() const& -> decltype(auto) {return ref::operator()();}

	using ref::take;
	constexpr auto take(difference_type n) && -> decltype(auto) {return ref::take(n).element_moved();}

	using ref::drop;
	constexpr auto drop(difference_type n) && -> decltype(auto) {return ref::drop(n).element_moved();}

 protected:
	static_array(static_array&& other, allocator_type const& alloc) noexcept          // 6b  TODO(correaa) move from array only
	: array_alloc{alloc}  // TODO(correaa) : handle allocation propagation here
	, ref{other.base_, other.extensions()} {
		other.layout_mutable() = {};
	//  other.ref::layout_t::operator=({});
		other.base_ = nullptr;
	}

	static_array(static_array&& other) noexcept
	: static_array(std::move(other), allocator_type{}) {}  // 6b

 public:
	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // decltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>
	// analogous to std::vector::vector (5) https://en.cppreference.com/w/cpp/container/vector/vector
	static_array(It first, It last, allocator_type const& alloc)
	: array_alloc{alloc}
	, ref {
		array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(layout_type {index_extension {adl_distance(first, last)}*multi::extensions(*first)}.num_elements())),
		index_extension {adl_distance(first, last)}*multi::extensions(*first)
	} {
		adl_alloc_uninitialized_copy(static_array::alloc(), first, last, ref::begin());
	}

	template<class It, class = typename std::iterator_traits<std::decay_t<It>>::difference_type>  // decltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>
	// analogous to std::vector::vector (5) https://en.cppreference.com/w/cpp/container/vector/vector
	static_array(It first, It last) : static_array(first, last, allocator_type{}) {}

	template<
		class Range, class = std::enable_if_t<not std::is_base_of<static_array, std::decay_t<Range>>{}>,
		class = decltype(/*static_array*/(std::declval<Range const&>().begin() - std::declval<Range const&>().end())),  // instantiation of static_array here gives a compiler error in 11.0, partially defined type?
		class = std::enable_if_t<not is_basic_array<Range const&>{}>
	>
	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions // NOLINTNEXTLINE(runtime/explicit)
	static_array(Range const& rng)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array{rng.begin(), rng.end()} {}

	template<class TT>
	auto uninitialized_fill_elements(TT const& value) {
		return array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), value);
	}

	// vvv TODO(correaa) : check if really necessary
	template<class TT, class... As>
	static_array(array_ref<TT, D, As...> const& other, allocator_type const& alloc)
	: array_alloc{alloc}
	, ref{
		array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(other.num_elements())),
		other.extensions()
	} {
		if constexpr(std::is_trivial_v<T>) {
			                    adl_copy_n(                       other.data_elements(), other.num_elements(), this->data_elements());
		} else {
			adl_alloc_uninitialized_copy_n(static_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
		}
	}

	template<class TT, class... As>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented  // NOLINTNEXTLINE(runtime/explicit)
	static_array(array_ref<TT, D, As...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array(other, allocator_type{}) {}
	// ^^^ TODO(correaa) : check if really necessary

	static_array(typename static_array::extensions_type extensions, typename static_array::element const& elem, allocator_type const& alloc)  // 2
	: array_alloc{alloc}
	, ref{array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), nullptr), extensions} {
		array_alloc::uninitialized_fill_n(this->data_elements(), static_cast<typename allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);
	}

	template<class Element, std::enable_if_t<std::is_convertible<Element, typename static_array::element>{} and (D == 0), int> = 0>
	explicit static_array(Element const& elem, allocator_type const& alloc)
	: static_array(typename static_array::extensions_type{}, elem, alloc) {}

	static_array(typename static_array::extensions_type extensions, typename static_array::element const& elem)  // 2
	: array_alloc{}
	, ref{array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), nullptr), extensions} {
		array_alloc::uninitialized_fill_n(this->base(), static_cast<typename allocator_traits<allocator_type>::size_type>(this->num_elements()), elem);
	}

	explicit static_array(typename static_array::extensions_type extensions, typename allocator_traits<Alloc>::const_void_pointer hint)
	: array_alloc{}
	, ref{array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), hint), extensions} {}

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)  // 3
	= delete;

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value)  // 3
	= delete;

//  analgous to std::vector::vector ((4)) https://en.cppreference.com/w/cpp/container/vector/vector
	explicit static_array(typename static_array::extensions_type extensions, allocator_type const& alloc)
	: array_alloc{alloc}
	, ref{array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements())), extensions} {
		uninitialized_default_construct();
	}

	explicit static_array(typename static_array::extensions_type extensions)
	: static_array(extensions, allocator_type{}) {}

	template<class TT, class... Args,
		class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::basic_array<TT, D, Args...>::element>{}>,
		class = decltype(adl_copy(std::declval<multi::basic_array<TT, D, Args...> const&>().begin(), std::declval<multi::basic_array<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	static_array(multi::basic_array<TT, D, Args...> const& other, allocator_type const& alloc)
	: static_array(other.extensions(), alloc) {
		adl_uninitialized_copy(other.begin(), other.end(), this->begin());  // TODO(correaa): call this conditionally on T properties
	}

	template<class TT, class... Args,
	//  class = std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::basic_array<TT, D, Args...>::element>{}>,
		class = decltype(adl_copy(std::declval<multi::basic_array<TT, D, Args...> const&>().begin(), std::declval<multi::basic_array<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented  // NOLINTNEXTLINE(runtime/explicit)
	static_array(multi::basic_array<TT, D, Args...> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	: static_array(other, allocator_type{}) {}

	template<class TT, class... Args>
	explicit static_array(array_ref<TT, D, Args...>&& other)
	: array_alloc{}
	, ref{array_alloc::allocate(other.num_elements()), other.extensions()} {
		static_array::uninitialized_copy_elements(std::move(other).data_elements());
	}

	static_array(static_array const& other)                               // 5b
	: array_alloc{allocator_traits<Alloc>::select_on_container_copy_construction(other.alloc())}
	, ref{array_alloc::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), extensions(other)} {
		uninitialized_copy_elements(other.data_elements());
	}

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	static_array(std::initializer_list<typename static_array<T, D>::value_type> values)
	: static_array{static_array<T, D>(values.begin(), values.end())} {}  // construct all with default constructor and copy to special memory at the end

	static_array(
		std::initializer_list<typename static_array<T, D>::value_type> values,
		allocator_type const& alloc
	)
	: static_array{static_array<T, D>(values.begin(), values.end()), alloc} {}

	template<class TT, std::size_t N>
	constexpr explicit static_array(TT(&array)[N])  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backward compatibility
	: static_array(std::begin(array), std::end(array)) {}

	// template<class It> static auto distance(It a, It b) {
	// 	using std::distance;
	// 	return distance(a, b);
	// }

	constexpr auto begin() const& -> typename static_array::const_iterator {return ref:: begin();}
	constexpr auto end  () const& -> typename static_array::const_iterator {return ref:: end  ();}

	constexpr auto begin()     && -> typename static_array:: move_iterator {return ref::mbegin();}
	constexpr auto end  ()     && -> typename static_array:: move_iterator {return ref::mend  ();}

	constexpr auto begin()      & -> typename static_array::      iterator {return ref:: begin();}
	constexpr auto end  ()      & -> typename static_array::      iterator {return ref:: end  ();}

	constexpr auto operator[](index idx) const& -> typename static_array::const_reference {return ref::operator[](idx);}
	constexpr auto operator[](index idx)     && -> decltype(auto)                         {
		if constexpr(D == 1) {return std::move(ref::operator[](idx)       );}
		else                 {return           ref::operator[](idx).moved();}  // NOLINT(readability/braces)
	}
	constexpr auto operator[](index idx)      & -> typename static_array::      reference {return ref::operator[](idx);}

 protected:
	void deallocate() {
		if(this->num_elements()) {
			alloc_traits::deallocate(this->alloc(), this->base_, static_cast<typename alloc_traits::size_type>(this->num_elements()));
		}
	}
	void clear() noexcept {
		this->destroy();
		deallocate();
		this->layout_mutable() = {};
	}
	template<class... Indices>
	constexpr auto reindex(Indices... idxs) & -> static_array&  {
		static_array::layout_t::reindex(idxs...);
		return *this;
	}
	template<class... Indices>
	constexpr auto reindex(Indices... idxs) && -> static_array&& {
		reindex(idxs...);
		return std::move(*this);
	}

 public:
	static_array() = default;
	~static_array() noexcept {destroy(); deallocate();}

	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element const>;
	using element_move_ptr  = multi::move_ptr<typename static_array::element_ptr>;

	using reference = typename std::conditional<
		(D > 1),
		basic_array<typename static_array::element, D - 1, typename static_array::element_ptr>,
		typename std::conditional<
			D == 1,
			typename std::iterator_traits<typename static_array::element_ptr>::reference,
			void
		>::type
	>::type;
	using const_reference = typename std::conditional<
		(D > 1),
		basic_array<typename static_array::element, D - 1, typename static_array::element_const_ptr>,  // TODO(correaa) should be const_reference, but doesn't work witn rangev3?
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
	auto get_allocator(static_array const& self) -> allocator_type {return self.get_allocator();}

	       HD constexpr auto data_elements()            const& ->                        element_const_ptr {return this->base_;}
	       HD constexpr auto data_elements()                 & -> typename static_array::element_ptr       {return this->base_;}
	       HD constexpr auto data_elements()                && -> typename static_array::element_move_ptr  {return std::make_move_iterator(this->base_);}
	friend constexpr auto data_elements(static_array const& self) {return           self .data_elements();}
	friend constexpr auto data_elements(static_array      & self) {return           self .data_elements();}
	friend constexpr auto data_elements(static_array     && self) {return std::move(self).data_elements();}

	       constexpr auto base()                 &       -> typename static_array::element_ptr       {return ref::base();}
	       constexpr auto base()            const&       -> typename static_array::element_const_ptr {return typename static_array::element_const_ptr{ref::base()};}
	friend constexpr auto base(static_array      & self) -> typename static_array::element_ptr       {return self.base();}
	friend constexpr auto base(static_array const& self) -> typename static_array::element_const_ptr {return self.base();}

	       constexpr auto origin()                 &       -> typename static_array::element_ptr       {return ref::origin();}
	       constexpr auto origin()            const&       -> typename static_array::element_const_ptr {return ref::origin();}
	friend constexpr auto origin(static_array      & self) -> typename static_array::element_ptr       {return self.origin();}
	friend constexpr auto origin(static_array const& self) -> typename static_array::element_const_ptr {return self.origin();}

 private:
	constexpr auto rotated_aux() const {
		typename static_array::layout_t new_layout = this->layout();
		new_layout.rotate();
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}

 public:
	constexpr auto rotated() const& {return std::move(*this).rotated_aux();}
	constexpr auto rotated()      & {return std::move(*this).rotated_aux();}
	constexpr auto rotated()     && {return std::move(*this).rotated_aux();}

	friend constexpr auto rotated(static_array&       self) -> decltype(auto) {return self.rotated();}
	friend /*constexpr*/ auto rotated(static_array const& self) -> decltype(auto) {return self.rotated();}

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

	template<class TT, class... Args>
	auto operator=(multi::basic_array<TT, D, Args...> const& other) -> static_array& {
		ref::operator=(other);  // TODO(correaa) : protect for self assigment
		return *this;
	}
	auto operator=(static_array const& other) & -> static_array& {
		if(std::addressof(other) == this) {return *this;}  // cert-oop54-cpp
		assert( extensions(other) == static_array::extensions() );
		if(&other == this) {return *this;}  // lints (cert-oop54-cpp) : handle self-assignment properly
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
	void serialize(Archive& arxiv, const unsigned int version) {
		ref::serialize(arxiv, version);
	}
};

template<class T, class Alloc>
struct static_array<T, 0, Alloc>  // NOLINT(fuchsia-multiple-inheritance) : design
: protected array_allocator<Alloc>
, public array_ref<T, 0, typename allocator_traits<typename array_allocator<Alloc>::allocator_type>::pointer> {
	static_assert( std::is_same<typename allocator_traits<Alloc>::value_type, typename static_array::element>{},
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
	using alloc_traits = typename multi::allocator_traits<allocator_type>;
	using ref = array_ref<T, 0, typename allocator_traits<typename allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;

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
	constexpr explicit static_array(allocator_type const& alloc) : array_alloc{alloc} {}

 protected:
	constexpr static_array(static_array&& other, allocator_type const& alloc)  // 6b
	: array_alloc{alloc}
	, ref{other.base_, other.extensions()} {
		other.ref::layout_t::operator=({});
	}

 public:
	using ref::operator==;
	using ref::operator!=;

	static_array(typename static_array::extensions_type extensions, typename static_array::element const& elem, allocator_type const& alloc)  // 2
	: array_alloc{alloc}
	, ref(static_array::allocate(typename static_array::layout_t{extensions}.num_elements()), extensions) {
		uninitialized_fill(elem);
	}

	static_array(typename static_array::element_type const& elem, allocator_type const& alloc)
	: static_array(typename static_array::extensions_type{}, elem, alloc) {}

	auto uninitialized_fill(typename static_array::element const& elem) {
		array_alloc::uninitialized_fill_n(
			this->base_,
			static_cast<typename allocator_traits<allocator_type>::size_type>(this->num_elements()),
			elem
		);
	}

	static_array(
		typename static_array::extensions_type const& extensions,
		typename static_array::element const& elem
	)  // 2
	: array_alloc{}
	, ref(static_array::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(typename static_array::layout_t{extensions}.num_elements()), nullptr), extensions) {
		uninitialized_fill(elem);
	}

	static_array() : static_array(multi::iextensions<0>{}) {}

	explicit static_array(typename static_array::element const& elem)       // 2
	: static_array(multi::iextensions<0>{}, elem) {}

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, value_type>{}>>
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value, allocator_type const& alloc)  // 3
	: static_array(extension*extensions(value), alloc) {
		using std::fill; fill(this->begin(), this->end(), value);
	}
	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, value_type>{}>>
	explicit static_array(typename static_array::index_extension const& extension, ValueType const& value)  // 3  // TODO(correaa) : call other constructor (above)
	: static_array(extension*extensions(value)) {
		using std::fill; fill(this->begin(), this->end(), value);
	}

	explicit static_array(typename static_array::extensions_type const& extensions, allocator_type const& alloc)  // 3
	: array_alloc{alloc}
	, ref{static_array::allocate(typename static_array::layout_t{extensions}.num_elements()), extensions} {
		uninitialized_value_construct();
	}
	explicit static_array(typename static_array::extensions_type const& extensions)  // 3
	: static_array(extensions, allocator_type{}) {}

	template<class TT, class... Args>
	explicit static_array(multi::basic_array<TT, 0, Args...> const& other, allocator_type const& alloc)
	: array_alloc{alloc}
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

	static_array(static_array const& other, allocator_type const& alloc)   // 5b
	: array_alloc{alloc}
	, ref{static_array::allocate(other.num_elements()), extensions(other)} {
		uninitialized_copy_(other.data_elements());
	}

	static_array(static_array const& other)                                // 5b
	: array_alloc{other.get_allocator()}
	, ref{static_array::allocate(other.num_elements(), other.data_elements()), {}} {
		uninitialized_copy(other.data_elements());
	}

	static_array(static_array&& other) noexcept  // it is private because it is a valid operation for derived classes //5b
	: array_alloc{other.get_allocator()}
	, ref{static_array::allocate(static_cast<typename allocator_traits<allocator_type>::size_type>(other.num_elements()), other.data_elements()), other.extensions()} {
		uninitialized_move(other.data_elements());
	}
//  template<class It> static auto distance(It a, It b) {using std::distance; return distance(a, b);}

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
	auto get_allocator(static_array const& self) -> allocator_type {return self.get_allocator();}

	       constexpr auto base()                 &       -> typename static_array::element_ptr       {return ref::base();}
	       constexpr auto base()            const&       -> typename static_array::element_const_ptr {return ref::base();}
	friend constexpr auto base(static_array      & self) -> typename static_array::element_ptr       {return self.base();}
	friend constexpr auto base(static_array const& self) -> typename static_array::element_const_ptr {return self.base();}

	       constexpr auto origin()                 &       -> typename static_array::element_ptr       {return ref::origin();}
	       constexpr auto origin()            const&       -> typename static_array::element_const_ptr {return ref::origin();}
	friend constexpr auto origin(static_array      & self) -> typename static_array::element_ptr       {return self.origin();}
	friend constexpr auto origin(static_array const& self) -> typename static_array::element_const_ptr {return self.origin();}

	constexpr explicit operator typename std::iterator_traits<typename static_array::element_const_ptr>::reference() const& {
		return *(this->base_);
	}
	constexpr explicit operator typename std::add_rvalue_reference<typename std::iterator_traits<typename static_array::element_ptr>::reference>::type()&& {
		return *(this->base_);
	}
	constexpr explicit operator typename std::iterator_traits<typename static_array::element_ptr>::reference()& {
		return *(this->base_);
	}

	constexpr auto rotated() const& {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated()  & {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

	constexpr auto rotated() && {
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate();
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}

	friend constexpr auto rotated(static_array&       self) -> decltype(auto) {return self.rotated();}
	friend constexpr auto rotated(static_array const& self) -> decltype(auto) {return self.rotated();}

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

	constexpr explicit operator basic_array<value_type, 0, typename static_array::element_const_ptr, typename static_array::layout_t>()& {
		return this->template static_array_cast<value_type, typename static_array::element_const_ptr>();
		  // return static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}

	template<class Archive>
	void serialize(Archive& arxiv, const unsigned int version) {
		ref::serialize(arxiv, version);
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
	auto operator&()     && -> array      * = delete;       // NOLINT(google-runtime-operator) : delete operator&& defined in base class to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&()      & -> array      * {return this;}  // NOLINT(google-runtime-operator) : delete operator&& defined in base class to avoid taking address of temporary
	// NOLINTNEXTLINE(runtime/operator)
	auto operator&() const& -> array const* {return this;}  // NOLINT(google-runtime-operator) : delete operator&& defined in base class to avoid taking address of temporary

	friend auto sizes(array const& self) -> typename array::sizes_type {return self.sizes();}

	template<class Archive>
	void serialize(Archive& arxiv, const unsigned int version) {
		using AT = multi::archive_traits<Archive>;
		auto extensions_ = this->extensions();
		arxiv &                   AT::make_nvp("extensions", extensions_);
	//  arxiv & boost::serialization::make_nvp("extensions", extensions );
	//  arxiv &        cereal       ::make_nvp("extensions", extensions );
	//  arxiv &        BOOST_SERIALIZATION_NVP(              extensions );
	//  arxiv &                     CEREAL_NVP(              extensions );
	//  arxiv &                                              extensions  ;
		if(this->extensions() != extensions_) {clear(); this->reextent(extensions_);}
		static_::serialize(arxiv, version);
	}

	using static_::static_;
	using typename static_::value_type;

	// cppcheck-suppress noExplicitConstructor ; to allow assignment-like construction of nested arrays
	array(std::initializer_list<typename static_array<T, D>::value_type> ilv)
	: static_{static_array<T, D>(ilv.begin(), ilv.end())} {}  // construct all with default constructor and copy to special memory at the end

	array() = default;
	array(array const&) = default;

	auto reshape(typename array::extensions_type extensions) & -> array& {
		typename array::layout_t new_layout{extensions};  // TODO(correaa) implement move-reextent in terms of reshape
		assert( new_layout.num_elements() == this->num_elements() );
		this->layout_mutable() = new_layout;
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

	auto move() & -> basic_array<typename array::element, D, multi::move_ptr<typename array::element>> {
		basic_array<typename array::element, D, multi::move_ptr<typename array::element>>
		ret = multi::static_array_cast<typename array::element, multi::move_ptr<typename array::element>>(*this);
		layout_t<D>::operator=({});
		return ret;
	}
	friend auto move(array& self) -> basic_array<typename array::element, D, multi::move_ptr<typename array::element> >{
		return self.move();
	}

	array(array&& other, typename array::allocator_type const& alloc) noexcept : static_{std::move(other), alloc} {}
	array(array&& other) noexcept : array{std::move(other), other.get_allocator()} {}

	friend
	#if not defined(__NVCC__) and not defined(__NVCOMPILER) and not defined(__INTEL_COMPILER)
	constexpr
	#endif
	auto get_allocator(array const& self) -> typename array::allocator_type {return self.get_allocator();}

	void swap(array& other) noexcept {
		using std::swap;
		if constexpr(allocator_traits<typename array::allocator_type>::propagate_on_container_swap::value) {
			swap(this->alloc(), other.alloc());
		}
		swap(this->base_, other.base_);
		swap(
			this->layout_mutable(),
			other.layout_mutable()
		);
	}

#ifndef NOEXCEPT_ASSIGNMENT
	auto operator=(array&& other) noexcept -> array& {
		clear();
		this->base_ = other.base_;
		if constexpr(allocator_traits<typename array::allocator_type>::propagate_on_container_move_assignment::value) {
			this->alloc() = std::move(other.alloc());
		}
		this->layout_mutable() = std::exchange(other.layout_mutable(), {});
		return *this;
	}

	auto operator=(array const& other) -> array& {
		if(array::extensions() == other.extensions()) {
			if(this == &other) {return *this;}  // required by cert-oop54-cpp
			if constexpr(allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			static_::operator=(other);
		} else {
			clear();
			if constexpr(allocator_traits<typename array::allocator_type>::propagate_on_container_copy_assignment::value) {
				this->alloc() = other.alloc();
			}
			this->layout_mutable() = other.layout();
			array::allocate();
			array::uninitialized_copy_elements(other.data_elements());
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
		} else if(this->num_elements() == other.extensions().num_elements()) {
			reshape(other.extensions());
			static_::operator=(other);
		} else {
			operator=(static_cast<array>(std::forward<Range>(other)));
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

	friend void swap(array& self, array& other) {self.swap(other);}

	void assign(typename array::extensions_type extensions, typename array::element const& elem) {
		if(array::extensions() == extensions) {
			adl_fill_n(this->base_, this->num_elements(), elem);
		} else {
			this->clear();
			(*this).array::layout_t::operator=(layout_t<D>{extensions});
			this->base_ = this->static_::array_alloc::allocate(this->num_elements(), nullptr);
			adl_alloc_uninitialized_fill_n(this->alloc(), this->base_, this->num_elements(), elem);  // recursive_uninitialized_fill<dimensionality>(alloc(), begin(), end(), e);
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
	void assign(std::initializer_list<value_type> values) {assign(values.begin(), values.end());}

	template<class Range> auto assign(Range&& other) &
	->decltype(assign(adl_begin(other), adl_end(other))) {  // TODO(correaa) use forward
		return assign(adl_begin(other), adl_end(other)); }

	auto operator=(std::initializer_list<value_type> values) -> array& {
		assign(values.begin(), values.end());
		return *this;
	}

	template<class... TTs>
	[[deprecated]] auto reextent(std::tuple<TTs...> const& other) -> array& {
		return reextent(
			std::apply([](auto const&... extensions) {return typename array::extensions_type(extensions...);}, other)
		);  // paren is important here ext_type(...) for allow narrowing casts ^^^
	}

	auto reextent(typename array::extensions_type const& extensions) && -> array& {
		if(extensions == this->extensions()) {return *this;}
		this->destroy();
		this->deallocate();
		this->layout_mutable() = typename array::layout_t{extensions};
		this->base_ = this->static_::array_alloc::allocate(
			static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(
				typename array::layout_t{extensions}.num_elements()
			),
			this->data_elements()  // used as hint
		);
		if constexpr(not std::is_trivially_default_constructible<typename array::element>{}) {  // TODO(correaa) convert into constexpr if
			adl_alloc_uninitialized_value_construct_n(this->alloc(), this->base_, this->num_elements());
		}

		return *this;
	}

	auto reextent(typename array::extensions_type const& extensions)  & -> array& {
		if(extensions == this->extensions()) {return *this;}
#if 0
		array tmp(x, this->get_allocator());  // TODO(correaa) opportunity missed to use hint allocation
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
#else
		auto&& tmp = typename array::ref{
			this->static_::array_alloc::allocate(
				static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(
					typename array::layout_t{extensions}.num_elements()
				),
				this->data_elements()  // used as hint
			),
			extensions
		};
		if constexpr(not std::is_trivially_default_constructible<typename array::element>{}) {  // TODO(correaa) convert into constexpr if
			adl_alloc_uninitialized_value_construct_n(this->alloc(), tmp.data_elements(), tmp.num_elements());
		}
		auto const is = intersection(this->extensions(), extensions);
		tmp.apply(is) = this->apply(is);  // TODO(correaa) : use (and implement) `.move();`
		this->destroy();
		this->deallocate();
		this->base_ = tmp.base();
		this->layout_mutable() = tmp.layout();
#endif
		return *this;
	}

	auto reextent(typename array::extensions_type const& extensions, typename array::element const& elem) && -> array& {
		if(extensions == this->extensions()) {return *this;}
		this->destroy();
		this->deallocate();
		this->layout_mutable() = typename array::layout_t{extensions};
		this->base_ = this->static_::array_alloc::allocate(
			static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(
				typename array::layout_t{extensions}.num_elements()
			),
			this->data_elements()  // used as hint
		);
		this->uninitialized_fill_n(this->base_, static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(this->num_elements()), elem);

		return *this;
	}

	auto reextent(typename array::extensions_type const& exs, typename array::element const& elem) & -> array& {
		if(exs == this->extensions()) {
			return *this;
		}
#if 0
		array tmp(x, e, this->get_allocator());  // TODO(correaa) opportunity missed to use hint allocation
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
#else  // implementation with hint
		auto&& tmp = typename array::ref{this->static_::array_alloc::allocate(
			static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(typename array::layout_t{exs}.num_elements()),
			this->data_elements()  // use as hint
		), exs};
		this->uninitialized_fill_n(tmp.data_elements(), static_cast<typename allocator_traits<typename array::allocator_type>::size_type>(tmp.num_elements()), elem);
		auto const is = intersection(this->extensions(), exs);
		tmp.apply(is) = this->apply(is);
		this->destroy();
		this->deallocate();
		this->base_ = tmp.base();  // TODO(correaa) : use (and implement) `.move();`
		this->layout_mutable() = tmp.layout();
	//  (*this).array::layout_t::operator=(tmp.layout());
#endif
		return *this;
	}
	template<class... Indices> constexpr auto reindex(Indices... idxs)&& -> array&& {this->layout_mutable().reindex(idxs...); return std::move(*this);}
	template<class... Indices> constexpr auto reindex(Indices... idxs) & -> array & {this->layout_mutable().reindex(idxs...); return           *this ;}

	~array() = default;
};

#if defined(__cpp_deduction_guides)

#define IL std::initializer_list  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing TODO(correaa) remove

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
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<0>, T) -> array<T, 0>;
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<1>, T) -> array<T, 1>;
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<2>, T) -> array<T, 2>;
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<3>, T) -> array<T, 3>;
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<4>, T) -> array<T, 4>;
template<class T, class = std::enable_if_t<not boost::multi::is_allocator_v<T>>> array(iextensions<5>, T) -> array<T, 5>;

// generalization, will not work with naked {n, m, ...} notation (or single integer notation)
template<dimensionality_type D, class T, class = std::enable_if_t<!boost::multi::is_allocator_v<T>> >
array(iextensions<D>, T) -> array<T, D>;

template<class MatrixRef, class DT = typename MatrixRef::decay_type, class T = typename DT::element, dimensionality_type D = DT::rank_v, class Alloc = typename DT::allocator_type>
array(MatrixRef)->array<T, D, Alloc>;

template<typename T, dimensionality_type D, typename P> array(basic_array<T, D, P>)->array<T, D>;

#endif  // ends defined(__cpp_deduction_guides)

template <class T, std::size_t N>
auto decay(const T(&arr)[N]) noexcept -> multi::array<typename std::remove_all_extents<T[N]>::type, std::rank_v<T[N]>> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	return multi::array_cref<typename std::remove_all_extents<T[N]>::type, std::rank_v<T[N]>>(data_elements(arr), extensions(arr));  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
}

template<class T, std::size_t N>
struct array_traits<T[N], void, void> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using reference = T&;
	using element = std::remove_all_extents_t<T[N]>;  //  NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	using decay_type = multi::array<T, 1>;
};

}  // end namespace boost::multi

#if defined(__cpp_lib_memory_resource) && (__cpp_lib_memory_resource >= 201603)
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
