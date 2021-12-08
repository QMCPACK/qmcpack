// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP

#include "./memory/pointer_traits.hpp"
#include "utility.hpp"

#include "./config/ASSERT.hpp"
#include "./config/DELETE.hpp"
#include "./config/MARK.hpp"
#include "./config/NODISCARD.hpp"

#include "./detail/layout.hpp"
#include "./detail/memory.hpp"     // for pointer_traits
#include "./detail/operators.hpp"  // for random_iterable
#include "./detail/types.hpp"      // for dimensionality_type

#if defined(__NVCC__)
#define HD __host__ __device__
#else
#define HD
#endif

#include<algorithm>   // fpr copy_n
#include<cstring>     // for memset in reinterpret_cast
#include<functional>  // for invoke
#include<iterator>    // for next
#include<memory>      // for pointer_traits
#include<utility>     // for forward

namespace std {

template<class T>
struct pointer_traits<std::move_iterator<T*>> : std::pointer_traits<T*> {
	template<class U> using rebind =
		std::conditional_t<
			std::is_const<U>{},
			U*,
			std::pointer_traits<std::move_iterator<U*>>
		>;
};

}  // end namespace std

namespace boost {
namespace multi {

template<class T> auto modify(T const& t) -> T& {return const_cast<T&>(t);}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) see what is this used for

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>>
struct basic_array;

template<typename T, dimensionality_type D, class A = std::allocator<T>> struct array;

template<class To, class From, std::enable_if_t<std::is_convertible<From, To>{}, int> = 0>
constexpr auto _implicit_cast(From&& f) -> To {return static_cast<To>(f);}

template<class To, class From, std::enable_if_t<std::is_constructible<To, From>{} and not std::is_convertible<From, To>{}, int> = 0>
constexpr auto _explicit_cast(From&& f) -> To {return static_cast<To>(f);}

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D> >
struct array_types : Layout{  // cppcheck-suppress syntaxError ; false positive in cppcheck
	using element = T;
	using element_type = element;  // this follows more closely https://en.cppreference.com/w/cpp/memory/pointer_traits
	using element_ptr = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element const>;
	using element_ref = typename std::iterator_traits<element_ptr>::reference;

	using layout_t = Layout;

	using value_type = typename std::conditional<
		(D > 1),
		array<element, D-1, typename multi::pointer_traits<element_ptr>::default_allocator_type>,
		element
	>::type;

	using reference = typename std::conditional<
		(D > 1),
		basic_array<element, D-1, element_ptr>,
		typename std::iterator_traits<element_ptr>::reference
	>::type;

	using const_reference = typename std::conditional<
		(D > 1),
		basic_array<element, D-1, element_const_ptr>,
		typename std::iterator_traits<element_const_ptr>::reference
	>::type;

	NODISCARD("") HD constexpr auto  base() const -> element_ptr       {return base_;}
	NODISCARD("")    constexpr auto cbase() const -> element_const_ptr {return base_;}

	NODISCARD("")    constexpr auto mbase() const&    -> element_ptr& {return base_;}

	friend           auto  base(array_types const& s) -> element_ptr  {return s.base();}

	       constexpr auto layout() const               -> layout_t const& {return *this;}
	friend constexpr auto layout(array_types const& s) -> layout_t const& {return s.layout();}

	       constexpr auto origin()           const&    -> decltype(auto) {return base_+Layout::origin();}
	friend constexpr auto origin(array_types const& s) -> decltype(auto) {return s.origin();}

 protected:
	using derived = basic_array<T, D, ElementPtr, Layout>;
	element_ptr base_;  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes,misc-non-private-member-variables-in-classes) : TODO(correaa) try to make it private, [static_]array needs mutation
	constexpr explicit array_types(std::nullptr_t np) : Layout{}, base_{np} {}

 public:
	array_types() = default;

//  #if defined(__NVCC__)
//  	__host__ __device__ // TODO check why this is necessary (nvcc 11), removing this gives a, trivial_device_copy D->H failed: cudaErrorLaunchFailure: unspecified launch failure
//  #endif
	constexpr array_types(layout_t const& l, element_ptr const& data)
	: Layout{l}, base_{data} {}

 protected:  // TODO(correaa) : find why this needs to be public and not protected or friend
	template<class ArrayTypes, typename = std::enable_if_t<not std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>
		, decltype(_explicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	>
	constexpr explicit array_types(ArrayTypes const& a) : Layout{a.layout()}, base_{a.base_} {}

	template<
		class ArrayTypes,
		typename = std::enable_if_t<not std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>,
		decltype(_implicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointers are implicitly convertible
	constexpr/*implct*/array_types(ArrayTypes const& a)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : inherit behavior of underlying pointer
	: Layout{a.layout()}, base_{a.base_} {}
	// ^^^ TODO(correaa) : call explicit from implicit, careful with infinite recursion

	template<
		typename ElementPtr2,
		typename = decltype(Layout{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().layout()}),
		typename = decltype(element_ptr{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().base_})
	>
	constexpr explicit array_types(array_types<T, D, ElementPtr2, Layout> const& other)
	: Layout{other.layout()}, base_{other.base_} {}

	template<class T2, dimensionality_type D2, class E2, class L2> friend struct array_types;
};

template<class Ref, class Layout>
struct basic_array_ptr  // NOLINT(fuchsia-multiple-inheritance) : to allow mixin CRTP
: private Ref  // TODO(correaa) : remove inheritance from Ref??
, boost::multi::iterator_facade<
	basic_array_ptr<Ref, Layout>, void, std::random_access_iterator_tag,
	Ref const&, typename Layout::difference_type
> {  //, boost::multi::totally_ordered2<basic_array_ptr<Ref, Layout>, void>
	~basic_array_ptr() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	auto operator=(basic_array_ptr&& other)  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	-> basic_array_ptr& {
		operator=(other);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}
	using pointer = Ref const*;
	using element_type = typename Ref::decay_type;
	using difference_type = typename Layout::difference_type;

	using value_type = element_type;
	using reference = Ref;
	using iterator_category = std::random_access_iterator_tag;

	constexpr explicit basic_array_ptr(std::nullptr_t p) : Ref{p} {}
	constexpr basic_array_ptr() : basic_array_ptr{nullptr} {}

	template<class, class> friend struct basic_array_ptr;

	constexpr basic_array_ptr(typename Ref::element_ptr p, layout_t<typename Ref::rank{}-1> l) : Ref{l, p} {}
	constexpr basic_array_ptr(typename Ref::element_ptr p, index_extensions<typename Ref::rank{}> e) : Ref{p, e} {}

	basic_array_ptr(basic_array_ptr      &&) noexcept = default;
	basic_array_ptr(basic_array_ptr const& )          = default;

	auto operator=(basic_array_ptr const& other) -> basic_array_ptr& {
		if(this == &other) {return *this;}  // lints(cert-oop54-cpp)
		this->base_ = other.base_;
		static_cast<Layout&>(*this) = other.layout();
		return *this;
	}
	constexpr explicit operator bool() const {return this->base_;}

	constexpr auto  dereference() const -> Ref {return Ref{this->layout(), this->base_};}

	HD constexpr auto  operator* () const -> Ref{return Ref{*this};}

	constexpr auto operator->() const -> Ref* {return  const_cast<basic_array_ptr*>(this);}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) find a better way without const_cast
	constexpr auto operator->()       -> Ref* {return  this;}

	constexpr auto  operator[](difference_type n) const -> Ref {return *(*this + n);}

	constexpr auto operator<(basic_array_ptr const& o) const -> bool{return distance_to(o) > 0;}

	constexpr basic_array_ptr(typename Ref::element_ptr p, Layout const& l) : Ref{l, p} {}

	template<typename T, dimensionality_type D, typename ElementPtr, class LLayout>
	friend struct basic_array;

	constexpr auto base() const {return this->base_;}

	friend constexpr auto base(basic_array_ptr const& self) {return self.base();}

	using Ref::base_;
	using Ref::layout;

	constexpr auto operator==(basic_array_ptr const& o) const -> bool{return base_==o.base_ and layout()==o.layout();}

	template<class O> constexpr auto operator==(O const& o) const -> bool{return base()==o->base() and layout() == o->layout();}
	template<class O> constexpr auto operator!=(O const& o) const -> bool{return not ((*this)==o);}

	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr auto operator==(O const& o, basic_array_ptr const& s) -> bool{return s.operator==(o);}
	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr auto operator!=(O const& o, basic_array_ptr const& s) -> bool{return not(o==s);}

 protected:
	constexpr void increment() {base_ += Ref::nelems();}
	constexpr void decrement() {base_ -= Ref::nelems();}

	constexpr void advance(difference_type n) {base_ += Ref::nelems()*n;}
	constexpr auto distance_to(basic_array_ptr const& other) const -> difference_type {
		assert( Ref::nelems() == other.Ref::nelems() and Ref::nelems() != 0 );
		assert( (other.base_ - base_)%Ref::nelems() == 0);
		assert( layout() == other.layout() );
		return (other.base_ - base_)/Ref::nelems();
	}

 public:
	constexpr auto operator+=(difference_type n) -> basic_array_ptr&{advance(n); return *this;}
};

template<class Element, dimensionality_type D, typename ElementPtr>
struct array_iterator;

template<class Element, dimensionality_type D, typename ElementPtr>
struct array_iterator
: boost::multi::iterator_facade<
	array_iterator<Element, D, ElementPtr>, void, std::random_access_iterator_tag,
	basic_array<Element, D-1, ElementPtr> const&, typename layout_t<D-1>::difference_type
>
, multi::decrementable<array_iterator<Element, D, ElementPtr>>
, multi::incrementable<array_iterator<Element, D, ElementPtr>>
, multi::affine<array_iterator<Element, D, ElementPtr>, multi::difference_type>
, multi::totally_ordered2<array_iterator<Element, D, ElementPtr>, void> {
	~array_iterator() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	auto operator=(array_iterator&&)  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	-> array_iterator& = default;

	array_iterator(array_iterator&&) noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	= default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	using difference_type = typename layout_t<D>::difference_type;
	using element = Element;
	using element_ptr = ElementPtr;
	using value_type = typename basic_array<element, D-1, element_ptr>::decay_type;

	using pointer   = basic_array<element, D-1, element_ptr>*;
	using reference = basic_array<element, D-1, element_ptr>;

	using iterator_category = std::random_access_iterator_tag;

	constexpr static dimensionality_type rank_v = D;
	using rank = std::integral_constant<dimensionality_type, D>;

	using ptr_type = basic_array_ptr<basic_array<element, D-1, element_ptr>, layout_t<D-1>>;
	using stride_type = index;

	explicit constexpr array_iterator(std::nullptr_t p) : ptr_{p} {}/*, stride_{1}*/
	constexpr array_iterator() : array_iterator{nullptr} {}

	template<class, dimensionality_type, class> friend struct array_iterator;

	template<
		class EElement, typename PPtr,
		decltype(_explicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().ptr_.base()))* = nullptr
	>
	constexpr explicit array_iterator(array_iterator<EElement, D, PPtr> const& o)
	: ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_} {}

	template<class EElement, typename PPtr,
		decltype(_implicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().ptr_.base()))* = nullptr
	>
	constexpr/*implct*/array_iterator(array_iterator<EElement, D, PPtr> const& o)   // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : propagate implicitness of pointer
	: ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_} {}
	// ^^^ TODO(correaa) : implement implcit in terms of explicit? be careful of infinite recursion

	array_iterator(array_iterator const&) = default;
	auto operator=(array_iterator const& other) -> array_iterator& = default;

	explicit constexpr operator bool() const {return static_cast<bool>(ptr_.base_);}
	HD constexpr auto operator*() const -> basic_array<element, D-1, element_ptr> {/*assert(*this);*/ return {*ptr_};}

	constexpr auto operator->() const -> decltype(auto) {/*assert(*this);*/ return ptr_;}

	HD constexpr auto operator+ (difference_type n) const -> array_iterator {array_iterator ret{*this}; ret+=n; return ret;}
	HD constexpr auto operator[](difference_type n) const -> basic_array<element, D-1, element_ptr> {return *((*this) + n);}

	constexpr auto operator==(array_iterator const& o) const -> bool {return ptr_==o.ptr_ and stride_==o.stride_ and ptr_.layout() == o.ptr_.layout();}
	constexpr auto operator< (array_iterator const& o) const -> bool {return distance_to(o) > 0;}

	explicit constexpr array_iterator(typename basic_array<element, D-1, element_ptr>::element_ptr p, layout_t<D-1> l, index stride)
	: ptr_{p, l}, stride_{stride} {}

	template<class, dimensionality_type, class, class> friend struct basic_array;

	template<class... As>
	HD constexpr auto operator()(index i, As... as) const -> decltype(auto) {return this->operator[](i)(as...); }
	HD constexpr auto operator()(index i          ) const -> decltype(auto) {return this->operator[](i)       ; }

 private:
	template<typename Tuple, std::size_t ... I>
	HD constexpr auto apply_impl(Tuple const& t, std::index_sequence<I...>/*012*/) const -> decltype(auto) {
		return this->operator()(std::get<I>(t)...);
	}

 public:
	template<typename Tuple> HD constexpr auto apply(Tuple const& t) const -> decltype(auto) {return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}

 private:
	ptr_type ptr_;
	stride_type stride_ = {1};  // nice non-zero default

	constexpr auto equal(array_iterator const& o) const -> bool {return ptr_==o.ptr_ and stride_==o.stride_;}
	constexpr void decrement() {ptr_.base_ -= stride_;}
	constexpr void advance(difference_type n) {ptr_.base_ += stride_*n;}
	constexpr auto distance_to(array_iterator const& other) const -> difference_type {
		assert( stride_ == other.stride_); assert( stride_ != 0 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		return (other.ptr_.base_ - ptr_.base_)/stride_;
	}

 public:
	       constexpr auto base()              const&    -> element_ptr {return ptr_.base_;}
	friend constexpr auto base(array_iterator const& s) -> element_ptr {return s.base();}

	       constexpr auto stride()              const&    -> stride_type {return   stride_;}
	friend constexpr auto stride(array_iterator const& s) -> stride_type {return s.stride_;}

	constexpr auto operator++() -> array_iterator& {ptr_.base_ += stride_; return *this;}
	constexpr auto operator--() -> array_iterator& {decrement(); return *this;}

	friend constexpr auto operator-(array_iterator const& self, array_iterator const& other) -> difference_type {
		assert(self.stride_ == other.stride_);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		assert(self.stride_ != 0);              // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		return (self.ptr_.base_ - other.ptr_.base_)/self.stride_;
	}

	constexpr auto operator+=(difference_type d) -> array_iterator& {advance(+d); return *this;}
	constexpr auto operator-=(difference_type d) -> array_iterator& {advance(-d); return *this;}
};

template<class It>
struct biiterator
: boost::multi::iterator_facade<
	biiterator<It>,
	typename std::iterator_traits<It>::value_type, std::random_access_iterator_tag,
	decltype(*(std::move((*std::declval<It>())).begin())), multi::difference_type
>
,	multi::affine<biiterator<It>, multi::difference_type>
,	multi::decrementable<biiterator<It>>
,	multi::incrementable<biiterator<It>>
,	multi::totally_ordered2<biiterator<It>, void> {
private:
	It me_ = {};
	std::ptrdiff_t pos_ = 0;
	std::ptrdiff_t stride_ = 1;

public:
	constexpr biiterator(It me, std::ptrdiff_t pos, std::ptrdiff_t stride)
	: me_{me}, pos_{pos}, stride_{stride} {}

	constexpr auto operator++() -> decltype(auto) {
		++pos_;
		if(pos_==stride_) {
			++me_;
			pos_ = 0;
		}
		return *this;
	}

	constexpr auto operator==(biiterator const& o) const -> bool {return me_==o.me_ and pos_==o.pos_;}
	constexpr auto operator+=(multi::difference_type n) -> biiterator& {me_ += n/stride_; pos_ += n%stride_; return *this;}

	constexpr auto operator*() const -> decltype(auto) {
		auto meb = std::move(*me_).begin();
		return meb[pos_];
	}

	using difference_type = std::ptrdiff_t;
	using reference = decltype(*std::declval<biiterator>());
	using value_type = std::decay_t<reference>;
	using pointer = value_type*;
	using iterator_category = std::random_access_iterator_tag;
};

template<class It>
auto ref(It begin, It end)
->multi::basic_array<typename It::element, It::rank_v, typename It::element_ptr> {
	return multi::basic_array<typename It::element, It::rank_v, typename It::element_ptr>{begin, end};
}

template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
struct basic_array
: multi::partially_ordered2<basic_array<T, D, ElementPtr, Layout>, void>
, array_types<T, D, ElementPtr, Layout> {
	using types = array_types<T, D, ElementPtr, Layout>;

	friend struct basic_array<typename types::element, Layout::rank_v + 1, typename types::element_ptr >;
	friend struct basic_array<typename types::element, Layout::rank_v + 1, typename types::element_ptr&>;

	using types::layout;
	using layout_type = Layout;

	constexpr auto layout() const -> layout_type {return array_types<T, D, ElementPtr, Layout>::layout();}

	using basic_const_array = basic_array<T, D, typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>, Layout>;

	basic_array() = default;

	constexpr basic_array(layout_type const& layout, ElementPtr const& p)
	: array_types<T, D, ElementPtr, Layout>{layout, p} {}

 protected:
	using types::types;

	template<typename, dimensionality_type, class Alloc> friend struct static_array;
	basic_array(basic_array const&) = default;

	template<class, class> friend struct basic_array_ptr;

 public:
	using typename types::element_ptr;
	using typename types::element_const_ptr;

	~basic_array() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	basic_array(basic_array&&)  // in C++ < 17 this is necessary to return references from functions
		noexcept = default;  // lints(readability-redundant-access-specifiers)

	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array&& a) {
		return std::move(a).template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array const& a) {
		return a.template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	friend constexpr auto dimensionality(basic_array const&/*self*/) {return D;}

	using typename types::reference;

	using default_allocator_type = typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type;

	constexpr auto get_allocator() const -> default_allocator_type {
		using multi::get_allocator;
		return get_allocator(this->base());
	}

	friend auto get_allocator(basic_array const& s) -> default_allocator_type {return s.get_allocator();}
	template<class P>
	static constexpr auto get_allocator_(P const& p) -> default_allocator_type {
		return multi::default_allocator_of(p);
	}

	template<class Archive>
	auto serialize(Archive& ar, unsigned int /*version*/) {
		std::for_each(this->begin(), this->end(), [&](auto&& e){ar & multi::archive_traits<Archive>::make_nvp("item", e);});
	}

	using decay_type = array<typename types::element_type, D, typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type>;

	friend constexpr auto decay(basic_array const& s) -> decay_type {return s.decay();}
	       constexpr auto decay()           const&    -> decay_type {
		decay_type ret{std::move(modify(*this))};
		return ret;
	}

	friend           auto operator+(basic_array const& s) -> decay_type {return s.decay();}

	using typename types::const_reference;

 private:
	HD constexpr auto at_(index i) const {  // MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"out of bounds");
		return reference(this->layout().sub(), this->base() + Layout::operator()(i));  // cppcheck-suppress syntaxError ; bug in cppcheck 2.5
	}

 public:
	HD constexpr auto operator[](index i) const& -> const_reference {return at_(i);}
	HD constexpr auto operator[](index i)     && ->       reference {return at_(i);}
	HD constexpr auto operator[](index i)      & ->       reference {return at_(i);}

	template<class Tuple = std::array<index, static_cast<std::size_t>(D)>, typename = std::enable_if_t<(std::tuple_size<Tuple>{} > 1)> >
	HD constexpr auto operator[](Tuple const& t) const
	->decltype(operator[](std::get<0>(t))[detail::tuple_tail(t)]) {
		return operator[](std::get<0>(t))[detail::tuple_tail(t)]; }

	template<class Tuple, typename = std::enable_if_t<std::tuple_size<Tuple>::value==1> >
	HD constexpr auto operator[](Tuple const& t) const
	->decltype(operator[](std::get<0>(t))) {
		return operator[](std::get<0>(t)); }

	template<class Tuple, std::enable_if_t<std::tuple_size<std::decay_t<Tuple>>::value==0, int> = 0>
	constexpr auto operator[](Tuple const& /*no indices*/) const -> basic_const_array {
		return *this;
	}

	using typename types::index;

	constexpr auto reindexed(typename basic_array::index first) const& -> basic_const_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	constexpr auto reindexed(typename basic_array::index first)& -> basic_array{
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	constexpr auto reindexed(typename basic_array::index first)&& -> basic_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	template<class... Indexes>
	constexpr auto reindexed(typename basic_array::index first, Indexes... idxs) const& -> basic_const_array {
		return ((reindexed(first)<<1).reindexed(idxs...))>>1;
	}
	template<class... Indexes>
	constexpr auto reindexed(typename basic_array::index first, Indexes... idxs) & -> basic_array {
		return ((reindexed(first)<<1).reindexed(idxs...))>>1;
	}
	template<class... Indexes>
	constexpr auto reindexed(typename basic_array::index first, Indexes... idxs)&& -> basic_array {
		return ((std::move(*this).reindexed(first)<<1).reindexed(idxs...))>>1;
	}

 private:
	HD constexpr auto sliced_aux(index first, index last) const -> basic_array {
		MULTI_ACCESS_ASSERT(((first==last) or this->extension().contains(first   ))&&"sliced first out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		MULTI_ACCESS_ASSERT(((first==last) or this->extension().contains(last - 1))&&"sliced last  out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		typename types::layout_t new_layout = this->layout();
		new_layout.nelems() = this->stride()*(last - first);  // TODO(correaa) : reconstruct layout instead of mutating it
		return {new_layout, types::base_ + Layout::operator()(first)};
	}

 public:
	HD constexpr auto sliced(index first, index last) const& -> basic_const_array {return sliced_aux(first, last);}
	HD constexpr auto sliced(index first, index last)      & -> basic_array       {return sliced_aux(first, last);}
	HD constexpr auto sliced(index first, index last)     && -> basic_array       {return sliced_aux(first, last);}

	constexpr auto blocked(typename basic_array::index first, typename basic_array::index last) const& -> basic_const_array {return sliced(first, last).reindexed(first);}
	constexpr auto blocked(typename basic_array::index first, typename basic_array::index last)      & -> basic_array       {return sliced(first, last).reindexed(first);}

	using iextension = typename basic_array::index_extension;

	NODISCARD("no side effects")
	constexpr auto stenciled(iextension x)                                             & -> basic_array{return blocked(x.start(), x.finish());}
	constexpr auto stenciled(iextension x, iextension x1)                              & -> basic_array{return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2)               & -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3)& -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}
	template<class... Xs>
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)& -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;}

	NODISCARD("no side effects")
	constexpr auto stenciled(iextension x)                                             && -> basic_array{return blocked(x.start(), x.finish());}
	constexpr auto stenciled(iextension x, iextension x1)                              && -> basic_array{return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2)               && -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3)&& -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}
	template<class... Xs>
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)&& -> basic_array{return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;}

	NODISCARD("no side effects")
	constexpr auto stenciled(iextension x)                                              const& -> basic_const_array {return blocked(x.start(), x.finish());}
	constexpr auto stenciled(iextension x, iextension x1)                               const& -> basic_const_array {return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2)                const& -> basic_const_array {return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3) const& -> basic_const_array {return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}

	template<class... Xs>
	constexpr auto stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)const& -> basic_const_array {
		return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;
	}

	constexpr auto elements_at(size_type n) const& -> decltype(auto) {
		assert(n < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}
	constexpr auto elements_at(size_type n) && -> decltype(auto) {
		assert(n < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}
	constexpr auto elements_at(size_type n) & -> decltype(auto) {
		assert(n < this->num_elements());
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}

	constexpr auto strided(typename types::index s) const -> basic_array{
		typename types::layout_t new_layout = *this;
		new_layout.stride_*=s;
		return {new_layout, types::base_};
	}
	constexpr auto sliced(
		typename types::index first, typename types::index last, typename types::index stride_
	) const -> basic_array {
		return sliced(first, last).strided(stride_);
	}

	using index_range = typename basic_array::index_range;

	HD constexpr auto range(index_range ir) const& -> decltype(auto) {return                  sliced(ir.front(), ir.front() + ir.size());}
	HD constexpr auto range(index_range ir)     && -> decltype(auto) {return std::move(*this).sliced(ir.front(), ir.front() + ir.size());}
	HD constexpr auto range(index_range ir)      & -> decltype(auto) {return                  sliced(ir.front(), ir.front() + ir.size());}

	constexpr auto range(typename types::index_range const& ir, dimensionality_type n) const {
		return rotated(n).range(ir).rotated(-n);
	}

	friend constexpr auto flattened(basic_array&& s) -> decltype(auto) {return std::move(s).flattened();}
	       constexpr auto flattened()&& -> decltype(auto) {
		multi::biiterator<std::decay_t<decltype(std::move(*this).begin())>> biit{std::move(*this).begin(), 0, size(*(std::move(*this).begin()))};
		return basic_array<typename std::iterator_traits<decltype(biit)>::value_type, 1, decltype(biit)>{
			multi::layout_t<1>(1, 0, this->size()*size(*(std::move(*this).begin()))),
			biit
		};
	}
	constexpr auto is_flattable() const -> bool{return this->stride() == this->layout().sub().nelems();}

	friend constexpr auto flatted(basic_array const& s) {return s.flatted();}
	       constexpr auto flatted()           const& {
		assert(is_flattable() && "flatted doesn't work for all layouts!");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		multi::layout_t<D-1> new_layout{this->layout().sub()};
		new_layout.nelems() *= this->size();  // TODO(correaa) : use immutable layout
		return basic_array<T, D-1, ElementPtr>{new_layout, types::base_};
	}

	NODISCARD("because it has no side-effect")
	constexpr auto diagonal()&& {return this->diagonal();}

	NODISCARD("because it has no side-effect")
	constexpr auto diagonal()& -> basic_array<T, D-1, typename basic_array::element_ptr> {
		auto L = std::min(std::get<0>(this->sizes()), std::get<1>(this->sizes()));
		multi::layout_t<D-1> new_layout{(*this)({0, L}, {0, L}).layout().sub()};
		new_layout.nelems() += (*this)({0, L}, {0, L}).layout().nelems();  // TODO(correaa) : don't use mutation
		new_layout.stride() += (*this)({0, L}, {0, L}).layout().stride();  // TODO(correaa) : don't use mutation
		return {new_layout, types::base_};
	}

	NODISCARD("because it has no side-effect")
	constexpr auto diagonal() const& -> basic_array<T, D-1, typename basic_array::element_const_ptr> {
		auto L = std::min(std::get<0>(this->sizes()), std::get<1>(this->sizes()));
		multi::layout_t<D-1> new_layout{(*this)({0, L}, {0, L}).layout().sub_};
		new_layout.nelems_ += (*this)({0, L}, {0, L}).layout().nelems_;
		new_layout.stride_ += (*this)({0, L}, {0, L}).layout().stride_;
		return {new_layout, types::base_};
	}

	friend constexpr auto diagonal(basic_array const& s) {return           s .diagonal();}
	friend constexpr auto diagonal(basic_array&       s) {return           s .diagonal();}
	friend constexpr auto diagonal(basic_array&&      s) {return std::move(s).diagonal();}

	using partitioned_type       = basic_array<T, D+1, element_ptr      >;
	using partitioned_const_type = basic_array<T, D+1, element_const_ptr>;

 private:
	constexpr auto partitioned_aux(size_type s) const -> partitioned_type {
		assert(s != 0);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		assert( (this->layout().nelems() % s) == 0);  // if you get an assertion here it means that you are partitioning an array with an incommunsurate partition // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : : normal in a constexpr function
		multi::layout_t<D+1> new_layout{this->layout(), this->layout().nelems()/s, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= s;
		return {new_layout, types::base_};
	}

 public:
	       constexpr auto partitioned(size_type n) const& -> partitioned_const_type {return partitioned_aux(n);}
	       constexpr auto partitioned(size_type n)      & -> partitioned_type       {return partitioned_aux(n);}
	       constexpr auto partitioned(size_type n)     && -> partitioned_type       {return partitioned_aux(n);}
	friend constexpr auto partitioned(basic_array const& s, size_type n) -> partitioned_const_type {return           s .partitioned(n);}
	friend constexpr auto partitioned(basic_array      & s, size_type n) -> partitioned_type       {return           s .partitioned(n);}
	friend constexpr auto partitioned(basic_array     && s, size_type n) -> partitioned_type       {return std::move(s).partitioned(n);}

 private:
	constexpr auto reversed_aux() const -> basic_array{
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}

 public:
	       constexpr auto reversed()           const&    -> basic_const_array {return reversed_aux();}
	       constexpr auto reversed()                &    -> basic_array       {return reversed_aux();}
	       constexpr auto reversed()               &&    -> basic_array       {return reversed_aux();}
	friend constexpr auto reversed(basic_array const& s) -> basic_const_array {return           s .reversed();}
	friend constexpr auto reversed(basic_array      & s) -> basic_array       {return           s .reversed();}
	friend constexpr auto reversed(basic_array     && s) -> basic_array       {return std::move(s).reversed();}

	constexpr auto transposed() const& -> basic_array{
		return {this->layout().transpose(), types::base_};
	}
	friend constexpr auto transposed(basic_array const& s) -> basic_array {return s.transposed();}
	friend
#if not((defined(__INTEL_COMPILER) and (__INTEL_COMPILER < 1911)) or defined(__NVCC__))
	constexpr
#endif
	auto operator~ (basic_array const& s) -> basic_array {return s.transposed();}

	HD constexpr auto rotated()      &  -> basic_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	HD constexpr auto rotated()      && -> basic_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	HD constexpr auto rotated() const& -> basic_const_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate();
		typename basic_const_array::element_ptr new_base_{types::base_};
		return basic_const_array{new_layout, new_base_};
	}

	friend constexpr auto rotated(basic_array const&  s) -> basic_const_array {return           s .rotated();}
	friend constexpr auto rotated(basic_array      &  s) -> basic_array       {return           s .rotated();}
	friend constexpr auto rotated(basic_array      && s) -> basic_array       {return std::move(s).rotated();}

	HD constexpr auto unrotated()      & {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	HD constexpr auto unrotated()     && {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	HD constexpr auto unrotated() const& {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate();
		return basic_const_array{new_layout, types::base_};
	}
	friend constexpr auto unrotated(basic_array const& self) {return self.unrotated();}

	constexpr auto rotated(dimensionality_type i) & -> basic_array{
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}
	constexpr auto rotated(dimensionality_type i) && -> basic_array{return rotated(i);}
	constexpr auto rotated(dimensionality_type i) const& -> basic_const_array{
		typename types::layout_t new_layout = this->layout();
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}

	constexpr auto unrotated(dimensionality_type i) & -> basic_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}
	constexpr auto unrotated(dimensionality_type i)     && -> basic_array       {return unrotated(i);}
	constexpr auto unrotated(dimensionality_type i) const& -> basic_const_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}

	constexpr auto operator<<(dimensionality_type i)      & -> decltype(auto) {return                    rotated(i);}
	constexpr auto operator>>(dimensionality_type i)      & -> decltype(auto) {return                  unrotated(i);}
	constexpr auto operator<<(dimensionality_type i)     && -> decltype(auto) {return std::move(*this).  rotated(i);}
	constexpr auto operator>>(dimensionality_type i)     && -> decltype(auto) {return std::move(*this).unrotated(i);}
	constexpr auto operator<<(dimensionality_type i) const& -> decltype(auto) {return                    rotated(i);}
	constexpr auto operator>>(dimensionality_type i) const& -> decltype(auto) {return                  unrotated(i);}

	constexpr auto operator|(typename basic_array::size_type n)      & -> decltype(auto) {return partitioned(n);}
	constexpr auto operator|(typename basic_array::size_type n)     && -> decltype(auto) {return std::move(*this).partitioned(n);}
	constexpr auto operator|(typename basic_array::size_type n) const& -> decltype(auto) {return partitioned(n);}

	HD constexpr auto operator()()      & -> basic_array       {return *this;}
	HD constexpr auto operator()()     && -> basic_array       {return this->operator()();}
	HD constexpr auto operator()() const& -> basic_const_array {return {this->layout(), this->base()};}

 private:
	template<typename, dimensionality_type, typename, class> friend struct basic_array;

	HD constexpr auto paren_()      & -> basic_array       {return *this;}
	HD constexpr auto paren_()     && -> basic_array       {return this->operator()();}
	HD constexpr auto paren_() const& -> basic_const_array {return {this->layout(), this->base()};}

	template<class... As> HD constexpr auto paren_(index_range a, As... as)      & {
	// return range(a).rotated().paren_(as...).unrotated();
		auto&& tmp = range(a);
		auto&& tmp2 =
			std::move(tmp).
			rotated();
		auto&& tmp3 = std::move(tmp2).paren_(as...);
		auto&& ret = std::move(tmp3).unrotated();
		return std::move(ret);
	}
	template<class... As> HD constexpr auto paren_(index_range a, As... as)     && {
		auto&& tmp = std::move(*this).range(a);
		auto&& tmp2 = std::move(tmp).rotated().paren_(as...);
		auto&& ret = std::move(tmp2).unrotated();
		return std::move(ret);
	}
	template<class... As> HD constexpr auto paren_(index_range a, As... as) const& {return range(a).rotated().paren_(as...).unrotated();}

	template<class... As> HD constexpr auto paren_(intersecting_range<index> inr, As... as)      & -> decltype(auto) {return paren_(intersection(this->extension(), inr), as...);}
	template<class... As> HD constexpr auto paren_(intersecting_range<index> inr, As... as)     && -> decltype(auto) {return paren_(intersection(this->extension(), inr), as...);}
	template<class... As> HD constexpr auto paren_(intersecting_range<index> inr, As... as) const& -> decltype(auto) {return paren_(intersection(this->extension(), inr), as...);}

	template<class... As> HD constexpr auto paren_(index i, As... as)      & -> decltype(auto) {return operator[](i).paren_(as...);}
	template<class... As> HD constexpr auto paren_(index i, As... as)     && -> decltype(auto) {return operator[](i).paren_(as...);}
	template<class... As> HD constexpr auto paren_(index i, As... as) const& -> decltype(auto) {return operator[](i).paren_(as...);}

 public:
	// vvv DO NOT remove default parameter `= irange` : the default template parameters below help interpret for {first, last} simple syntax as index ranges
	template<class B1 = irange>                                                                       HD constexpr auto operator()(B1 b1)                                const& -> decltype(auto)                          {return paren_(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr auto operator()(B1 b1, B2 b2)                         const& -> decltype(paren_(b1, b2))                {return paren_(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr auto operator()(B1 b1, B2 b2, B3 b3)                  const& -> decltype(auto)                          {return paren_(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr auto operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) const& -> decltype(paren_(b1, b2, b3, b4, as...)) {return paren_(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       HD constexpr auto operator()(B1 b1)                                     & -> decltype(auto)                          {return paren_(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr auto operator()(B1 b1, B2 b2)                              & -> decltype(paren_(b1, b2))                {return paren_(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr auto operator()(B1 b1, B2 b2, B3 b3)                       & -> decltype(auto)                          {return paren_(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr auto operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as)      & -> decltype(paren_(b1, b2, b3, b4, as...)) {return paren_(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       HD constexpr auto operator()(B1 b1)                                    && -> decltype(auto)                                                        {return std::move(*this).paren_(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr auto operator()(B1 b1, B2 b2)                             && -> decltype(std::declval<basic_array&&>().paren_(b1, b2))                {return std::move(*this).paren_(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr auto operator()(B1 b1, B2 b2, B3 b3)                      && -> decltype(auto)                                                        {return std::move(*this).paren_(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr auto operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as)     && -> decltype(std::declval<basic_array&&>().paren_(b1, b2, b3, b4, as...)) {return std::move(*this).paren_(b1, b2, b3, b4, as...);}

	template<class... Args>
	constexpr auto operator()(Args&&... args) &
	->decltype(paren(*this, std::forward<Args>(args)...)) {
		return paren(*this, std::forward<Args>(args)...); }

	template<class... Args>
	constexpr auto operator()(Args&&... args) &&
	->decltype(paren(std::move(*this), std::forward<Args>(args)...)) {
		return paren(std::move(*this), std::forward<Args>(args)...); }

	template<class... Args>
	constexpr auto operator()(Args&&... args) const&
	->decltype(paren(*this, std::forward<Args>(args)...)) {
		return paren(*this, std::forward<Args>(args)...); }

 private:
	template<typename Tuple, std::size_t ... I> constexpr auto apply_impl(Tuple const& t, std::index_sequence<I...>/*012*/) const& -> decltype(auto) {return            this->operator()(std::get<I>(t)...);}
	template<typename Tuple, std::size_t ... I> constexpr auto apply_impl(Tuple const& t, std::index_sequence<I...>/*012*/)      & -> decltype(auto) {return            this->operator()(std::get<I>(t)...);}
	template<typename Tuple, std::size_t ... I> constexpr auto apply_impl(Tuple const& t, std::index_sequence<I...>/*012*/)     && -> decltype(auto) {return std::move(*this).operator()(std::get<I>(t)...);}

 public:
	template<typename Tuple> constexpr auto apply(Tuple const& t) const& -> decltype(auto) {return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	template<typename Tuple> constexpr auto apply(Tuple const& t)     && -> decltype(auto) {return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	template<typename Tuple> constexpr auto apply(Tuple const& t)      & -> decltype(auto) {return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}

	using Layout::nelems;
	using Layout::stride;
	using Layout::sub;

	using       iterator = array_iterator<typename types::element, D, typename types::element_ptr      >;
	using const_iterator = array_iterator<typename types::element, D, typename types::element_const_ptr>;

 private:
	explicit constexpr basic_array(iterator begin, iterator end)
	: basic_array{
		layout_type{begin->layout(), begin.stride(), 0, begin.stride()*(end - begin)},
		begin.base()
	} {
		assert(begin.stride() == end.stride());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		assert(begin->layout() == end->layout());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	}
	friend auto ref<iterator>(iterator begin, iterator end) -> multi::basic_array<typename iterator::element, iterator::rank_v, typename iterator::element_ptr>;

	template<class Iterator>
	struct basic_reverse_iterator
	: std::reverse_iterator<Iterator>
	, boost::multi::totally_ordered2<basic_reverse_iterator<Iterator>, void> {
		template<class O, typename = decltype(std::reverse_iterator<Iterator>{base(std::declval<O const&>())})>
		constexpr explicit basic_reverse_iterator(O const& o) : std::reverse_iterator<Iterator>{base(o)} {}
		constexpr basic_reverse_iterator() : std::reverse_iterator<Iterator>{} {}
		constexpr explicit basic_reverse_iterator(Iterator it) : std::reverse_iterator<Iterator>(std::prev(it)) {}
		constexpr explicit operator Iterator() const {
			auto ret = this->base();
			if(ret!=Iterator{}) {return ++ret;}
			return Iterator{};
		}
		constexpr explicit operator bool() const {return static_cast<bool>(this->base());}
		constexpr auto operator==(basic_reverse_iterator const& other) const -> bool {return (this->base() == other.base());}
		constexpr auto operator*()  const -> typename Iterator::reference {return this->current;}
		constexpr auto operator->() const -> typename Iterator::pointer   {return &this->current;}
		constexpr auto operator[](typename Iterator::difference_type n) const -> typename Iterator::reference {return *(this->current - n);}
		constexpr auto operator<(basic_reverse_iterator const& o) const -> bool {return o.base() < this->base();}
	};

 public:
	using reverse_iterator = basic_reverse_iterator<iterator>;
	using ptr = basic_array_ptr<basic_array, Layout>;
	using const_ptr = basic_array_ptr<basic_const_array, Layout>;

	constexpr auto addressof() &&{return ptr{this->base_, this->layout()};}

	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&()     && {return       ptr {this->base_, this->layout()};}  // NOLINT(google-runtime-operator) // gives compiler crash in g++-7 (Ubuntu 7.5.0-6ubuntu4) 7.5.0
//  constexpr auto operator&()      & {return       ptr {this->base_, this->layout()};}
//  constexpr auto operator&() const& {return const_ptr {this->base_, this->layout()};}

	constexpr auto begin(dimensionality_type d) && -> iterator {
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(0       ), l.sub_, l.stride_};
	}
	constexpr auto end(dimensionality_type d) && -> iterator {
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(l.size()), l.sub_, l.stride_};
	}

 private:
	constexpr auto begin_aux() const {return iterator {types::base_           , sub(), stride()};}
	constexpr auto end_aux()   const {return iterator {types::base_ + nelems(), sub(), stride()};}

 public:
	       constexpr auto begin()          &    -> iterator {return begin_aux();}
	       constexpr auto end  ()          &    -> iterator {return end_aux()  ;}
	friend constexpr auto begin(basic_array& s) -> iterator {return s.begin();}
	friend constexpr auto end  (basic_array& s) -> iterator {return s.end  ();}

	       constexpr auto begin()          &&    -> iterator {return              begin();}
	       constexpr auto end  ()          &&    -> iterator {return              end()  ;}
	friend constexpr auto begin(basic_array&& s) -> iterator {return std::move(s).begin();}
	friend constexpr auto end  (basic_array&& s) -> iterator {return std::move(s).end()  ;}

	       constexpr auto begin()           const&    -> const_iterator {return begin_aux();}
	       constexpr auto end  ()           const&    -> const_iterator {return end_aux()  ;}
	friend constexpr auto begin(basic_array const& s) -> const_iterator {return s.begin();}
	friend constexpr auto end  (basic_array const& s) -> const_iterator {return s.end()  ;}

	       constexpr auto cbegin()           const& -> const_iterator {return begin();}
	       constexpr auto cend()             const& -> const_iterator {return end()  ;}
	friend constexpr auto cbegin(basic_array const& s) {return s.cbegin();}
	friend constexpr auto cend  (basic_array const& s) {return s.cend()  ;}

	template<class It> constexpr auto assign(It first) & -> It {adl_copy_n(first, this->size(), begin()); std::advance(first, this->size()); return first;}
	template<class It> constexpr auto assign(It first)&& -> It {return assign(first);}

	template<
		class Range,
		class = std::enable_if_t<not std::is_base_of<basic_array, Range>{}>,
		class = decltype(adl_copy_n(adl_begin(std::declval<Range const&>()), std::declval<typename basic_array::size_type>(), std::declval<typename basic_array::iterator>()))
	>
//  constexpr
	auto operator=(Range const& r)&  // check that you LHS is not read-only
	-> basic_array& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		assert(this->size() == r.size());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D="}+std::to_string(D)+" from range to "+typeid(T).name() );
		adl_copy_n(adl_begin(r), this->size(), begin());
		return *this;
	}

	template<class Range, class = std::enable_if_t<not std::is_base_of<basic_array, Range>{}> >
	auto operator=(Range const& r)&&
	-> basic_array&  {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		operator=(r);
		return *this;
	}

	template<class TT, class... As>
//  constexpr
	auto operator=(basic_array<TT, D, As...> const& o)& -> basic_array& {
		assert(this->extension() == o.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  MULTI_MARK_SCOPE( std::string{"multi::operator= (D="}+std::to_string(D)+") from "+typeid(TT).name()+" to "+typeid(T).name() );
		if(this->is_empty()) {return *this;}
		if(this->num_elements() == this->nelems() and o.num_elements() == this->nelems() and this->layout() == o.layout()) {
			adl_copy_n(o.base(), o.num_elements(), this->base());
		} else if(o.stride() < (~o).stride()) {
			adl_copy_n( (~o).begin(), (~o).size(), (~(*this)).begin() );
		} else {
			assign(o.begin());
		}
		return *this;
	}

	constexpr auto operator=(basic_array&& o)&&  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
	noexcept  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor) // TODO(correaa) : make conditionally noexcept
	-> basic_array& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		assert(this->extensions() == o.extensions());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		if(this->is_empty()) {return *this;}
		basic_array::operator=(o);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	template<class TT, class... As>
	constexpr auto operator=(basic_array<TT, D, As...> const& o)&&
	-> basic_array& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		operator=(o);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

//  constexpr
	auto operator=(basic_array               const& o) & -> basic_array& {
		if(this == &o) {return *this;}  // lints(cert-oop54-cpp)
		assert(this->extension() == o.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  MULTI_MARK_SCOPE("multi::operator= [D="+std::to_string(D)+"] from "+typeid(T).name()+" to "+typeid(T).name() );
		if(this->num_elements() == this->nelems() and o.num_elements() == this->nelems() and this->layout() == o.layout()) {
			adl_copy_n(o.base(), o.num_elements(), this->base());
		} else if(o.stride() < (~o).stride()) {
			adl_copy_n( (~o).begin(), (~o).size(), (~(*this)).begin() );
		} else {
			assign(o.begin());
		}
		return *this;
	}

	constexpr auto operator=(basic_array const& o) &&
	-> basic_array& {  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
		if(this == &o) {return *this;}  // lints(cert-oop54-cpp)
		operator=(o);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	template<class Array> void swap(Array&& o) && {
		assert( std::move(*this).extension() == std::forward<Array>(o).extension() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> constexpr void swap(A&& o) & {return swap(std::forward<A>(o));}

	friend constexpr void swap(basic_array&& a, basic_array&& b) {std::move(a).swap(std::move(b));}

	template<class Array> constexpr void swap(basic_array const& s, Array&& a) {s.swap(a);}
	template<class Array> constexpr void swap(Array&& a, basic_array const& s) {s.swap(a);}

	template<class Array>
	constexpr auto operator==(Array const& o) const&
	->decltype(this->extension()==o.extension() and adl_equal(this->begin(), this->end(), adl_begin(o))) {
		return this->extension()==o.extension() and adl_equal(this->begin(), this->end(), adl_begin(o)); }

	template<class Array>
	constexpr auto operator!=(Array const& o) const&
	->decltype(not (this->extension()==o.extension() and adl_equal(this->begin(), this->end(), adl_begin(o)))) {
		return not (this->extension()==o.extension() and adl_equal(this->begin(), this->end(), adl_begin(o))); }

	template<class TT, class... As>
	constexpr auto operator==(basic_array<TT, D, As...> const& o) const& -> bool {
		return (this->extension()==o.extension()) and adl_equal(this->begin(), this->end(), adl_begin(o));
	}

	template<class It>
	constexpr auto equal(It begin) const& -> bool {
		return adl_equal(
			std::move(modify(*this)).begin(),  // TODO(correaa) : what is this?
			std::move(modify(*this)).end(),
			begin
		);
	}

 private:
	friend constexpr auto lexicographical_compare(basic_array const& a1, basic_array const& a2) -> bool {
		if(a1.extension().first() > a2.extension().first()) {return true ;}
		if(a1.extension().first() < a2.extension().first()) {return false;}
		return adl_lexicographical_compare(
			a1.begin(), a1.end(),
			a2.begin(), a2.end()
		);
	}

 public:
	template<class O> constexpr auto operator<(O&& o)&& -> bool {return lexicographical_compare(std::move(*this), std::forward<O>(o));}
	template<class O> constexpr auto operator>(O&& o)&& -> bool {return lexicographical_compare(std::forward<O>(o), std::move(*this));}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() const -> basic_array<T2, D, P2> {
		P2 p2{this->base_};
		return basic_array<T2, D, P2>{this->layout(), p2};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM pm) const& -> basic_array<T2, D, P2> {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");

		return basic_array<T2, D, P2>{this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM pm) & -> basic_array<T2, D, P2> {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");

		return basic_array<T2, D, P2>{this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr auto member_cast(PM pm) && -> basic_array<T2, D, P2> {
		return this->member_cast<T2, P2, Element, PM>(pm);
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr auto reinterpret_array_cast() const& -> basic_array<std::decay_t<T2>, D, P2> {
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		auto thisbase = this->base();
		return {
			this->layout().scale(sizeof(T)/sizeof(T2)),
			static_cast<P2>(static_cast<void*>(thisbase))  // TODO(correaa) : use true reinterpret_cast
		};
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast()& -> basic_array<std::decay_t<T2>, D, P2> {
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );

		return {
			this->layout().scale(sizeof(T)/sizeof(T2)),  // NOLINT(bugprone-sizeof-expression) : sizes are compatible according to static assert above
			reinterpret_cast<P2>(this->base())  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) use customized multi/std::reinterpret_pointer_cast
		};
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast()&& ->  basic_array<std::decay_t<T2>, D, P2> {
		return this->template reinterpret_array_cast<T2, P2>();
	}

	template<class T2 = std::remove_const_t<T>, class P2 = T2*>
	constexpr auto const_array_cast()&& {
		return basic_array<std::decay_t<T2>, D, P2>{this->layout(), const_cast<P2>(this->base())};  // NOLINT(cppcoreguidelines-pro-type-const-cast) : to implement consts cast
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(multi::size_type n) & -> basic_array<std::decay_t<T2>, D + 1, P2> {
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");

		assert( n > 0 );
		assert( sizeof(T) == sizeof(T2)*static_cast<std::size_t>(n) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		typename basic_array::element_ptr const thisbase = this->base();
		P2 new_base; std::memcpy(static_cast<void*>(&new_base), static_cast<void const*>(&thisbase), sizeof(P2));  // reinterpret_cast<P2 const&>(thisbase) // TODO find a better way, fancy pointers wouldn't need reinterpret_cast
		return {
			layout_t<D+1>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}.rotate(),  // NOLINT(bugprone-sizeof-expression) T and T2 are size compatible (see static_assert above)
			new_base
		};
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(multi::size_type n)     && -> basic_array<std::decay_t<T2>, D + 1, P2>{return reinterpret_array_cast<T2, P2>(n);}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr auto reinterpret_array_cast(size_type n) const& -> basic_array<std::decay_t<T2>, D + 1, P2>{
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");

		assert( sizeof(T) == sizeof(T2)*static_cast<std::size_t>(n) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : checck implicit size compatibility
		return {
			layout_t<D+1>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}.rotate(),
			static_cast<P2>(static_cast<void*>(this->base()))
		};
	}
};

template<class Element, typename Ptr> struct array_iterator<Element, 0, Ptr>{};

template<class Element, typename Ptr>
struct array_iterator<Element, 1, Ptr>
: boost::multi::iterator_facade<
	array_iterator<Element, 1, Ptr>,
	Element, std::random_access_iterator_tag,
	typename std::iterator_traits<Ptr>::reference, multi::difference_type
>
, multi::affine          <array_iterator<Element, 1, Ptr>, multi::difference_type>
, multi::decrementable   <array_iterator<Element, 1, Ptr>>
, multi::incrementable   <array_iterator<Element, 1, Ptr>>
, multi::totally_ordered2<array_iterator<Element, 1, Ptr>, void> {
	using affine = multi::affine<array_iterator<Element, 1, Ptr>, multi::difference_type>;
	using difference_type = typename affine::difference_type;

	array_iterator() = default;

	template<
		class Other,
		decltype(_implicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().data_)* = nullptr
	>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	constexpr/*implct*/array_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_} {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to reproduce the implicitness of the argument

	template<
		class Other,
		decltype(_explicit_cast<Ptr>(typename Other::pointer{}))* = nullptr,
		decltype(std::declval<Other const&>().data_)* = nullptr
	>
	constexpr explicit array_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_} {}

	template<class, dimensionality_type, class> friend struct array_iterator;

	constexpr explicit array_iterator(std::nullptr_t nu)  : data_{nu} {}
	constexpr explicit array_iterator(Ptr const& p) : data_{p} {}

	template<
		class EElement, typename PPtr,
		typename = decltype(_implicit_cast<Ptr>(std::declval<array_iterator<EElement, 1, PPtr>>().data_))
	>
	constexpr/*implicit*/array_iterator(array_iterator<EElement, 1, PPtr> const& other)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to reproduce the implicitness of original pointer
	: data_{other.data_}, stride_{other.stride_} {}

	explicit constexpr operator bool() const {return static_cast<bool>(this->data_);}

	HD constexpr auto operator[](typename array_iterator::difference_type n) const -> typename std::iterator_traits<Ptr>::reference {
		return *((*this) + n);
	}

	constexpr auto operator->() const -> Ptr {return data_;}

	using element = Element;
	using element_ptr = Ptr;
	using pointer = element_ptr;
	using stride_type = multi::index;

	static constexpr dimensionality_type rank_v = 1;
	using rank = std::integral_constant<dimensionality_type, rank_v>;

	constexpr auto operator<(array_iterator const& o) const -> bool {return distance_to(o) > 0;}

	explicit constexpr array_iterator(Ptr d, typename basic_array<Element, 1, Ptr>::index s)
	: data_{d}, stride_{s} {}

 private:
	friend struct basic_array<Element, 1, Ptr>;

	element_ptr data_{nullptr};  // TODO(correaa) : consider uninitialized pointer
	stride_type stride_ = {1};

	constexpr auto distance_to(array_iterator const& other) const -> difference_type {
		assert(stride_==other.stride_ and (other.data_-data_)%stride_ == 0);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		return (other.data_ - data_)/stride_;
	}

 public:
	HD constexpr auto operator+(difference_type n) const -> array_iterator {array_iterator ret{*this}; ret+=n; return ret;}

	[[deprecated("use base() for iterator")]]
	constexpr auto data() const -> element_ptr {return data_;}

	constexpr auto base()              const& -> element_ptr {return data_;}

	friend  // TODO(correaa) : defined FRIEND_CONSTEXPR ?
#if not((defined(__INTEL_COMPILER) and (__INTEL_COMPILER < 1911)) or defined(__NVCC__))
	constexpr  // this generates a problem with intel compiler 19 "a constexpr function cannot have a nonliteral return type"
#endif
	auto base(array_iterator const& s) -> element_ptr {return s.base();}

	       constexpr auto stride()              const&    -> stride_type {return   stride_;}
	friend constexpr auto stride(array_iterator const& s) -> stride_type {return s.stride_;}

	constexpr auto operator++() -> array_iterator& {data_+=stride_; return *this;}
	constexpr auto operator--() -> array_iterator& {data_-=stride_; return *this;}

	friend constexpr auto operator==(array_iterator const& a, array_iterator const& b) -> bool {return    (a.data_ == b.data_);}
	friend constexpr auto operator!=(array_iterator const& a, array_iterator const& b) -> bool {return not(a.data_ == b.data_);}

	HD constexpr auto operator*() const -> typename std::iterator_traits<element_ptr>::reference {return *data_;}

	constexpr auto operator-(array_iterator const& o) const -> difference_type {return -distance_to(o);}

	constexpr auto operator+=(difference_type d) -> array_iterator& {data_+=stride_*d; return *this;}
	constexpr auto operator-=(difference_type d) -> array_iterator& {data_-=stride_*d; return *this;}
};

template<class Element, dimensionality_type D, typename... Ts>
using iterator = array_iterator<Element, D, Ts...>;

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{0}, ElementPtr, Layout>
: array_types<T, dimensionality_type(0), ElementPtr, Layout> {
	using types = array_types<T, dimensionality_type{0}, ElementPtr, Layout>;
	using types::types;

	using element      = typename types::element;
	using element_ref  = typename std::iterator_traits<typename basic_array::element_ptr>::reference;
	using element_cref = typename std::iterator_traits<typename basic_array::element_const_ptr>::reference;
	using iterator = array_iterator<T, 0, ElementPtr>;

//  constexpr
	auto operator=(element const& e) & -> basic_array& {
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D=0 from "}+typeid(T).name()+" to "+typeid(T).name() );
		adl_copy_n(&e, 1, this->base_);
		return *this;
	}
	constexpr auto operator= (element const& e) && -> basic_array& {
		operator=(e);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	constexpr auto operator==(element const& e) const& -> bool {return adl_equal(&e, std::next(&e), this->base_);}
	constexpr auto operator!=(element const& e) const& -> bool {return not operator==(e);}

	template<class TT, class=decltype(std::declval<typename basic_array::element>()==std::declval<TT>())>
	constexpr auto operator==(TT const& e) const&
	->decltype(adl_equal(&e, std::next(&e), this->base_)) {
		return adl_equal(&e, std::next(&e), this->base_); }

	template<class TT>
	constexpr auto operator!=(TT const& e) const&->decltype(!operator==(e)) {return !operator==(e);}

	template<class Range0>
	auto operator=(Range0&& r)& -> basic_array& {
		adl_copy_n(&r, 1, this->base_);
		return *this;
	}

	auto elements_at(size_type n) const& -> element_cref {assert(n < this->num_elements()); return *(this->base_);}
	auto elements_at(size_type n)     && -> element_ref  {assert(n < this->num_elements()); return *(this->base_);}
	auto elements_at(size_type n)      & -> element_ref  {assert(n < this->num_elements()); return *(this->base_);}

	constexpr auto operator!=(basic_array const& o) const& -> bool {return not adl_equal(o.base_, o.base_ + 1, this->base_);}
	constexpr auto operator==(basic_array const& o) const& -> bool {return     adl_equal(o.base_, o.base_ + 1, this->base_);}

	using decay_type = typename types::element;

	constexpr auto operator()() const& -> element_ref {return *(this->base_);}

	constexpr operator element_ref ()                            && {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	constexpr operator element_ref ()                             & {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax
	constexpr operator element_cref()                        const& {return *(this->base_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax

	template<class Archive>
	auto serialize(Archive& ar, const unsigned int /*version*/){
		ar & multi::archive_traits<Archive>::make_nvp("element", *(this->base_));
	}
};

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{1}, ElementPtr, Layout>  // NOLINT(fuchsia-multiple-inheritance) : to define operators via CRTP
: multi::partially_ordered2<basic_array<T, dimensionality_type(1), ElementPtr, Layout>, void>
, multi::random_iterable<basic_array<T, dimensionality_type(1), ElementPtr, Layout> >
, array_types<T, dimensionality_type(1), ElementPtr, Layout> {
	~basic_array() = default;  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)

	auto operator=(basic_array&& other)&
	noexcept(std::is_nothrow_copy_assignable<T>{})  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	-> basic_array& {  // lints(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
		operator=(other);
		return *this;  // lints([cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

	static constexpr dimensionality_type rank_v = 1;
	using rank = std::integral_constant<dimensionality_type, rank_v>;

	using types = array_types<T, dimensionality_type{1}, ElementPtr, Layout>;
	using types::types;

	using default_allocator_type = typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type;

	constexpr auto get_allocator()           const&       -> default_allocator_type {return default_allocator_of(basic_array::base());}
	friend    auto get_allocator(basic_array const& self) -> default_allocator_type {return self.get_allocator();}

	using decay_type = array<typename types::element, dimensionality_type{1}, typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type>;

	       constexpr auto decay()           const&       -> decay_type {return decay_type{*this};}
	friend constexpr auto decay(basic_array const& self) -> decay_type {return self.decay();}

	using basic_const_array = basic_array<
		T, 1,
		typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>,
		Layout
	>;

	using typename types::element_ptr;
	using typename types::element_const_ptr;

 protected:
	template<class A>
	constexpr void intersection_assign_(A&& other)& {
		for(auto idx : intersection(types::extension(), extension(other))) {
			operator[](idx) = std::forward<A>(other)[idx];
		}
	}
	template<class A> constexpr void intersection_assign_(A&& o)&& {intersection_assign_(std::forward<A>(o));}

	basic_array(basic_array const&) = default;

	template<class TT, dimensionality_type DD, typename EP, class LLayout> friend struct basic_array;
	template<class TT, dimensionality_type DD, class Alloc>                friend struct static_array;

	template<class T2, class P2, class TT, dimensionality_type DD, class PP>
	friend constexpr auto static_array_cast(basic_array<TT, DD, PP> const&) -> decltype(auto);

	friend auto operator+(basic_array const& self) -> decay_type {return self.decay();}

	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array&& a) {
		return std::move(a).template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array const& a) {
		return a.template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}

 public:
	basic_array(basic_array&&) noexcept = default;  // in C++ 14 this is necessary to return array references from functions
// in c++17 things changed and non-moveable non-copyable types can be returned from functions and captured by auto

 protected:
	template<class, class> friend struct basic_array_ptr;
	template<class, dimensionality_type D, class> friend struct array_iterator;

 public:
	friend constexpr auto dimensionality(basic_array const&/*self*/) -> dimensionality_type {return 1;}

	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&() && -> basic_array_ptr<basic_array, Layout> {  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed
		return {this->base_, this->layout()};
	}

	// NOLINTNEXTLINE(runtime/operator)
	constexpr auto operator&()  & -> basic_array_ptr<basic_array, Layout> {  // NOLINT(google-runtime-operator) : taking address of a reference-like object should be allowed
		return {this->base_, this->layout()};
	}

	constexpr void assign(std::initializer_list<typename basic_array::value_type> il) const {assert( il.size() == static_cast<std::size_t>(this->size()) );
		assign(il.begin(), il.end());
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

	template<class Archive>
	void serialize(Archive& ar, unsigned /*version*/) {
		std::for_each(this->begin(), this->end(), [&](auto&& e){ar& multi::archive_traits<Archive>::make_nvp("item", e);});
	}
	auto operator=(basic_array const& o)    & -> basic_array& {  // TODO(correaa) : make sfinae friendly
		if(this == &o) {return *this;}
		assert(this->extension() == o.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D=1 from "}+typeid(T).name()+" to "+typeid(T).name() );
		this->assign(o.begin(), o.end());  // TODO(correaa) : improve performance by rotating
		return *this;
	}
	template<class TT, dimensionality_type DD, class... As>
	constexpr auto operator=(basic_array const& o) && -> basic_array& {  // TODO(correaa) : make sfinae friendly
		this->operator=(o);
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}

 private:
	HD constexpr auto bracket_aux(index i) const& -> typename basic_array::reference {  // TODO(correaa): consider removing HD
		MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"out of bounds");  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		return *(this->base() + Layout::operator()(i));  // in C++17 this is allowed even with syntethic references
	}

 public:
	HD constexpr auto operator[](index i) const& -> typename basic_array::const_reference {return bracket_aux(i);}
	HD constexpr auto operator[](index i)      & -> typename basic_array::      reference {return bracket_aux(i);}
	HD constexpr auto operator[](index i)     && -> typename basic_array::      reference {return bracket_aux(i);}

 private:
	template<class Self, typename Tuple, std::size_t ... I, basic_array* = nullptr>
	friend HD constexpr auto apply_impl(Self&& self, Tuple const& t, std::index_sequence<I...> /*012*/) -> decltype(auto) {
		return std::forward<Self>(self)(std::get<I>(t)...);
	}

 public:
	template<typename Tuple> HD constexpr auto apply(Tuple const& t) const& -> decltype(auto) {return apply_impl(          *this , t, std::make_index_sequence<std::tuple_size<Tuple>{}>());}
	template<typename Tuple> HD constexpr auto apply(Tuple const& t)     && -> decltype(auto) {return apply_impl(std::move(*this), t, std::make_index_sequence<std::tuple_size<Tuple>{}>());}
	template<typename Tuple> HD constexpr auto apply(Tuple const& t)      & -> decltype(auto) {return apply_impl(          *this , t, std::make_index_sequence<std::tuple_size<Tuple>{}>());}

	template<class Tuple, std::enable_if_t< std::tuple_size<Tuple>{} == 0 , int> = 0> HD constexpr auto operator[](Tuple const& /*empty*/) const& -> decltype(auto) {return *this;}
	template<class Tuple, std::enable_if_t< std::tuple_size<Tuple>{} == 1 , int> = 0> HD constexpr auto operator[](Tuple const& indices  ) const& -> decltype(auto) {return operator[](std::get<0>(indices));}
	template<class Tuple, std::enable_if_t<(std::tuple_size<Tuple>{} >  1), int> = 0> HD constexpr auto operator[](Tuple const& indices  ) const&
	->decltype(operator[](std::get<0>(indices))[detail::tuple_tail(indices)]) {
		return operator[](std::get<0>(indices))[detail::tuple_tail(indices)]; }

	HD constexpr auto elements_at(size_type n) const& -> decltype(auto) {assert(n < this->num_elements()); return operator[](n);}
	HD constexpr auto elements_at(size_type n)     && -> decltype(auto) {assert(n < this->num_elements()); return operator[](n);}
	HD constexpr auto elements_at(size_type n)      & -> decltype(auto) {assert(n < this->num_elements()); return operator[](n);}

	using typename types::index;

	constexpr auto reindexed(typename basic_array::index first) && {return reindexed(first);}
	constexpr auto reindexed(typename basic_array::index first)  & {
		typename types::layout_t new_layout = this->layout();
		new_layout.reindex(first);
		return basic_array{new_layout, types::base_};
	}

 private:
	HD constexpr auto sliced_aux(index first, index last) const -> basic_array {
		typename types::layout_t new_layout = this->layout();
		if(this->is_empty()) {
			assert(first == last);  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
			new_layout.nelems() = 0;  // TODO(correaa) : don't use mutation
		} else {
			(new_layout.nelems() /= Layout::size())*=(last - first);
		}
		return {new_layout, types::base_ + Layout::operator()(first)};
	}

 public:
	HD constexpr auto sliced(index first, index last) const& -> basic_const_array {return sliced_aux(first, last);}
	HD constexpr auto sliced(index first, index last)      & -> basic_array       {return sliced_aux(first, last);}
	HD constexpr auto sliced(index first, index last)     && -> basic_array       {return sliced_aux(first, last);}

	constexpr auto blocked(typename basic_array::index first, typename basic_array::index last)& -> basic_array {
		return sliced(first, last).reindexed(first);
	}
	constexpr auto stenciled(typename basic_array::index_extension x) -> basic_array {
		return blocked(x.start(), x.finish());
	}

	constexpr auto strided(typename types::index s) const -> basic_array {
		typename types::layout_t new_layout = this->layout();
		new_layout.stride_*=s;
		return {new_layout, types::base_};
	}

	HD constexpr auto sliced(typename types::index first, typename types::index last, typename types::index stride) const
	-> basic_array {
		return sliced(first, last).strided(stride);
	}

	HD constexpr auto range(index_range const& ir)      & {return                  sliced(ir.front(), ir.last());}
	HD constexpr auto range(index_range const& ir)     && {return std::move(*this).sliced(ir.front(), ir.last());}
	HD constexpr auto range(index_range const& ir) const& {return                  sliced(ir.front(), ir.last());}

	constexpr auto operator()() const& -> basic_const_array {return {this->layout(), this->base()};}
	constexpr auto operator()()     && -> basic_array       {return *this;}
	constexpr auto operator()()      & -> basic_array       {return *this;}

	constexpr auto operator()(index_range const& ir)      & {return                  range(ir);}
	constexpr auto operator()(index_range const& ir)     && {return std::move(*this).range(ir);}
	constexpr auto operator()(index_range const& ir) const& {return                  range(ir);}

	HD constexpr auto operator()(index i)      & -> decltype(auto) {return                  operator[](i);}
	HD constexpr auto operator()(index i)     && -> decltype(auto) {return std::move(*this).operator[](i);}
	HD constexpr auto operator()(index i) const& -> decltype(auto) {return                  operator[](i);}

 private:
	HD constexpr auto paren_()      & {return operator()();}
	HD constexpr auto paren_()     && {return operator()();}
	HD constexpr auto paren_() const& {return operator()();}

	HD constexpr auto paren_(index_range const& ir)      & {return range(ir);}
	HD constexpr auto paren_(index_range const& ir)     && {return range(ir);}
	HD constexpr auto paren_(index_range const& ir) const& {return range(ir);}

	HD constexpr auto paren_(index i)      & -> decltype(auto) {return operator[](i);}
	HD constexpr auto paren_(index i)     && -> decltype(auto) {return operator[](i);}
	HD constexpr auto paren_(index i) const& -> decltype(auto) {return operator[](i);}

	constexpr auto paren_(intersecting_range<index> const& inr)      & -> decltype(auto) {return                  paren_(intersection(this->extension(), inr));}
	constexpr auto paren_(intersecting_range<index> const& inr)     && -> decltype(auto) {return std::move(*this).paren_(intersection(this->extension(), inr));}
	constexpr auto paren_(intersecting_range<index> const& inr) const& -> decltype(auto) {return                  paren_(intersection(this->extension(), inr));}

 public:
	constexpr auto operator()(intersecting_range<index> const& ir)      & -> decltype(auto) {return                  paren_(ir);}
	constexpr auto operator()(intersecting_range<index> const& ir)     && -> decltype(auto) {return std::move(*this).paren_(ir);}
	constexpr auto operator()(intersecting_range<index> const& ir) const& -> decltype(auto) {return                  paren_(ir);}

	template<class... Args>
	constexpr auto operator()(Args&&... args) &
	->decltype(paren(*this, std::forward<Args>(args)...)) {
		return paren(*this, std::forward<Args>(args)...); }

	template<class... Args>
	constexpr auto operator()(Args&&... args) &&
	->decltype(paren(std::move(*this), std::forward<Args>(args)...)) {
		return paren(std::move(*this), std::forward<Args>(args)...); }

	template<class... Args>
	constexpr auto operator()(Args&&... args) const&
	->decltype(paren(*this, std::forward<Args>(args)...)) {
		return paren(*this, std::forward<Args>(args)...); }

	using partitioned_type       = basic_array<T, 2, element_ptr      >;
	using partitioned_const_type = basic_array<T, 2, element_const_ptr>;

 private:
	constexpr auto partitioned_aux(size_type s) const -> partitioned_type {
		assert( s != 0 );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		assert( (this->layout().nelems() %s) == 0 );  // TODO(correaa) remove assert? truncate left over? (like mathematica) // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		multi::layout_t<2> new_layout{this->layout(), this->layout().nelems()/s, 0, this->layout().nelems()};
		new_layout.sub().nelems() /= s;  // TODO(correaa) : don't use mutation
		return {new_layout, types::base_};
	}

 public:
	constexpr auto partitioned(size_type s) const& -> partitioned_const_type {return partitioned_aux(s);}
	constexpr auto partitioned(size_type s)      & -> partitioned_type       {return partitioned_aux(s);}
	constexpr auto partitioned(size_type s)     && -> partitioned_type       {return partitioned_aux(s);}

 private:
	constexpr auto reversed_aux() const -> basic_array {
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}

 public:
	constexpr auto reversed() const& -> basic_const_array {return reversed_aux();}
	constexpr auto reversed()      & -> basic_array       {return reversed_aux();}
	constexpr auto reversed()     && -> basic_array       {return reversed_aux();}

	friend constexpr auto reversed(basic_array const& s) -> basic_const_array {return           s .reversed();}
	friend constexpr auto reversed(basic_array      & s) -> basic_array       {return           s .reversed();}
	friend constexpr auto reversed(basic_array     && s) -> basic_array       {return std::move(s).reversed();}

	friend constexpr auto   rotated(basic_array const& s) -> decltype(auto) {return s.  rotated();}
	friend constexpr auto unrotated(basic_array const& s) -> decltype(auto) {return s.unrotated();}

	constexpr auto   rotated()      & -> decltype(auto) {return operator()();}
	constexpr auto   rotated()     && -> decltype(auto) {return operator()();}
	constexpr auto   rotated() const& -> decltype(auto) {return operator()();}

	constexpr auto   rotated(dimensionality_type d)      & -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}
	constexpr auto   rotated(dimensionality_type d)     && -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}
	constexpr auto   rotated(dimensionality_type d) const& -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}

	HD constexpr auto unrotated() const& -> decltype(auto) {return operator()();}
	HD constexpr auto unrotated()     && -> decltype(auto) {return operator()();}
	HD constexpr auto unrotated()      & -> decltype(auto) {return operator()();}

	constexpr auto unrotated(dimensionality_type d)      & -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}
	constexpr auto unrotated(dimensionality_type d)     && -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}
	constexpr auto unrotated(dimensionality_type d) const& -> decltype(auto) {
		assert(d == 1); (void)d;
		return operator()();
	}

	constexpr auto operator<<(dimensionality_type i) const -> decltype(auto) {return   rotated(i);}
	constexpr auto operator>>(dimensionality_type i) const -> decltype(auto) {return unrotated(i);}

	using       iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_ptr      >;
	using const_iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_const_ptr>;
	using reverse_iterator = std::reverse_iterator<iterator>;

 private:
	constexpr       auto begin_aux() const{return iterator{this->base_                 , this->stride()};}
	constexpr       auto end_aux  () const{return iterator{this->base_ + types::nelems(), this->stride()};}

 public:
	constexpr auto begin() const& -> const_iterator {return begin_aux();}
	constexpr auto begin()      & ->       iterator {return begin_aux();}
	constexpr auto begin()     && ->       iterator {return begin_aux();}

	constexpr auto end  ()const& -> const_iterator {return end_aux();}
	constexpr auto end  ()     & ->       iterator {return end_aux();}
	constexpr auto end  ()    && ->       iterator {return end_aux();}

	friend auto begin(basic_array const& s) -> const_iterator {return           s .begin();}
	friend auto begin(basic_array      & s) ->       iterator {return           s .begin();}
	friend auto begin(basic_array     && s) ->       iterator {return std::move(s).begin();}

	friend auto end  (basic_array const& s) -> const_iterator {return           s .end()  ;}
	friend auto end  (basic_array      & s) ->       iterator {return           s .end()  ;}
	friend auto end  (basic_array     && s) ->       iterator {return std::move(s).end()  ;}

	constexpr auto cbegin() const -> const_iterator {return begin();}
	constexpr auto cend  () const -> const_iterator {return end()  ;}

	friend constexpr auto cbegin(basic_array const& s) {return s.cbegin();}
	friend constexpr auto cend  (basic_array const& s) {return s.cend()  ;}

	template<class TT, class... As, class = decltype(
		adl_copy(
			std::declval<basic_array<TT, 1, As...> const&>().begin(),
			std::declval<basic_array<TT, 1, As...> const&>().end()  ,
			std::declval<iterator>()
		)
	)>
//  constexpr
	auto operator=(basic_array<TT, 1, As...> const& other)&& -> basic_array& {
		assert( this->extensions() == other.extensions() );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D=1 from "}+typeid(TT).name()+" to "+typeid(T).name() ); // this is not the place for benchmark, benchmark implementations
		if(this->is_empty()) {return *this;}
		adl_copy(other.begin(), other.end(), this->begin());
		return *this;
	}

	template<class TT, class... As>
	auto operator=(basic_array<TT, 1, As...> const& other)& -> basic_array& {
		assert(this->extensions() == other.extensions());
		if(this->is_empty()) {return *this;}
		adl_copy(other.begin(), other.end(), this->begin());
		return *this;
	}

	template<class It> constexpr auto assign(It f)&&
	->decltype(adl_copy_n(f, this->size(), std::declval<iterator>()), void()) {
		return adl_copy_n(f, this->size(), std::move(*this).begin()), void(); }

	template<typename Array>
	constexpr auto operator==(Array const& o) const& -> bool {
		return (this->extension()==extension(o)) and adl_equal(this->begin(), this->end(), adl_begin(o));
	}

	constexpr auto operator<(basic_array const& o) const& -> bool {return lexicographical_compare(*this, o);}

	template<class Array> constexpr void swap(Array&& o)&& {
		assert(this->extension() == o.extension());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> constexpr void swap(A&& o)& {return swap(std::forward<A>(o));}

	friend constexpr void swap(basic_array&& a, basic_array&& b) {std::move(a).swap(std::move(b));}

	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend constexpr void swap(basic_array&& s, A&& a) {s.swap(a);}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend constexpr void swap(A&& a, basic_array&& s) {s.swap(a);}

 private:
	template<class A1, class A2>
	static constexpr auto lexicographical_compare(A1 const& a1, A2 const& a2) {
		if(extension(a1).first() > extension(a2).first()) {return true ;}
		if(extension(a1).first() < extension(a2).first()) {return false;}
		return adl_lexicographical_compare(adl_begin(a1), adl_end(a1), adl_begin(a2), adl_end(a2));
	}

 public:
	template<class O> constexpr auto operator<(O const& o) const -> bool {return lexicographical_compare(*this, o);}
	template<class O> constexpr auto operator>(O const& o) const -> bool {return lexicographical_compare(o, *this);}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	constexpr auto static_array_cast() const -> basic_array<T2, 1, P2> {  // name taken from std::static_pointer_cast
		return {this->layout(), static_cast<P2>(this->base())};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>, class... Args>
	constexpr auto static_array_cast(Args&&... args) const -> basic_array<T2, 1, P2> {  // name taken from std::static_pointer_cast
		return {this->layout(), P2{this->base(), std::forward<Args>(args)...}};
	}

	template<
		class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 std::decay_t<Element>::*
	> NODISCARD("")
	constexpr auto member_cast(PM pm) const -> basic_array<T2, 1, P2> {
		static_assert(sizeof(T)%sizeof(T2) == 0,
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");

#if defined(__GNUC__) and (not defined(__INTEL_COMPILER))
		auto&& r1 = (*(reinterpret_cast<typename basic_array::element_type* const&>(basic_array::base_))).*pm;  // ->*pm;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : reinterpret is what the function does. alternative for GCC/NVCC
		auto* p1 = &r1; P2 p2 = reinterpret_cast<P2&>(p1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : TODO(correaa) : find a better way
		return {this->layout().scale(sizeof(T)/sizeof(T2)), p2};
#else
		return {this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};  // this crashes nvcc 11.2-11.4 and some? gcc compiler
#endif
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	NODISCARD("")
	          auto reinterpret_array_cast() const& -> basic_array<std::decay_t<T2>, 1, P2> {
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

		static_assert( sizeof(P2) == sizeof(typename basic_array::element_ptr), "reinterpret on equal size pointers?");  // NOLINT(bugprone-sizeof-expression) : check that pointers (not pointees) are the same size
		typename basic_array::element_ptr const thisbase = this->base();
		P2 new_base; std::memcpy(static_cast<void*>(&new_base), static_cast<void const*>(&thisbase), sizeof(P2));  // NOLINT(bugprone-sizeof-expression) : take size of pointer not pointee // TODO(correaa) : find a better way, fancy pointers wouldn't need reinterpret_cast
		return {this->layout().scale(sizeof(T)/sizeof(T2)), new_base};
	}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr auto reinterpret_array_cast(size_type n) const& -> basic_array<std::decay_t<T2>, 2, P2> {
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

		auto thisbase = this->base();
		return basic_array<std::decay_t<T2>, 2, P2>{
			layout_t<2>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n},
			static_cast<P2>(static_cast<void*>(thisbase))
		}.rotated();
	}

	// TODO(correaa) : rename to reinterpret_pointer_cast?
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(size_type n)& -> basic_array<std::decay_t<T2>, 2, P2> {
		static_assert( sizeof(T)%sizeof(T2)== 0,
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases");

		typename basic_array::element_ptr const thisbase = this->base();
		return basic_array<std::decay_t<T2>, 2, P2>{
			layout_t<2>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n},
			static_cast<P2>(static_cast<void*>(thisbase))
		}.rotated();
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr auto reinterpret_array_cast(size_type n)&& -> basic_array<std::decay_t<T2>, 2, P2> {
		return this->reinterpret_array_cast(n);
	}

	template<class TT = typename basic_array::element_type>
	constexpr auto fill(TT const& value)& -> decltype(auto) {
		return adl_fill_n(this->begin(), this->size(), value), *this;
	}
	constexpr auto fill()& -> decltype(auto) {return fill(typename basic_array::element_type{});}

	template<class TT = typename basic_array::element_type>
	constexpr auto fill(TT const& value)&& -> decltype(auto) {return std::move(this->fill(value));}
	constexpr auto fill()&& -> decltype(auto) {
		return std::move(*this).fill(typename basic_array::element_type{});
	}
};

template<class T2, class P2, class Array, class... Args>
constexpr auto static_array_cast(Array&& a, Args&&... args) -> decltype(auto) {
	return a.template static_array_cast<T2, P2>(std::forward<Args>(args)...);
}

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref
// TODO(correaa) : inheredit from multi::partially_ordered2<array_ref<T, D, ElementPtr>, void>?
: basic_array<T, D, ElementPtr> {
	~array_ref() = default;  // lints(cppcoreguidelines-special-member-functions)

#if not defined(__NVCC__)  // crashes nvcc 11.3 !!!!
	constexpr auto operator=(array_ref&& other)  // lints(cppcoreguidelines-special-member-functions)
	noexcept(std::is_nothrow_copy_assignable<T>{})  // lints(hicpp-noexcept-move,performance-noexcept-move-constructor)
	-> array_ref& {
		operator=(other);  // array_refs do not rebind! elements must be copied
		return *this;  // lints(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
	}
#else
//  constexpr auto operator=(array_ref&&);
#endif

 protected:
	constexpr array_ref() noexcept
	: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{}, nullptr} {}

	using iterator = typename basic_array<T, D, ElementPtr>::iterator;

 public:  // lints(hicpp-use-equals-delete,modernize-use-equals-delete)
	array_ref(iterator, iterator) = delete;

 protected:
	[[deprecated("references are not copyable, use auto&&")]]
	array_ref(array_ref const&) = default;  // don't try to use `auto` for references, use `auto&&` or explicit value type

 public:
	array_ref(array_ref&&) noexcept = default;  // this needs to be public in c++14

	template<class OtherPtr, class=std::enable_if_t<not std::is_same<OtherPtr, ElementPtr>{}> >
	constexpr explicit array_ref(array_ref<T, D, OtherPtr>&& other)
	: basic_array<T, D, ElementPtr>{other.layout(), ElementPtr{other.base()}} {}

	constexpr explicit array_ref(typename array_ref::element_ptr p, typename array_ref::extensions_type e) noexcept  // TODO(correa) eliminate this ctor
	: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p} {}

	constexpr array_ref(typename array_ref::extensions_type e, typename array_ref::element_ptr p) noexcept
	: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p} {}

	template<class TT, std::size_t N>
	// cppcheck-suppress noExplicitConstructor ; to allow terse syntax and because a reference to c-array can be represented as an array_ref
	constexpr array_ref(  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow terse syntax and because a reference to c-array can be represented as an array_ref
		TT(&t)[N]  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	)
	: array_ref(static_cast<typename array_ref::element_ptr>(t), extensions(t)) {}

	using basic_array<T, D, ElementPtr>::operator=;
	using basic_array<T, D, ElementPtr>::operator==;

 private:
	template<class It> constexpr auto copy_elements(It first) {
		return adl_copy_n(first, array_ref::num_elements(), array_ref::data_elements());
	}
	template<class It> constexpr auto equal_elements(It first) const {
		return adl_equal(first, first + this->num_elements(), this->data_elements());
	}

 public:
	NODISCARD("")
	constexpr auto data_elements() const& -> typename array_ref::element_ptr {return array_ref::base_;}

	constexpr auto operator=(array_ref const& other) & -> array_ref& {
		if(this == &other) {return *this;}  // lints(cert-oop54-cpp)
		assert(this->num_elements() == other.num_elements());  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : normal in a constexpr function
		array_ref::copy_elements(other.data_elements());
		return *this;
	}
	constexpr auto operator=(array_ref const& other) && -> array_ref& {
		if(this == &other) {return *this;}  // lints(cert-oop54-cpp)
		operator=(other);
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
//  constexpr
	auto operator=(array_ref<TT, DD, As...> const& o)& -> array_ref& {
		assert( this->extensions() == o.extensions() );
	//  MULTI_MARK_SCOPE(std::string{"multi::operator= D="}+std::to_string(D)+" from "+typeid(TT).name()+" to "+typeid(T).name() );
		adl_copy_n(o.data_elements(), o.num_elements(), this->data_elements());
		return *this;
	}

	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr auto operator=(array_ref<TT, DD, As...> const& o)&& -> array_ref& {
		this->operator=(o);
		return *this;  // lints (cppcoreguidelines-c-copy-assignment-signature)
	}

	using  elements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_ptr      >;
	using celements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_const_ptr>;

 private:
	NODISCARD("")
	constexpr auto elements_aux() const {
		return elements_type{
			this->data_elements(),
			typename elements_type::extensions_type{multi::iextension{this->num_elements()}}
		};
	}

 public:
	NODISCARD("")
	       constexpr auto  elements()         const&       -> celements_type {return elements_aux();}
	       constexpr auto  elements()              &       ->  elements_type {return elements_aux();}
	       constexpr auto  elements()             &&       ->  elements_type {return elements_aux();}

	friend constexpr auto elements(array_ref      & s) ->  elements_type {return           s . elements();}
	friend constexpr auto elements(array_ref     && s) ->  elements_type {return std::move(s). elements();}
	friend constexpr auto elements(array_ref const& s) -> celements_type {return           s . elements();}

	NODISCARD("")
	       constexpr auto celements()         const&    {return celements_type{array_ref::data_elements(), array_ref::num_elements()};}
	friend constexpr auto celements(array_ref const& s) {return s.celements();}

	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr auto operator==(array_ref<TT, DD, As...>&& o) const& -> bool {
		if(this->extensions() != o.extensions()) {return false;}  // TODO(correaa) : or assert?
		return equal_elements(std::move(o).data_elements());
	}

	       constexpr auto data_elements()        &&    -> typename array_ref::element_ptr {return array_ref::base_;}
	friend constexpr auto data_elements(array_ref&& s) -> typename array_ref::element_ptr {return std::move(s).data_elements();}

	// data() is here for compatibility with std::vector
	template<class Dummy = void, std::enable_if_t<(D == 1) and sizeof(Dummy*), int> = 0> NODISCARD("") constexpr auto data() const& {return data_elements();}
	template<class Dummy = void, std::enable_if_t<(D == 1) and sizeof(Dummy*), int> = 0> NODISCARD("") constexpr auto data()     && {return data_elements();}
	template<class Dummy = void, std::enable_if_t<(D == 1) and sizeof(Dummy*), int> = 0> NODISCARD("") constexpr auto data()      & {return data_elements();}

	// TODO(correaa) : find a way to use [[deprecated("use data_elements()")]] for friend functions
	friend constexpr auto data(array_ref const& s) -> typename array_ref::element_ptr {return s.data_elements();}
	friend constexpr auto data(array_ref      & s) -> typename array_ref::element_ptr {return s.data_elements();}
	friend constexpr auto data(array_ref     && s) -> typename array_ref::element_ptr {return std::move(s).data_elements();}

	using decay_type = typename array_ref::decay_type;

	NODISCARD("")
	       constexpr auto decay()         const&    -> decay_type const& {return static_cast<decay_type const&>(*this);}
	friend constexpr auto decay(array_ref const& s) -> decay_type const& {return s.decay();}

	template<class Archive>
	auto serialize(Archive& ar, const unsigned int v){  // TODO(correaa) : consider small and large implementations
//  	using boost::serialization::make_nvp;
//  	if(this->num_elements() < (2<<8) )
			basic_array<T, D, ElementPtr>::serialize(ar, v);
//  	else{
//  		using boost::serialization::make_binary_object;
//  		using boost::serialization::make_array;
//  		if(std::is_trivially_copy_assignable<typename array_ref::element>{})
//  			ar & multi::archive_traits<Archive>::make_nvp("binary_data", multi::archive_traits<Archive>::make_binary_object(this->data(), sizeof(typename array_ref::element)*this->num_elements())); //#include<boost/serialization/binary_object.hpp>
//  		else ar & multi::archive_traits<Archive>::make_nvp("data", multi::archive_traits<Archive>::make_array(this->data(), this->num_elements()));
//  	}
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

template<class T, dimensionality_type D, typename Ptr = T*>
struct array_ptr
: basic_array_ptr<basic_array<T, D, Ptr>
, typename array_ref<T, D, Ptr>::layout_t> {
	using basic_ptr = basic_array_ptr<basic_array<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;

	constexpr array_ptr(Ptr p, multi::extensions_t<D> x) : basic_ptr{p, multi::layout_t<D>{x}} {}

	constexpr explicit array_ptr(std::nullptr_t p) : array_ptr{p, multi::extensions_t<D>{}} {}

	template<typename CArray>
	constexpr explicit array_ptr(CArray* p) : array_ptr{data_elements(*p), extensions(*p)} {}

	constexpr auto operator*() const {
		return array_ref<T, D, Ptr>{this->base(), (*this)->extensions()};
	}
};

template<class T, typename Ptr>
class array_ptr<T, 0, Ptr> : multi::array_ref<T, 0, Ptr>{
 public:
	constexpr explicit array_ptr(Ptr p, typename multi::array_ref<T, 0, Ptr>::extensions_type x) : multi::array_ref<T, 0, Ptr>(p, x) {}
	constexpr explicit array_ptr(Ptr p) : array_ptr(p, typename multi::array_ref<T, 0, Ptr>::extensions_type{}) {}

	constexpr explicit operator bool() const {return this->base();}
	constexpr explicit operator Ptr () const {return this->base();}

	friend constexpr auto operator==(array_ptr const& self, array_ptr const& other) -> bool {return self.base() == other.base();}
	friend constexpr auto operator!=(array_ptr const& self, array_ptr const& other) -> bool {return self.base() != other.base();}

	constexpr auto operator* () const -> multi::array_ref<T, 0, Ptr>& {return const_cast<array_ptr&>(*this);}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) : find a way to avoid using const_cast
	constexpr auto operator->() const -> multi::array_ref<T, 0, Ptr>* {return const_cast<array_ptr*>( this);}  // NOLINT(cppcoreguidelines-pro-type-const-cast) : TODO(correaa) : find a way to avoid using const_cast
};

template<class TT, std::size_t N>
constexpr auto addressof(TT(&t)[N]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return array_ptr<
		std::decay_t<std::remove_all_extents_t<TT[N]>>, std::rank<TT[N]>{}, std::remove_all_extents_t<TT[N]>*  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	>(&t);
}

template<class T, dimensionality_type D, typename Ptr = T*>
using array_cptr = array_ptr<T, D, 	typename std::pointer_traits<Ptr>::template rebind<T const>>;

template<dimensionality_type D, class P>
constexpr auto make_array_ref(P p, multi::extensions_t<D> x) {
	return array_ref<typename std::iterator_traits<P>::value_type, D, P>(p, x);
}

template<class P> auto make_array_ref(P p, extensions_t<0> x) {return make_array_ref<0>(p, x);}
template<class P> auto make_array_ref(P p, extensions_t<1> x) {return make_array_ref<1>(p, x);}
template<class P> auto make_array_ref(P p, extensions_t<2> x) {return make_array_ref<2>(p, x);}
template<class P> auto make_array_ref(P p, extensions_t<3> x) {return make_array_ref<3>(p, x);}
template<class P> auto make_array_ref(P p, extensions_t<4> x) {return make_array_ref<4>(p, x);}
template<class P> auto make_array_ref(P p, extensions_t<5> x) {return make_array_ref<5>(p, x);}

// In ICC you need to specify the dimensionality in make_array_ref<D>
// #if defined(__INTEL_COMPILER)
// template<dimensionality_type D, class P>
// auto make_array_ref(P p, std::initializer_list<index_extension> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
// template<dimensionality_type D, class P>
// auto make_array_ref(P p, std::initializer_list<index> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
// #endif

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
	typename V = typename std::remove_all_extents<T[N]>::type, std::size_t D = std::rank<T[N]>{}   // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
>
array_ptr(T(*)[N])->array_ptr<V, static_cast<multi::dimensionality_type>(D)>;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility

template<class Ptr> array_ref(Ptr, index_extensions<0>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 0, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<1>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 1, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<2>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 2, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<3>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 3, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<4>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 4, Ptr>;
template<class Ptr> array_ref(Ptr, index_extensions<5>)->array_ref<typename std::iterator_traits<Ptr>::value_type, 5, Ptr>;

template<class It, class Tuple> array_ref(It, Tuple)->array_ref<typename std::iterator_traits<It>::value_type, std::tuple_size<Tuple>::value, It>;
#endif

// TODO(correaa) move to utility
template<class T, std::size_t N>
constexpr auto rotated(const T(&t)[N]) noexcept {                                                  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(t))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		base(t), extensions(t)
	).rotated();
}
template<class T, std::size_t N>
constexpr auto rotated(T(&t)[N]) noexcept {                                                        // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(t))>(  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : backwards compatibility
		base(t), extensions(t)
	).rotated();
}

template<class TD, class Ptr> struct Array_aux;
template<class T, std::size_t D, class Ptr> struct Array_aux<   T[D], Ptr>{using type = array    <T, static_cast<multi::dimensionality_type>(D), Ptr>  ;};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : enable syntax
template<class T, std::size_t D, class Ptr> struct Array_aux<T(&)[D], Ptr>{using type = array_ref<T, static_cast<multi::dimensionality_type>(D), Ptr>&&;};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : enable syntax
template<class T, std::size_t D, class Ptr> struct Array_aux<T(*)[D], Ptr>{using type = array_ptr<T, static_cast<multi::dimensionality_type>(D), Ptr>  ;};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : enable syntax

template<class TD, class Second =
	std::conditional_t<
		std::is_reference<TD>{} or std::is_pointer<TD>{},
		std::add_pointer_t<std::remove_all_extents_t<std::remove_reference_t<std::remove_pointer_t<TD>>>>,
		std::allocator<std::remove_all_extents_t<TD>>
	>
> using Array = typename Array_aux<TD, Second>::type;

template<class RandomAccessIterator, dimensionality_type D>
constexpr auto operator/(RandomAccessIterator data, multi::extensions_t<D> x)
-> multi::array_ptr<typename std::iterator_traits<RandomAccessIterator>::value_type, D, RandomAccessIterator>
{return {data, x};}

template<class T, dimensionality_type D, class... Ts>
constexpr auto is_basic_array_aux(basic_array<T, D, Ts...> const&) -> std::true_type;
constexpr auto is_basic_array_aux(...                            ) -> std::false_type;

template<class A> struct is_basic_array: decltype(is_basic_array_aux(std::declval<A>())){};

template<class In, class T, dimensionality_type N, class TP, class = std::enable_if_t<(N > 1)>, class = decltype(adl_begin(*In{}), adl_end(*In{}))>
constexpr auto uninitialized_copy
// require N>1 (this is important because it forces calling placement new on the pointer
(In first, In last, multi::array_iterator<T, N, TP> dest) {
	using std::begin; using std::end;
	while(first != last) {
		adl_uninitialized_copy(adl_begin(*first), adl_end(*first), adl_begin(*dest));
		++first;
		++dest;
	}
	return dest;
}

// begin and end for forwarding reference are needed in this namespace
// to overwrite the behavior of std::begin and std::end
// which take rvalue-references as const-references.

template<class T> auto begin(T&& t)
->decltype(std::forward<T>(t).begin()) {
	return std::forward<T>(t).begin(); }

template<class T> auto end(T&& t)
->decltype(std::forward<T>(t).end()) {
	return std::forward<T>(t).end(); }

}  // end namespace multi
}  // end namespace boost
#endif

