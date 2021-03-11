#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0x&&$0x&&rm $0x&&(rm -rf test/build&&mkdir -p test/build&&cd test/build&&time cmake ..&&make -j 8&&time ctest --output-on-failure -j 8);exit
#endif
// © Alfredo Correa 2018-2020

//#if (defined(__clang__) and defined(__CUDA__)) or defined(__NVCC__)
//#define BOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK // see comments https://www.boost.org/doc/libs/1_72_0/boost/utility/result_of.hpp
//#endif

#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP

#include "./memory/pointer_traits.hpp"
#include "utility.hpp" 

#include "./detail/layout.hpp"
#include "./detail/types.hpp"     // dimensionality_type
#include "./detail/operators.hpp" // random_iterable
#include "./detail/memory.hpp"    // pointer_traits

#include "./config/NODISCARD.hpp"
#include "./config/DELETE.hpp"
#include "./config/ASSERT.hpp"
#include "./config/MARK.hpp"

#if defined(__NVCC__)
#define HD __host__ __device__
#else
#define HD
#endif

//#include<iostream> // debug

#include<algorithm>  // copy_n
#include<cstring>    // for memset in reinterpret_cast
#include<functional> // invoke
#include<memory>     // pointer_traits

namespace std{
	template<class T>
	struct pointer_traits<std::move_iterator<T*>> : std::pointer_traits<T*>{
		template<class U> using rebind = 
			std::conditional_t<std::is_const<U>{}, 
				U*,
				std::pointer_traits<std::move_iterator<U*>>
			>;
	};
}

namespace boost{
namespace multi{

template<class T> T& modify(T const& t){return const_cast<T&>(t);}

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>> 
struct basic_array;

template<typename T, dimensionality_type D, class A = std::allocator<T>> struct array;

template<class To, class From, std::enable_if_t<std::is_convertible<From, To>{},int> =0>
constexpr To _implicit_cast(From&& f){return static_cast<To>(f);}

template<class To, class From, std::enable_if_t<std::is_constructible<To, From>{} and not std::is_convertible<From, To>{},int> =0>
constexpr To _explicit_cast(From&& f){return static_cast<To>(f);}

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>>
struct array_types : Layout{
	using element = T;
	using element_type = element; // this follows more closely https://en.cppreference.com/w/cpp/memory/pointer_traits
	constexpr static dimensionality_type dimensionality = D;
	using element_ptr = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element_type const>; //multi::const_iterator<ElementPtr>; 
	using element_ref = typename std::iterator_traits<element_ptr>::reference;
	using layout_t = Layout;
	using value_type = typename std::conditional<
		(dimensionality>1),
		array<element, dimensionality-1, typename multi::pointer_traits<element_ptr>::default_allocator_type>, 
		typename std::conditional<
			dimensionality == 1,
			element,
			element
		>::type
	>::type;
	using reference = typename std::conditional<(dimensionality > 1), 
		basic_array<element, dimensionality-1, element_ptr>,
		typename std::conditional<(dimensionality == 1), 
			  typename std::iterator_traits<element_ptr>::reference   // this seems more correct but it doesn't work with cuda fancy reference
			, typename std::iterator_traits<element_ptr>::reference
		>::type
	>::type;

	using const_reference = typename std::conditional<(dimensionality > 1), 
		basic_array<element, dimensionality-1, element_const_ptr>,
	//	decltype(*std::declval<element_const_ptr>())&
		typename std::iterator_traits<element_const_ptr>::reference
	//	typename std::pointer_traits<element_const_ptr>::reference   // this seems more correct but it doesn't work with cuda fancy reference
	>::type;

	HD constexpr element_ptr        base() const{return base_;}
       constexpr element_const_ptr cbase() const{return base_;}

	constexpr element_ptr& mbase() const{return base_;}
	friend constexpr element_ptr base(array_types const& s){return s.base();}
	constexpr layout_t const& layout() const{return *this;}
	friend constexpr layout_t const& layout(array_types const& s){return s.layout();}
	constexpr element_ptr            origin() const{return base_+Layout::origin();} //	element_const_ptr     corigin() const{return origin();}
	friend constexpr decltype(auto)  origin(array_types const& s){return s.origin();} //	friend decltype(auto) corigin(array_types const& s){return s.corigin();}
protected:
	using derived = basic_array<T, D, ElementPtr, Layout>;
	element_ptr base_;
	constexpr explicit array_types(std::nullptr_t np) : Layout{}, base_{np}{}
public:
	array_types() = default;
//#if defined(__NVCC__) 
//	__host__ __device__ // TODO check why this is necessary (nvcc 11), removing this gives a, trivial_device_copy D->H failed: cudaErrorLaunchFailure: unspecified launch failure
//#endif
	constexpr array_types(layout_t const& l, element_ptr const& data): Layout{l}, base_{data}{}
	array_types(array_types const&) = default;
//	template<class T2, class P2, class Array> friend decltype(auto) static_array_cast(Array&&);
public://TODO find why this needs to be public and not protected or friend
	template<class ArrayTypes, typename = std::enable_if_t<not std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>
		, decltype(_implicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	> // cppcheck-suppress noExplicitConstructor ; because underlying pointers are implicitly convertible
	constexpr array_types(ArrayTypes const& a) : Layout{a}, base_{a.base_}{}
	template<class ArrayTypes, typename = std::enable_if_t<not std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>
		, decltype(_explicit_cast<element_ptr>(std::declval<ArrayTypes const&>().base_))* = nullptr
	>
	constexpr explicit array_types(ArrayTypes const& a) : Layout{a}, base_{a.base_}{}



	template<typename ElementPtr2, 
		typename = decltype(Layout{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().layout()}),
		typename = decltype(element_ptr{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().base_})
	>
	constexpr array_types(array_types<T, D, ElementPtr2, Layout> const& other) : Layout{other.layout()}, base_{other.base_}{}
	template<class T2, dimensionality_type D2, class E2, class L2> friend struct array_types;
};

//template<class T, dimensionality_type D, class ElementPtr = T*>
//struct BasicArrayPtr{
//	using element_ptr = ElementPtr;
//	static constexpr dimensionality_type dimensionality = D;
//	using reference = basic_array<T, D, ElementPtr>;
//private:
//	using layout_type = layout_t<dimensionality>;
//	element_ptr base_;
//	layout_type layout_;
//	constexpr BasicArrayPtr(element_ptr base, layout_type layout) : base_{base}, layout_{layout}{}
//	template<class, dimensionality_type, class, class> friend struct basic_array;
//public:
//	BasicArrayPtr() = default;
//	constexpr BasicArrayPtr(BasicArrayPtr const& o) : base_{o.base_}, layout_{o.layout_}{}
//	constexpr BasicArrayPtr(std::nullptr_t p) : base_{p}{}
//	constexpr BasicArrayPtr& operator=(BasicArrayPtr const&) = default;
//	constexpr bool operator==(BasicArrayPtr const& o) const{return base_==o.base_ and layout_==o.layout_;}
//	constexpr bool operator!=(BasicArrayPtr const& o) const{return base_!=o.base_ or  layout_!=o.layout_;}
//	constexpr explicit operator bool() const{return base_;}
//	constexpr reference Dereference() const{return {layout_, base_};}
//};

template<class Ref, class Layout>
struct basic_array_ptr : 
	private Ref,
	boost::multi::iterator_facade<
		basic_array_ptr<Ref, Layout>, void, std::random_access_iterator_tag, 
		Ref const&, typename Layout::difference_type
	>//, boost::multi::totally_ordered2<basic_array_ptr<Ref, Layout>, void>
{
	using pointer = Ref const*;
	using element_type = typename Ref::decay_type;
	using difference_type = typename Layout::difference_type;

	using value_type = element_type;
	using reference = Ref;// const&;
	using iterator_category = std::random_access_iterator_tag;

	constexpr basic_array_ptr(std::nullptr_t p = nullptr) : Ref{p}{} // TODO remove default argument, add default ctor
	template<class, class> friend struct basic_array_ptr;
	constexpr basic_array_ptr(typename Ref::element_ptr p, layout_t<Ref::dimensionality-1> l) : Ref{l, p}{}
	constexpr basic_array_ptr(typename Ref::element_ptr p, index_extensions<Ref::dimensionality> e) : Ref{p, e}{}

//	template<class Array, typename = decltype(typename Ref::element_ptr{typename Array::element_ptr{}})> 
//	constexpr basic_array_ptr(Array const& o) : Ref{o->layout(), o->base()}{}//, stride_{o.stride_}{}
//	constexpr basic_array_ptr(basic_array_ptr const& o) : Ref{static_cast<Layout const&>(o), o.base_}{}//, stride_{o.stride_}{}
	basic_array_ptr(basic_array_ptr&& o) = default;//: Ref{static_cast<Layout const&>(o), o.base_}{}//, stride_{o.stride_}{}
	basic_array_ptr(basic_array_ptr const& o) = default;//: Ref{static_cast<Layout const&>(o), o.base_}{}//, stride_{o.stride_}{}
	basic_array_ptr& operator=(basic_array_ptr const& other){
		this->base_ = other.base_;
		static_cast<Layout&>(*this) = other;
		return *this;
	}
	constexpr explicit operator bool() const{return this->base_;}
	constexpr Ref  dereference() const{return Ref{this->layout(), this->base_};}
	HD constexpr Ref  operator* () const{return Ref{*this};}
	constexpr Ref* operator->() const{return  const_cast<basic_array_ptr*>(this);}
	constexpr Ref* operator->(){return  this;}
	constexpr Ref  operator[](difference_type n) const{return *(*this + n);}
//	template<class O> bool operator==(O const& o) const{return equal(o);}
	constexpr bool operator<(basic_array_ptr const& o) const{return distance_to(o) > 0;}
	constexpr basic_array_ptr(typename Ref::element_ptr p, Layout const& l) : Ref{l, p}{}
	template<typename T, dimensionality_type D, typename ElementPtr, class LLayout>
	friend struct basic_array;
	constexpr auto base() const{return this->base_;}
	friend constexpr auto base(basic_array_ptr const& self){return self.base();}
	using Ref::base_;
	using Ref::layout;
	constexpr bool operator==(basic_array_ptr const& o) const{return base_==o.base_ and layout()==o.layout();}
	template<class O> constexpr bool operator==(O const& o) const{return base()==o->base() and layout() == o->layout();}
	template<class O> constexpr bool operator!=(O const& o) const{return not ((*this)==o);}
	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr bool operator==(O const& o, basic_array_ptr const& s){return s==o;}
	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr bool operator!=(O const& o, basic_array_ptr const& s){return not(o==s);}
protected:
	constexpr void increment(){base_ += Ref::nelems();}
	constexpr void decrement(){base_ -= Ref::nelems();}
	constexpr void advance(difference_type n){base_ += Ref::nelems()*n;}
	constexpr difference_type distance_to(basic_array_ptr const& other) const{
		assert( Ref::nelems() == other.Ref::nelems() and Ref::nelems() != 0 );
		assert( (other.base_ - base_)%Ref::nelems() == 0); 
		assert( layout() == other.layout() );
		return (other.base_ - base_)/Ref::nelems();
	}
public:
	constexpr basic_array_ptr& operator+=(difference_type n){advance(n); return *this;}
};

template<class Element, dimensionality_type D, typename ElementPtr>
struct array_iterator;

template<class Element, dimensionality_type D, typename ElementPtr>
struct array_iterator : 
	boost::multi::iterator_facade<
		array_iterator<Element, D, ElementPtr>, void, std::random_access_iterator_tag, 
		basic_array<Element, D-1, ElementPtr> const&, typename layout_t<D-1>::difference_type
	>,
	multi::decrementable<array_iterator<Element, D, ElementPtr>>,
	multi::incrementable<array_iterator<Element, D, ElementPtr>>,
	multi::affine<array_iterator<Element, D, ElementPtr>, multi::difference_type>,
	multi::totally_ordered2<array_iterator<Element, D, ElementPtr>, void>
{
	using difference_type = typename layout_t<D>::difference_type;
	using element = Element;//typename Ref::element;
	using element_ptr = ElementPtr;//typename Ref::element_ptr;
	using value_type = typename basic_array<element, D-1, element_ptr>::decay_type;

	using pointer   = basic_array<element, D-1, element_ptr>*;
	using reference = basic_array<element, D-1, element_ptr>&&;//Ref const&;
//	using element_type = typename Ref::value_type;
	using iterator_category = std::random_access_iterator_tag;

	using rank = std::integral_constant<dimensionality_type, D>;
	
	using ptr_type = basic_array_ptr<basic_array<element, D-1, element_ptr>, layout_t<D-1>>;
	using stride_type = index;

	constexpr array_iterator(std::nullptr_t p = nullptr) : ptr_{p}, stride_{1}{}//Ref{p}{}
	template<class, dimensionality_type, class> friend struct array_iterator;
//	template<class Other, typename = decltype(typename basic_array<element, D-1, element_ptr>::types::element_ptr{typename Other::element_ptr{}})> 
//	constexpr array_iterator(Other const& o) : /*Ref{o.layout(), o.base()},*/ ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_}{}

	template<class EElement, typename PPtr, 
		decltype(_implicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().ptr_.base()))* = nullptr // .base() (instead of .base_) is needed due to a bug in nvcc 11.1 not seeing the friend declaration?
	>
	constexpr          array_iterator(array_iterator<EElement, D, PPtr> const& o) : ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_}{} // TODO refactor basic_array_ptr to not depend on Ref template parameter
	template<class EElement, typename PPtr, 
		decltype(_explicit_cast<ElementPtr>(std::declval<array_iterator<EElement, D, PPtr>>().ptr_.base()))* = nullptr
	>
	constexpr explicit array_iterator(array_iterator<EElement, D, PPtr> const& o) : ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_}{} 

	array_iterator(array_iterator const&) = default;
	array_iterator& operator=(array_iterator const& other) = default;

	explicit constexpr operator bool() const{return static_cast<bool>(ptr_.base_);}
	HD constexpr basic_array<element, D-1, element_ptr> operator*() const{/*assert(*this);*/ return {*ptr_};}//return *this;}
	constexpr decltype(auto) operator->() const{/*assert(*this);*/ return ptr_;}//return this;}
	HD constexpr array_iterator operator+(difference_type n) const{array_iterator ret{*this}; ret+=n; return ret;}
	HD constexpr basic_array<element, D-1, element_ptr> operator[](difference_type n) const{return *((*this) + n);}

	constexpr bool operator==(array_iterator const& o) const{return ptr_==o.ptr_ and stride_==o.stride_ and ptr_.layout() == o.ptr_.layout();}
//	template<class O> constexpr bool operator==(O const& o) const{return equal(o);}
	constexpr bool operator<(array_iterator const& o) const{return distance_to(o) > 0;}
	constexpr array_iterator(typename basic_array<element, D-1, element_ptr>::element_ptr p, layout_t<D-1> l, index stride) : /*Ref{l, p},*/
		ptr_{p, l}, 
		stride_{stride}
	{}
	template<class, dimensionality_type, class, class> friend struct basic_array;
	template<class... As> HD constexpr decltype(auto) operator()(index i, As... as) const{
		return this->operator[](i)(as...);
	}
	                      HD constexpr decltype(auto) operator()(index i          ) const{return this->operator[](i)       ;}

private:
	template<typename Tuple, std::size_t ... I> 
	HD constexpr decltype(auto) apply_impl(Tuple const& t, std::index_sequence<I...>) const{return this->operator()(std::get<I>(t)...);}
public:
	template<typename Tuple> HD constexpr decltype(auto) apply(Tuple const& t) const{return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
private:
	ptr_type ptr_;
	stride_type stride_ = {1}; // nice non-zero default
	constexpr bool equal(array_iterator const& o) const{return ptr_==o.ptr_ and stride_==o.stride_;}//base_==o.base_ && stride_==o.stride_ && ptr_.layout()==o.ptr_.layout();}
	constexpr void decrement(){ptr_.base_ -= stride_;}
	constexpr void advance(difference_type n){ptr_.base_ += stride_*n;}
	constexpr difference_type distance_to(array_iterator const& other) const{
		assert( stride_ == other.stride_); assert( stride_ != 0 );
	//	assert( this->stride()==stride(other) and this->stride() );// and (base(other.ptr_) - base(this->ptr_))%stride_ == 0
	//	assert( stride_ == other.stride_ and stride_ != 0 and (other.ptr_.base_-ptr_.base_)%stride_ == 0 and ptr_.layout() == other.ptr_.layout() );
	//	assert( stride_ == other.stride_ and stride_ != 0 and (other.base_ - base_)%stride_ == 0 and layout() == other.layout() );
		return (other.ptr_.base_ - ptr_.base_)/stride_;
	}
public:
	       constexpr element_ptr base()              const&   {return ptr_.base_;}
	friend constexpr element_ptr base(array_iterator const& s){return s.base();}
	       constexpr stride_type stride()              const&   {return   stride_;}
	friend constexpr stride_type stride(array_iterator const& s){return s.stride_;}
	constexpr array_iterator& operator++(){ptr_.base_ += stride_; return *this;}
	constexpr array_iterator& operator--(){decrement(); return *this;}

	friend constexpr difference_type operator-(array_iterator const& self, array_iterator const& other){
		assert(self.stride_ == other.stride_); assert(self.stride_ != 0);
		return (self.ptr_.base_ - other.ptr_.base_)/self.stride_;
	}
	constexpr array_iterator& operator+=(difference_type d){advance(+d); return *this;}
	constexpr array_iterator& operator-=(difference_type d){advance(-d); return *this;}
};

template<class It>
struct biiterator : 
	boost::multi::iterator_facade<
		biiterator<It>,
		typename std::iterator_traits<It>::value_type, std::random_access_iterator_tag, 
		decltype(*(std::move((*std::declval<It>())).begin())), multi::difference_type
	>,
	multi::affine<biiterator<It>, multi::difference_type>,
	multi::decrementable<biiterator<It>>,
	multi::incrementable<biiterator<It>>,
	multi::totally_ordered2<biiterator<It>, void>
{
	It me_;
	std::ptrdiff_t pos_;
	std::ptrdiff_t stride_;
	biiterator() = default;
	biiterator(biiterator const& other) = default;// : me{other.me}, pos{other.pos}, stride{other.stride}{}
	constexpr biiterator(It me, std::ptrdiff_t pos, std::ptrdiff_t stride) : me_{me}, pos_{pos}, stride_{stride}{}
	constexpr decltype(auto) operator++(){
		++pos_;
		if(pos_==stride_){
			++me_;
			pos_ = 0;
		}
		return *this;
	}
	constexpr bool operator==(biiterator const& o) const{return me_==o.me_ and pos_==o.pos_;}
	constexpr biiterator& operator+=(multi::difference_type n){me_ += n/stride_; pos_ += n%stride_; return *this;}
	constexpr decltype(auto) operator*() const{
		auto meb = std::move(*me_).begin();
		return meb[pos_];
	}
	using difference_type = std::ptrdiff_t;
	using reference = decltype(*std::declval<biiterator>());
	using value_type = std::decay_t<reference>;
	using pointer = value_type*;
	using iterator_category = std::random_access_iterator_tag;
};

template<typename T, dimensionality_type D, typename ElementPtr, class Layout>
struct basic_array : 
	multi::partially_ordered2<basic_array<T, D, ElementPtr, Layout>, void>,
	array_types<T, D, ElementPtr, Layout>
{
	using types = array_types<T, D, ElementPtr, Layout>;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr >;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr&>;
	using types::layout;
	using layout_type = Layout;
	constexpr layout_type layout() const{return array_types<T, D, ElementPtr, Layout>::layout();}
	using basic_const_array = basic_array<T, D, typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>, Layout>;
	basic_array() = default;
	constexpr basic_array(layout_type const& layout, ElementPtr const& p) : array_types<T, D, ElementPtr, Layout>{layout, p}{}
protected:
	using types::types;
	template<typename, dimensionality_type, class Alloc> friend struct static_array;
	basic_array(basic_array const&) = default;
	template<class, class> friend struct basic_array_ptr;
public:
	using typename types::element_ptr;
	using typename types::element_const_ptr;
//#if __cplusplus >= 201703L
//#if defined(__INTEL_COMPILER) or defined(__NVCC__)
//public:
//#else
//protected:
//#endif
//	constexpr basic_array(basic_array&&) = default; // if you need to generate a copy you can't use `auto` here, use `decay`. maybe you want to return `decltype(auto)`.
//#else
//public: 
//	constexpr basic_array(basic_array&&) = default; // in C++ < 17 this is necessary to return references from functions
//#endif
public:
	basic_array(basic_array&&) = default; // in C++ < 17 this is necessary to return references from functions
public:
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array&& a){
		return std::move(a).template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array const& a){
		return a.template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	friend constexpr auto dimensionality(basic_array const& self){return self.dimensionality;}
	using typename types::reference;

	using default_allocator_type = typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type;

	constexpr default_allocator_type get_allocator() const{
		using multi::get_allocator;
		return get_allocator(this->base());
	}
	
	friend constexpr default_allocator_type get_allocator(basic_array const& s){return s.get_allocator();}
	template<class P>
	static constexpr default_allocator_type get_allocator_(P const& p){
		return multi::default_allocator_of(p);
	}
	using decay_type = array<typename types::element_type, D, typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type>;//get_allocator_(std::declval<ElementPtr>()))>;
	template<class Archive>
	auto serialize(Archive& ar, unsigned int /*file version*/){
		std::for_each(this->begin(), this->end(), [&](auto&& e){ar & multi::archive_traits<Archive>::make_nvp("item", e);});
	}
	constexpr decay_type decay() const{
		decay_type ret = std::move(modify(*this));
		return ret;
	}
	friend constexpr decay_type decay(basic_array const& s){return s.decay();}
	friend decay_type operator+(basic_array const& s){return s.decay();}

	using typename types::const_reference;

private:
	HD constexpr auto at_(index i) const{//MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"out of bounds");
		return reference(this->layout().sub_, this->base() + Layout::operator()(i));
	}
public:
	HD constexpr const_reference operator[](index i) const&{return at_(i);}
	HD constexpr       reference operator[](index i)     &&{return at_(i);}
	HD constexpr       reference operator[](index i)      &{return at_(i);}

	template<class Tp = std::array<index, static_cast<std::size_t>(D)>, typename = std::enable_if_t<(std::tuple_size<std::decay_t<Tp>>{}>1)> >
	HD constexpr auto operator[](Tp&& t) const
	->decltype(operator[](std::get<0>(t))[detail::tuple_tail(t)]){
		return operator[](std::get<0>(t))[detail::tuple_tail(t)];}
	template<class Tp, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tp>>::value==1> >
	HD constexpr auto operator[](Tp&& t) const
	->decltype(operator[](std::get<0>(t))){
		return operator[](std::get<0>(t));}
	template<class Tp = std::tuple<>, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tp>>::value==0> >
	HD constexpr decltype(auto) operator[](Tp&&) const{return *this;}
	using typename types::index;
	constexpr basic_const_array reindexed(typename basic_array::index first) const&{
		typename types::layout_t new_layout = *this;
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	constexpr basic_array reindexed(typename basic_array::index first)&{
		typename types::layout_t new_layout = *this;
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	constexpr basic_array reindexed(typename basic_array::index first)&&{
		typename types::layout_t new_layout = *this;
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	template<class... Indexes>
	constexpr basic_const_array reindexed(typename basic_array::index first, Indexes... idxs) const&{
		return ((reindexed(first)<<1).reindexed(idxs...))>>1;
	}
	template<class... Indexes>
	constexpr basic_array reindexed(typename basic_array::index first, Indexes... idxs) &{
		return ((reindexed(first)<<1).reindexed(idxs...))>>1;
	}
	template<class... Indexes>
	constexpr basic_array reindexed(typename basic_array::index first, Indexes... idxs)&&{
		return ((std::move(*this).reindexed(first)<<1).reindexed(idxs...))>>1;
	}
private:
	constexpr basic_array sliced_aux(index first, index last) const{
		typename types::layout_t new_layout = *this;
		if((this->size())==0){
			assert(first == last);
			new_layout.nelems_ = 0;
		}else{
			(new_layout.nelems_/=(this->size()))*=(last - first);
		}
		return {new_layout, types::base_ + Layout::operator()(first)};
	}
public:
	constexpr basic_const_array sliced(index first, index last) const&{return sliced_aux(first, last);}
	constexpr basic_array       sliced(index first, index last)      &{return sliced_aux(first, last);}
	constexpr basic_array       sliced(index first, index last)     &&{return sliced_aux(first, last);}

	constexpr basic_const_array blocked(typename basic_array::index first, typename basic_array::index last) const&{return sliced(first, last).reindexed(first);}
	constexpr basic_array blocked(typename basic_array::index first, typename basic_array::index last)&{return sliced(first, last).reindexed(first);}

	using iextension = typename basic_array::index_extension;
	NODISCARD("no side effects")
	constexpr basic_array stenciled(iextension x)                                             &{return blocked(x.start(), x.finish());}
	constexpr basic_array stenciled(iextension x, iextension x1)                              &{return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2)               &{return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2, iextension x3)&{return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}
	template<class... Xs>
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)&{return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;}

	NODISCARD("no side effects")
	constexpr basic_array stenciled(iextension x)                                             &&{return blocked(x.start(), x.finish());}
	constexpr basic_array stenciled(iextension x, iextension x1)                              &&{return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2)               &&{return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2, iextension x3)&&{return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}
	template<class... Xs>
	constexpr basic_array stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)&&{return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;}

	NODISCARD("no side effects")
	constexpr basic_const_array stenciled(iextension x)                                             const&{return blocked(x.start(), x.finish());}
	constexpr basic_const_array stenciled(iextension x, iextension x1)                              const&{return ((stenciled(x)<<1).stenciled(x1))>>1;}
	constexpr basic_const_array stenciled(iextension x, iextension x1, iextension x2)               const&{return ((stenciled(x)<<1).stenciled(x1, x2))>>1;}
	constexpr basic_const_array stenciled(iextension x, iextension x1, iextension x2, iextension x3)const&{return ((stenciled(x)<<1).stenciled(x1, x2, x3))>>1;}
	template<class... Xs>
	constexpr basic_const_array stenciled(iextension x, iextension x1, iextension x2, iextension x3, Xs... xs)const&{return ((stenciled(x)<<1).stenciled(x1, x2, x3, xs...))>>1;}

	constexpr decltype(auto) elements_at(size_type n) const&{assert(n < this->num_elements()); 
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}
	constexpr decltype(auto) elements_at(size_type n) &&{assert(n < this->num_elements()); 
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}
	constexpr decltype(auto) elements_at(size_type n) &{assert(n < this->num_elements()); 
		auto const sub_num_elements = this->begin()->num_elements();
		return operator[](n / sub_num_elements).elements_at(n % sub_num_elements);
	}

	constexpr basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = *this; 
		new_layout.stride_*=s;
		return {new_layout, types::base_};
	}
	constexpr basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}
	using index_range = typename basic_array::index_range;
	constexpr decltype(auto) range(index_range ir) &     {return sliced(ir.front(), ir.front() + ir.size());}
	constexpr decltype(auto) range(index_range ir) &&    {return range(ir);}
	constexpr decltype(auto) range(index_range ir) const&{return sliced(ir.front(), ir.front() + ir.size());}

	constexpr auto range(typename types::index_range const& ir, dimensionality_type n) const{
		return rotated(n).range(ir).rotated(-n);
	}
	constexpr decltype(auto) flattened()&&{
		multi::biiterator<std::decay_t<decltype(std::move(*this).begin())>> biit{std::move(*this).begin(), 0, size(*(std::move(*this).begin()))};
		return basic_array<typename std::iterator_traits<decltype(biit)>::value_type, 1, decltype(biit)>{
			multi::layout_t<1>(1, 0, this->size()*size(*(std::move(*this).begin()))),
			biit
		};
	}
	friend constexpr decltype(auto) flattened(basic_array&& self){return std::move(self).flattened();}
	constexpr bool is_flattable() const{return this->stride() == this->layout().sub_.nelems_;}
	constexpr auto flatted() const{
		assert(is_flattable() && "flatted doesn't work for all layouts!");//this->nelems());
		multi::layout_t<D-1> new_layout{this->layout().sub_};
		new_layout.nelems_*=this->size();
		return basic_array<T, D-1, ElementPtr>{new_layout, types::base_};
	}
	friend constexpr auto flatted(basic_array const& self){return self.flatted();}

	NODISCARD("because it has no side-effect")
	constexpr auto diagonal()&&{return this->diagonal();}
	NODISCARD("because it has no side-effect")
	constexpr basic_array<T, D-1, typename basic_array::element_ptr> diagonal()&{
		auto L = std::min(std::get<0>(this->sizes()), std::get<1>(this->sizes()));
		multi::layout_t<D-1> new_layout{(*this)({0, L}, {0, L}).layout().sub_};
		new_layout.nelems_ += (*this)({0, L}, {0, L}).layout().nelems_;
		new_layout.stride_ += (*this)({0, L}, {0, L}).layout().stride_;
		return {new_layout, types::base_};
	}
	NODISCARD("because it has no side-effect")
	constexpr basic_array<T, D-1, typename basic_array::element_const_ptr> diagonal() const&{
		auto L = std::min(std::get<0>(this->sizes()), std::get<1>(this->sizes()));
		multi::layout_t<D-1> new_layout{(*this)({0, L}, {0, L}).layout().sub_};
		new_layout.nelems_ += (*this)({0, L}, {0, L}).layout().nelems_;
		new_layout.stride_ += (*this)({0, L}, {0, L}).layout().stride_;
		return {new_layout, types::base_};
	}
	friend constexpr auto diagonal(basic_array const& self){return           self .diagonal();}
	friend constexpr auto diagonal(basic_array&       self){return           self .diagonal();}
	friend constexpr auto diagonal(basic_array&&      self){return std::move(self).diagonal();}

	using partitioned_type       = basic_array<T, D+1, element_ptr      >;
	using partitioned_const_type = basic_array<T, D+1, element_const_ptr>;
private:
	constexpr partitioned_type partitioned_aux(size_type s) const{
		assert(s != 0);
		assert(this->layout().nelems_%s==0);
		multi::layout_t<D+1> new_layout{this->layout(), this->layout().nelems_/s, 0, this->layout().nelems_};
		new_layout.sub_.nelems_/=s;
		return {new_layout, types::base_};
	}
public:
	constexpr partitioned_const_type partitioned(size_type s) const&{return partitioned_aux(s);}
	constexpr partitioned_type       partitioned(size_type s)      &{return partitioned_aux(s);}
	constexpr partitioned_type       partitioned(size_type s)     &&{return partitioned_aux(s);}
	
	friend constexpr partitioned_const_type partitioned(basic_array const& self, size_type s){return           self .partitioned(s);}
	friend constexpr partitioned_type       partitioned(basic_array      & self, size_type s){return           self .partitioned(s);}
	friend constexpr partitioned_type       partitioned(basic_array     && self, size_type s){return std::move(self).partitioned(s);}

private:
	constexpr basic_array reversed_aux() const{
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}
public:
	constexpr basic_const_array reversed() const&{return reversed_aux();}
	constexpr basic_array       reversed()      &{return reversed_aux();}
	constexpr basic_array       reversed()     &&{return reversed_aux();}

	friend constexpr basic_const_array reversed(basic_array const& s){return           s .reversed();}
	friend constexpr basic_array       reversed(basic_array      & s){return           s .reversed();}
	friend constexpr basic_array       reversed(basic_array     && s){return std::move(s).reversed();}

	constexpr basic_array transposed() const&{//	typename types::layout_t new_layout = *this;
		return {this->layout().transpose(), types::base_};
	}
	friend constexpr basic_array transposed(basic_array const& s){return s.transposed();}
	friend constexpr basic_array operator~ (basic_array const& s){return s.transposed();}

	constexpr basic_array rotated()&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	constexpr basic_array rotated()&&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	constexpr basic_const_array rotated() const&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		typename basic_const_array::element_ptr new_base_{types::base_};
		return basic_const_array{new_layout, new_base_};
	}
	friend constexpr basic_const_array rotated(basic_array const&  self){return self.rotated();}
	friend constexpr basic_array       rotated(basic_array      && self){return std::move(self).rotated();}
	friend constexpr basic_array       rotated(basic_array      &  self){return self.rotated();}

	constexpr auto unrotated() &{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	constexpr auto unrotated() &&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	constexpr auto unrotated() const&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_const_array{new_layout, types::base_};
	}
	friend constexpr auto unrotated(basic_array const& self){return self.unrotated();}

	constexpr basic_array rotated(dimensionality_type i) &{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}
	constexpr basic_array       rotated(dimensionality_type i) &&{return rotated(i);}
	constexpr basic_const_array rotated(dimensionality_type i) const&{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}

	constexpr basic_array unrotated(dimensionality_type i) &{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}
	constexpr basic_array       unrotated(dimensionality_type i) &&{return unrotated(i);}
	constexpr basic_const_array unrotated(dimensionality_type i) const&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}

	constexpr decltype(auto) operator<<(dimensionality_type i)      &{return                    rotated(i);}
	constexpr decltype(auto) operator>>(dimensionality_type i)      &{return                  unrotated(i);}
	constexpr decltype(auto) operator<<(dimensionality_type i)     &&{return std::move(*this).  rotated(i);}
	constexpr decltype(auto) operator>>(dimensionality_type i)     &&{return std::move(*this).unrotated(i);}
	constexpr decltype(auto) operator<<(dimensionality_type i) const&{return                    rotated(i);}
	constexpr decltype(auto) operator>>(dimensionality_type i) const&{return                  unrotated(i);}

	constexpr decltype(auto) operator|(typename basic_array::size_type n) &{return partitioned(n);}
	constexpr decltype(auto) operator|(typename basic_array::size_type n) &&{return std::move(*this).partitioned(n);}
	constexpr decltype(auto) operator|(typename basic_array::size_type n) const&{return partitioned(n);}

	HD constexpr basic_array       operator()() &     {return *this;}
	HD constexpr basic_array       operator()() &&    {return this->operator()();}
	HD constexpr basic_const_array operator()() const&{return {this->layout(), this->base()};}

public:
	template<typename, dimensionality_type, typename, class> friend struct basic_array;
	constexpr basic_array       paren() &     {return *this;}
	constexpr basic_array       paren() &&    {return this->operator()();}
	constexpr basic_const_array paren() const&{return {this->layout(), this->base()};}

	template<class... As> HD constexpr auto paren(index_range a, As... as) &     {return                  range(a).rotated().paren(as...).unrotated();}
	template<class... As> HD constexpr auto paren(index_range a, As... as) &&    {return this->range(a).rotated().paren(as...).unrotated();}
	template<class... As> HD constexpr auto paren(index_range a, As... as) const&{return                  range(a).rotated().paren(as...).unrotated();}

	template<class... As> HD constexpr decltype(auto) paren(intersecting_range<index> inr, As... as) &     {return                  paren(intersection(this->extension(), inr), as...);}
	template<class... As> HD constexpr decltype(auto) paren(intersecting_range<index> inr, As... as) &&    {return 				paren(intersection(this->extension(), inr), as...);}
	template<class... As> HD constexpr decltype(auto) paren(intersecting_range<index> inr, As... as) const&{return                  paren(intersection(this->extension(), inr), as...);}

	template<class... As> HD constexpr decltype(auto) paren(index i, As... as) &     {return                  operator[](i).paren(as...);}
	template<class... As> HD constexpr decltype(auto) paren(index i, As... as) &&    {return                  operator[](i).paren(as...);}
	template<class... As> HD constexpr decltype(auto) paren(index i, As... as) const&{return                  operator[](i).paren(as...);}
public:

	// the default template parameters below help interpret for {first, last} simple syntax as iranges
	// do not remove default parameter = irange
	template<class B1 = irange>                                                                       HD constexpr decltype(auto) operator()(B1 b1)                                const&{return paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr decltype(auto) operator()(B1 b1, B2 b2)                         const&{return paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  const&{return paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) const&{return paren(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       HD constexpr decltype(auto) operator()(B1 b1)                                     &{return paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr decltype(auto) operator()(B1 b1, B2 b2)                              &{return paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                       &{return paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as)      &{return paren(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       HD constexpr decltype(auto) operator()(B1 b1)                                    &&{return paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    HD constexpr decltype(auto) operator()(B1 b1, B2 b2)                             &&{return this->paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                      &&{return paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> HD constexpr decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as)     &&{return paren(b1, b2, b3, b4, as...);}

//	template<class B1 = iextension>                                                                   decltype(auto) block(B1 b1)                                const&{return block_aux(b1);}
//	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         const&{return paren(b1, b2);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  const&{return paren(b1, b2, b3);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) const&{return paren(b1, b2, b3, b4, as...);}

//	template<class B1 = irange>                                                                       decltype(auto) block(B1 b1)                                     &{return block_aux(b1);}
//	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         &{return paren(b1, b2);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  &{return paren(b1, b2, b3);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) &{return paren(b1, b2, b3, b4, as...);}

//	template<class B1 = irange>                                                                       decltype(auto) block(B1 b1)                                    &&{return std::move(*this).block_aux(b1);}
//	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         &&{return std::move(*this).paren(b1, b2);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  &&{return std::move(*this).paren(b1, b2, b3);}
//	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) &&{return std::move(*this).paren(b1, b2, b3, b4, as...);}

private:
	template<typename Tuple, std::size_t ... I> constexpr decltype(auto) apply_impl(Tuple const& t, std::index_sequence<I...>) const&{return this->operator()(std::get<I>(t)...);}
	template<typename Tuple, std::size_t ... I> constexpr decltype(auto) apply_impl(Tuple const& t, std::index_sequence<I...>)      &{return this->operator()(std::get<I>(t)...);}
	template<typename Tuple, std::size_t ... I> constexpr decltype(auto) apply_impl(Tuple const& t, std::index_sequence<I...>)     &&{return std::move(*this).operator()(std::get<I>(t)...);}
public:
	template<typename Tuple> constexpr decltype(auto) apply(Tuple const& t) const&{return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	template<typename Tuple> constexpr decltype(auto) apply(Tuple const& t)     &&{return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	template<typename Tuple> constexpr decltype(auto) apply(Tuple const& t)      &{return apply_impl(t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}

private:
	using Layout::nelems_;
	using Layout::stride_;
	using Layout::sub_;
public:
	using       iterator = array_iterator<typename types::element, D, typename types::element_ptr      >;//, typename types::reference      >;
	using const_iterator = array_iterator<typename types::element, D, typename types::element_const_ptr>;//, typename types::const_reference>;
private:
	template<class Iterator>
	struct basic_reverse_iterator : 
		std::reverse_iterator<Iterator>,
		boost::multi::totally_ordered2<basic_reverse_iterator<Iterator>, void>
	{
		template<class O, typename = decltype(std::reverse_iterator<Iterator>{base(std::declval<O const&>())})>
		constexpr explicit basic_reverse_iterator(O const& o) : std::reverse_iterator<Iterator>{base(o)}{}
		constexpr basic_reverse_iterator() : std::reverse_iterator<Iterator>{}{}
		constexpr explicit basic_reverse_iterator(Iterator it) : std::reverse_iterator<Iterator>(std::prev(it)){}
		constexpr explicit operator Iterator() const{auto ret = this->base(); if(ret!=Iterator{}) return ++ret; else return Iterator{};}
		constexpr explicit operator bool() const{return bool(this->base());}
		constexpr bool operator==(basic_reverse_iterator const& other) const{return (this->base() == other.base());}
		constexpr typename Iterator::reference operator*() const{return this->current;}
		constexpr typename Iterator::pointer operator->() const{return &this->current;}
		constexpr typename Iterator::reference operator[](typename Iterator::difference_type n) const{return *(this->current - n);}
		constexpr bool operator<(basic_reverse_iterator const& o) const{return o.base() < this->base();}
	};
public:
	using reverse_iterator = basic_reverse_iterator<iterator>;
	using ptr = basic_array_ptr<basic_array, Layout>;

//	ptr operator&() const&{return {this->base_, this->layout()};}
//	constexpr BasicArrayPtr<typename basic_array::element, basic_array::dimensionality, typename basic_array::element_ptr> 
	constexpr ptr addressof() &&{return {this->base_, this->layout()};}
	constexpr ptr operator&() &&{return {this->base_, this->layout()};}
//	ptr operator&() &     {return {this->base_, this->layout()};}

	constexpr iterator begin(dimensionality_type d) &&{
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(0       ), l.sub_, l.stride_};
	}
	constexpr iterator end(dimensionality_type d) &&{
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(l.size()), l.sub_, l.stride_};
	}

private:
	constexpr iterator begin_aux() const{return {types::base_          , sub_, stride_};}
	constexpr iterator end_aux()   const{return {types::base_ + nelems_, sub_, stride_};}

public:
	constexpr iterator begin()          &{return begin_aux();}
	constexpr iterator end  ()          &{return end_aux()  ;}
	friend constexpr iterator begin(basic_array& s){return s.begin();}
	friend constexpr iterator end  (basic_array& s){return s.end  ();}

	constexpr iterator begin()          &&   {return              begin();}
	constexpr iterator end  ()          &&   {return              end()  ;}
	friend constexpr iterator begin(basic_array&& s){return std::move(s).begin();}
	friend constexpr iterator end  (basic_array&& s){return std::move(s).end()  ;}

	constexpr const_iterator begin()           const&{return begin_aux();}
	constexpr const_iterator end  ()           const&{return end_aux()  ;}
	friend constexpr const_iterator begin(basic_array const& s){return s.begin();}
	friend constexpr const_iterator end  (basic_array const& s){return s.end()  ;}
	
	constexpr const_iterator cbegin() const{return begin();}
	constexpr const_iterator cend()   const{return end()  ;}
	friend constexpr auto cbegin(basic_array const& s){return s.cbegin();}
	friend constexpr auto cend  (basic_array const& s){return s.cend()  ;}

	template<class It> constexpr It assign(It first) &{adl_copy_n(first, this->size(), begin()); std::advance(first, this->size()); return first;}
	template<class It> constexpr It assign(It first)&&{return assign(first);}

	template<class Range, class = std::enable_if_t<not std::is_base_of<basic_array, Range>{}> >
//	constexpr
	auto operator=(Range const& r)& // check that you LHS is not read-only
	->decltype(assign(adl_begin(r)), std::declval<basic_array&>()){assert(this->size() == r.size());
		MULTI_MARK_SCOPE(std::string{"multi::operator= D="}+std::to_string(D)+" from range to "+typeid(T).name() );
		assign(adl_begin(r));
		return *this;
	}

	template<class Range, class = std::enable_if_t<not std::is_base_of<basic_array, Range>{}> >
	basic_array&& operator=(Range const& r)&&{return std::move(operator=(r));}

	template<class TT, class... As>
//	constexpr 
	basic_array& operator=(basic_array<TT, D, As...> const& o)&{assert( this->extension() == o.extension() );
		MULTI_MARK_SCOPE(std::string{"multi::operator= "}+std::to_string(D)+" from "+typeid(TT).name()+" to "+typeid(T).name() );
		if(this->num_elements() == this->nelems() and o.num_elements() == this->nelems() and this->layout() == o.layout()){
			adl_copy_n(o.base(), o.num_elements(), this->base());
		}else if(o.stride() < (~o).stride()){
			~(*this) = ~o;
		}else{
			assign(o.begin());
		}
		return *this;
	}
	template<class TT, class... As>
	constexpr basic_array&& operator=(basic_array<TT, D, As...>&& o)&&{return std::move(basic_array::operator=(std::move(o)));}

	template<class TT, class... As>
	constexpr basic_array&& operator=(basic_array<TT, D, As...> const& o)&&{return std::move(this->operator=(o));}

//	constexpr 
	basic_array&  operator=(basic_array               const& o) &{assert( this->extension() == o.extension() );
		MULTI_MARK_SCOPE("multi::operator= D="+std::to_string(D)+" from "+typeid(T).name()+" to "+typeid(T).name() );
		if(this->num_elements() == this->nelems() and o.num_elements() == this->nelems() and this->layout() == o.layout()){
			adl_copy_n(o.base(), o.num_elements(), this->base());
		}else if(o.stride() < (~o).stride()){
			~(*this) = ~o;
		}else{
			assign(o.begin());
		}
		return *this;
	}
	constexpr basic_array&& operator=(basic_array const& o) &&{return std::move(operator=(o));}
	template<class Array> void swap(Array&& o) &&{assert( std::move(*this).extension() == std::forward<Array>(o).extension() );
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> constexpr void swap(A&& o) &{return swap(std::forward<A>(o));}

	friend constexpr void swap(basic_array&& a, basic_array&& b){std::move(a).swap(std::move(b));}
	template<class Array> constexpr void swap(basic_array const& s, Array&& a){s.swap(a);}
	template<class Array> constexpr void swap(Array&& a, basic_array const& s){s.swap(a);}
	template<class Array>//, std::enable_if_t<std::is_same<Array, basic_array>{}, int> =0> 
	constexpr auto operator==(Array const& o) const&
	->decltype((this->extension()==o.extension()) and adl_equal(this->begin(), this->end(), adl_begin(o))){
		return (this->extension()==o.extension()) and adl_equal(this->begin(), this->end(), adl_begin(o));}
	template<class Array> constexpr bool operator!=(Array const& o) const&{return not((*this)==o);}
	template<class TT, class... As>
	constexpr bool operator==(basic_array<TT, D, As...> const& o) const&{
		return (this->extension()==o.extension()) and adl_equal(this->begin(), this->end(), adl_begin(o));		
	}
	template<class It>
	constexpr bool equal(It begin) const&{
		return adl_equal(
			std::move(modify(*this)).begin(), 
			std::move(modify(*this)).end(),
			begin
		);
	}
private:
	friend constexpr bool lexicographical_compare(basic_array&& a1, basic_array&& a2){
		if(a1.extension().first() > a2.extension().first()) return true;
		if(a1.extension().first() < a2.extension().first()) return false;
		return adl_lexicographical_compare(
			std::move(a1).begin(), std::move(a1).end(), 
			std::move(a2).begin(), std::move(a2).end()
		);
	}
public:
	template<class O> constexpr bool operator<(O&& o)&&{return lexicographical_compare(std::move(*this), std::move(o));}
	template<class O> constexpr bool operator>(O&& o)&&{return lexicographical_compare(std::move(o), std::move(*this));}
public:
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	constexpr basic_array<T2, D, P2> static_array_cast() const{
		P2 p2{this->base_};
		return basic_array<T2, D, P2>{this->layout(), p2};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr basic_array<T2, D, P2> member_cast(PM pm) const&{
		static_assert(sizeof(T)%sizeof(T2) == 0, 
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");
	//	return {this->layout().scale(sizeof(T)/sizeof(T2)), &(this->base_->*pm)};
		return basic_array<T2, D, P2>{this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr basic_array<T2, D, P2> member_cast(PM pm) &{
		static_assert(sizeof(T)%sizeof(T2) == 0, 
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");
	//	return {this->layout().scale(sizeof(T)/sizeof(T2)), &(this->base_->*pm)};
		return basic_array<T2, D, P2>{this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	constexpr basic_array<T2, D, P2> member_cast(PM pm) &&{return this->member_cast<T2, P2, Element, PM>(pm);}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr basic_array<std::decay_t<T2>, D, P2> reinterpret_array_cast() const&{
		static_assert( sizeof(T)%sizeof(T2)== 0, 
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
		auto thisbase = this->base();
		return {
			this->layout().scale(sizeof(T)/sizeof(T2)), 
			static_cast<P2>(static_cast<void*>(thisbase))
		};
	}
	
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, D, P2> reinterpret_array_cast()&{
		static_assert( sizeof(T)%sizeof(T2)== 0, 
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
	//	using void_ptr = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<void>;
		return {
			this->layout().scale(sizeof(T)/sizeof(T2)),
			reinterpret_cast<P2>(this->base())
		//	static_cast<P2>(static_cast<void*>(static_cast<void_ptr>(this->base())))
		};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, D, P2> reinterpret_array_cast()&&{return this->template reinterpret_array_cast<T2, P2>();}

	template<class T2, class P2 = T2*>
	constexpr basic_array<std::decay_t<T2>, D, P2> const_array_cast()&&{
		return {this->layout(), const_cast<P2>(this->base())};
	}
	
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, D + 1, P2> reinterpret_array_cast(size_type n) &{
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");
		assert( sizeof(T) == sizeof(T2)*n );
		auto const thisbase = this->base();
		P2 new_base; std::memcpy((void*)&new_base, (void const*)&thisbase, sizeof(P2)); //reinterpret_cast<P2 const&>(thisbase) // TODO find a better way, fancy pointers wouldn't need reinterpret_cast
		return { 
			layout_t<D+1>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}.rotate(), 
			new_base
		//	reinterpret_cast<P2>(this->base())
		//	static_cast<P2>(static_cast<void*>(this->base()))
		};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, D + 1, P2> reinterpret_array_cast(size_type n) &&{return reinterpret_array_cast<T2, P2>(n);}

	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr basic_array<std::decay_t<T2>, D + 1, P2> reinterpret_array_cast(size_type n) const&{
		static_assert( sizeof(T)%sizeof(T2) == 0,
			"error: reinterpret_array_cast is limited to integral stride values");
		assert( sizeof(T) == sizeof(T2)*n );
		return { 
			layout_t<D+1>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}.rotate(), 
			static_cast<P2>(static_cast<void*>(this->base()))
		};
	}

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Element, typename Ptr>//, typename Ref>
struct array_iterator<Element, 1, Ptr> ://, Ref> : 
	boost::multi::iterator_facade<
		array_iterator<Element, 1, Ptr>, 
		Element, std::random_access_iterator_tag, 
		typename std::iterator_traits<Ptr>::reference, multi::difference_type
	>,
	multi::affine<array_iterator<Element, 1, Ptr>, multi::difference_type>,
	multi::decrementable<array_iterator<Element, 1, Ptr>>,
	multi::incrementable<array_iterator<Element, 1, Ptr>>,
	multi::totally_ordered2<array_iterator<Element, 1, Ptr>, void>
{
	using affine = multi::affine<array_iterator<Element, 1, Ptr>, multi::difference_type>;
	using difference_type = typename affine::difference_type;

	array_iterator() = default;
	array_iterator(array_iterator const&) = default;

	template<class Other, decltype(_implicit_cast<Ptr>(typename Other::pointer{}))* = nullptr>
	// cppcheck-suppress noExplicitConstructor ; because underlying pointer is implicitly convertible
	constexpr           array_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_}{}
	template<class Other, decltype(_explicit_cast<Ptr>(typename Other::pointer{}))* = nullptr> 
	constexpr explicit array_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_}{}

	template<class, dimensionality_type, class> friend struct array_iterator;
	constexpr explicit array_iterator(std::nullptr_t nu)  : data_{nu}, stride_{1}{}
	constexpr explicit array_iterator(Ptr const& p) : data_{p}, stride_{1}{}
	template<class EElement, typename PPtr, 
		typename = decltype(_implicit_cast<Ptr>(std::declval<array_iterator<EElement, 1, PPtr>>().data_))
	>
	constexpr array_iterator(array_iterator<EElement, 1, PPtr> const& other) : data_{other.data_}, stride_{other.stride_}{} 
	explicit constexpr operator bool() const{return static_cast<bool>(this->data_);}
	HD constexpr typename std::iterator_traits<Ptr>::reference operator[](typename array_iterator::difference_type n) const{return *((*this) + n);}
	constexpr Ptr operator->() const{return data_;}
	using element = Element;
	using element_ptr = Ptr;
	using pointer = element_ptr;
	using stride_type = multi::index;
	using rank = std::integral_constant<dimensionality_type, 1>;
	constexpr bool operator<(array_iterator const& o) const{return distance_to(o) > 0;}
	constexpr array_iterator(Ptr d, typename basic_array<Element, 1, Ptr>::index s) : data_{d}, stride_{s}{} // TODO make explicit?
private:
	friend struct basic_array<Element, 1, Ptr>;
	element_ptr data_ = nullptr;
	stride_type stride_ = {1};
	constexpr difference_type distance_to(array_iterator const& other) const{
		assert(stride_==other.stride_ and (other.data_-data_)%stride_ == 0);
		return (other.data_ - data_)/stride_;
	}
public:
	HD constexpr array_iterator operator+(difference_type n) const{array_iterator ret{*this}; ret+=n; return ret;}
	[[deprecated("use base for iterator")]] constexpr element_ptr data() const{return data_;}
	// constexpr here creates problems with intel 19
	       constexpr element_ptr base()              const&   {return data_;}
	friend constexpr element_ptr base(array_iterator const& s){return s.base();}
	constexpr stride_type stride()              const&   {return   stride_;} friend
	constexpr stride_type stride(array_iterator const& s){return s.stride_;}
	constexpr array_iterator& operator++(){data_+=stride_; return *this;}
	constexpr array_iterator& operator--(){data_-=stride_; return *this;}
//	constexpr bool operator==(array_iterator const& o) const{return data_== o.data_;}
//	constexpr bool operator!=(array_iterator const& o) const{return data_!= o.data_;}
	friend constexpr bool operator==(array_iterator const& a, array_iterator const& b){return a.data_ == b.data_;}
	friend constexpr bool operator!=(array_iterator const& a, array_iterator const& b){return not(a==b);}
	HD constexpr typename std::iterator_traits<element_ptr>::reference operator*() const{return *data_;}
	constexpr difference_type operator-(array_iterator const& o) const{return -distance_to(o);}
	constexpr array_iterator& operator+=(difference_type d){data_+=stride_*d; return *this;}
	constexpr array_iterator& operator-=(difference_type d){data_-=stride_*d; return *this;}
};

template<class Element, dimensionality_type D, typename... Ts>
using iterator = array_iterator<Element, D, Ts...>;

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{0}, ElementPtr, Layout> :
	array_types<T, dimensionality_type(0), ElementPtr, Layout>
{
	using types = array_types<T, dimensionality_type{0}, ElementPtr, Layout>;
	using types::types;
	using element = typename types::element;
	using element_ref = typename std::iterator_traits<typename basic_array::element_ptr>::reference;//decltype(*typename basic_array::element_ptr{});
	using element_cref = typename std::iterator_traits<typename basic_array::element_const_ptr>::reference;
//	constexpr 
	basic_array& operator=(element const& e) &{
		MULTI_MARK_SCOPE(std::string{"multi::operator= D=0 from "}+typeid(T).name()+" to "+typeid(T).name() );
		adl_copy_n(&e, 1, this->base_); 
		return *this;
	}
	constexpr basic_array&& operator=(element const& e) &&{return std::move(operator=(e));}
	constexpr bool operator==(element const& e) const&{return adl_equal(&e, &e + 1, this->base_);}
	constexpr bool operator!=(element const& e) const&{return not operator==(e);}

	template<class TT, class=decltype(std::declval<typename basic_array::element>()==std::declval<TT>())>
	constexpr auto operator==(TT const& e) const&
	->decltype(adl_equal(&e, &e + 1, this->base_)){
		return adl_equal(&e, &e + 1, this->base_);}
	template<class TT> constexpr auto operator!=(TT const& e) const&->decltype(!operator==(e)){return !operator==(e);}
	
	template<class Range0>
	basic_array& operator=(Range0&& r)&{
	//	*this->base_ = std::forward<Range0>(r); 
		adl_copy_n(&r, 1, this->base_);
		return *this;
	}
	
	element_cref elements_at(size_type n) const&{assert(n < this->num_elements()); return *(this->base_);}
	element_ref  elements_at(size_type n)     &&{assert(n < this->num_elements()); return *(this->base_);}
	element_ref  elements_at(size_type n)      &{assert(n < this->num_elements()); return *(this->base_);}

	constexpr bool operator==(basic_array const& o) const&{assert(0);
		return adl_equal(o.base_, o.base_ + 1, this->base_);
	}
	constexpr bool operator!=(basic_array const& o) const&{return not operator==(o);}
	constexpr typename basic_array::element_ptr operator&() const{return this->base_;}
	using decay_type = typename types::element;
	constexpr element_ref operator()() const&{return *(this->base_);}
	constexpr operator element_ref()&&{return *(this->base_);}
	constexpr operator typename basic_array::element_type() const&{return *(this->base_);}
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int){
		ar & multi::archive_traits<Archive>::make_nvp("element", *(this->base_));
	}
};

template<typename T, typename ElementPtr, class Layout>
struct basic_array<T, dimensionality_type{1}, ElementPtr, Layout> : 
	multi::partially_ordered2<basic_array<T, dimensionality_type(1), ElementPtr, Layout>, void>,
	multi::random_iterable<basic_array<T, dimensionality_type(1), ElementPtr, Layout> >,
	array_types<T, dimensionality_type(1), ElementPtr, Layout>
{
	using types = array_types<T, dimensionality_type{1}, ElementPtr, Layout>;
	using types::types;
	
	using default_allocator_type = typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type;

	constexpr default_allocator_type get_allocator() const{return default_allocator_of(basic_array::base());}
	friend constexpr default_allocator_type get_allocator(basic_array const& self){return self.get_allocator();}
	using decay_type = array<typename types::element, dimensionality_type{1}, typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type>;
	       constexpr decay_type decay()           const&      {return decay_type{*this};}
	friend constexpr decay_type decay(basic_array const& self){return self.decay();}
	using basic_const_array = basic_array<
		T, 1, 
		typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>,
		Layout
	>;
	
	using typename types::element_ptr;
	using typename types::element_const_ptr;
protected:
	template<class A>
	constexpr void intersection_assign_(A&& other)&{
		for(auto idx : intersection(types::extension(), extension(other)))
			operator[](idx) = std::forward<A>(other)[idx];
	}
	template<class A> constexpr void intersection_assign_(A&& o)&&{intersection_assign_(std::forward<A>(o));}
protected:
	template<class TT, dimensionality_type DD, typename EP, class LLayout> friend struct basic_array;
	template<class TT, dimensionality_type DD, class Alloc> friend struct static_array;
	basic_array(basic_array const&) = default;
	template<class T2, class P2, class TT, dimensionality_type DD, class PP>
	friend constexpr decltype(auto) static_array_cast(basic_array<TT, DD, PP> const&);
public:
	friend decay_type operator+(basic_array const& self){return self.decay();}
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array&& a){
		return std::move(a).template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
	template<class T2> friend constexpr auto reinterpret_array_cast(basic_array const& a){
		return a.template reinterpret_array_cast<T2, typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>();
	}
public:
	basic_array(basic_array&&) = default; // in C++ < 17 this is necessary to return references from functions
// in c++17 things changed and non-moveable non-copyable types can be returned from functions and captured by auto
protected:
	template<class, class> friend struct basic_array_ptr;
	template<class, dimensionality_type D, class> friend struct array_iterator;
public:
//	using default_allocator_type = typename multi::pointer_traits<typename basic_array::element_ptr>::default_allocator_type;
	friend constexpr dimensionality_type dimensionality(basic_array const& self){return self.dimensionality;}
//	template<class BasicArray, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<BasicArray>>{}>, typename = decltype(types(std::declval<BasicArray&&>()))>
//	constexpr basic_array(BasicArray&& other) : types{std::forward<BasicArray>(other)}{}
//	basic_array_ptr<basic_array, Layout> operator&() const&{return {this->base_, this->layout()};}
	constexpr basic_array_ptr<basic_array, Layout> operator&() &&{return {this->base_, this->layout()};}
//	basic_array_ptr<basic_array, Layout> operator&() &{return {this->base_, this->layout()};}
	constexpr void assign(std::initializer_list<typename basic_array::value_type> il) const{assert( il.size() == static_cast<std::size_t>(this->size()) );
		assign(il.begin(), il.end());
	}
	
	template<class It> 
	constexpr It assign(It first) &{adl_copy_n(first, this->size(), this->begin()); std::advance(first, this->size()); return first;}
	template<class It> 
	constexpr It assign(It first)&&{return assign(first);}

	template<class It>
	constexpr void  assign(It first, It last) &{assert( std::distance(first, last) == this->size() );
		assign(first);
	}
	template<class It> 
	constexpr void assign(It first, It last)&&{assign(first, last);}

	template<class Archive>
	auto serialize(Archive& ar, unsigned){
		std::for_each(this->begin(), this->end(),[&](auto&& e){ar& multi::archive_traits<Archive>::make_nvp("item",e);});
	}
	basic_array& operator=(basic_array const& o)&{assert(this->extension() == o.extension()); 	// TODO make sfinae friendly
		MULTI_MARK_SCOPE(std::string{"multi::operator= D=1 from "}+typeid(T).name()+" to "+typeid(T).name() );
		this->assign(o.begin(), o.end()); // TODO improve performance by rotating
		return *this;
	} // TODO leave only r-value version?
	template<class TT, dimensionality_type DD, class... As>
	constexpr basic_array&& operator=(basic_array const& o)&&{return std::move(this->operator=(o));} 	// TODO make sfinae friendly

	HD constexpr typename basic_array::const_reference operator[](index i) const&{MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"out of bounds");
		return *(this->base() + Layout::operator()(i)); // in C++17 this is allowed even with syntethic references
	}
	HD constexpr typename basic_array::      reference operator[](index i)      &{//MULTI_ACCESS_ASSERT(this->extension().contains(i)&&"\nout of bounds");;
		return *(this->base() + Layout::operator()(i));
	}
	HD constexpr typename basic_array::reference operator[](index i)&&{return this->operator[](i);}

	template<class Self, typename Tuple, std::size_t ... I> 
	friend HD constexpr decltype(auto) apply_impl(Self&& self, Tuple const& t, std::index_sequence<I...>, basic_array* = 0){return std::forward<Self>(self)(std::get<I>(t)...);}
	template<typename Tuple> HD constexpr decltype(auto) apply(Tuple const& t) const&{return apply_impl(          *this , t, std::make_index_sequence<std::tuple_size<Tuple>::value>());} // TODO tuple_size_v in C++17
	template<typename Tuple> HD constexpr decltype(auto) apply(Tuple const& t)     &&{return apply_impl(std::move(*this), t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}
	template<typename Tuple> HD constexpr decltype(auto) apply(Tuple const& t)      &{return apply_impl(          *this , t, std::make_index_sequence<std::tuple_size<Tuple>::value>());}

	template<class Tuple, typename = std::enable_if_t<(std::tuple_size<std::decay_t<Tuple>>{}>1) > >
	HD constexpr auto operator[](Tuple&& t) const
	->decltype(operator[](std::get<0>(t))[detail::tuple_tail(t)]){
		return operator[](std::get<0>(t))[detail::tuple_tail(t)];}
	template<class Tuple, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tuple>>{}==1> >
	HD constexpr decltype(auto) operator[](Tuple&& t) const{return operator[](std::get<0>(t));}
	HD constexpr decltype(auto) operator[](std::tuple<>) const{return *this;}

	HD constexpr decltype(auto) elements_at(size_type n) const&{assert(n < this->num_elements()); return operator[](n);}
	HD constexpr decltype(auto) elements_at(size_type n)     &&{assert(n < this->num_elements()); return operator[](n);}
	HD constexpr decltype(auto) elements_at(size_type n)      &{assert(n < this->num_elements()); return operator[](n);}

	using typename types::index;
	constexpr basic_array reindexed(typename basic_array::index first)&&{
		typename types::layout_t new_layout = *this;
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
	constexpr basic_array reindexed(typename basic_array::index first)&{
		typename types::layout_t new_layout = *this;
		new_layout.reindex(first);
		return {new_layout, types::base_};
	}
private:
	constexpr basic_array sliced_aux(index first, index last) const{
		typename types::layout_t new_layout = *this;
		if(Layout::size()==0){
			assert(first == last);
			new_layout.nelems_ = 0;
		}else{
			(new_layout.nelems_/=Layout::size())*=(last - first);
		}
		return {new_layout, types::base_ + Layout::operator()(first)};
	}
public:
	constexpr basic_const_array sliced(index first, index last) const&{return sliced_aux(first, last);}
	constexpr basic_array       sliced(index first, index last)      &{return sliced_aux(first, last);}
	constexpr basic_array       sliced(index first, index last)     &&{return sliced_aux(first, last);}

	constexpr basic_array blocked(typename basic_array::index first, typename basic_array::index last)&{return sliced(first, last).reindexed(first);}
	constexpr basic_array stenciled(typename basic_array::index_extension x){return blocked(x.start(), x.finish());}
//	constexpr basic_array sliced(typename types::index first, typename types::index last)&&{return sliced(first, last);}
//	constexpr basic_const_array sliced(typename types::index first, typename types::index last) const&{
//		typename types::layout_t new_layout = *this; 
//		(new_layout.nelems_/=Layout::size())*=(last - first);
//		return {new_layout, types::base_ + Layout::operator()(first)};		
//	}
	constexpr basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = this->layout();
		new_layout.stride_*=s;
		return {new_layout, types::base_};//+ Layout::operator()(this->extension().front())};
	}

	constexpr basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}

	constexpr auto range(index_range const& ir)      &{return sliced(ir.front(), ir.last());}
	constexpr auto range(index_range const& ir)     &&{return std::move(*this).sliced(ir.front(), ir.last());}
	constexpr auto range(index_range const& ir) const&{return sliced(ir.front(), ir.last());}

	constexpr basic_const_array operator()() const&{return {this->layout(), this->base()};}
	constexpr basic_array       operator()()     &&{return *this;}
	constexpr basic_array       operator()()      &{return *this;}

	constexpr auto operator()(index_range const& ir) &{return range(ir);}
	constexpr auto operator()(index_range const& ir) &&{return range(ir);}
	constexpr auto operator()(index_range const& ir) const&{return range(ir);}

	HD constexpr decltype(auto) operator()(index i) &     {return operator[](i);}
	HD constexpr decltype(auto) operator()(index i) &&    {return operator[](i);}
	HD constexpr decltype(auto) operator()(index i) const&{return operator[](i);}

	HD constexpr auto paren() &{return operator()();}
	HD constexpr auto paren() &&{return operator()();}
	HD constexpr auto paren() const&{return operator()();}

	HD constexpr auto paren(index_range const& ir) &{return range(ir);}
	HD constexpr auto paren(index_range const& ir) &&{return range(ir);}
	HD constexpr auto paren(index_range const& ir) const&{return range(ir);}

	HD constexpr decltype(auto) paren(index i) &     {return operator[](i);}
	HD constexpr decltype(auto) paren(index i) &&    {return operator[](i);}
	HD constexpr decltype(auto) paren(index i) const&{return operator[](i);}

public:
	using partitioned_type       = basic_array<T, 2, element_ptr      >;
	using partitioned_const_type = basic_array<T, 2, element_const_ptr>;
private:
	constexpr partitioned_type partitioned_aux(size_type s) const{
		assert( s != 0 );
		assert( this->layout().nelems_%s==0 ); // TODO remove assert? truncate left over? (like mathematica)
		multi::layout_t<2> new_layout{this->layout(), this->layout().nelems_/s, 0, this->layout().nelems_};
		new_layout.sub_.nelems_/=s;
		return {new_layout, types::base_};
	}
public:
	constexpr partitioned_const_type partitioned(size_type s) const&{return partitioned_aux(s);}
	constexpr partitioned_type       partitioned(size_type s)      &{return partitioned_aux(s);}
	constexpr partitioned_type       partitioned(size_type s)     &&{return partitioned_aux(s);}
	
private:
	constexpr basic_array reversed_aux() const{
		auto new_layout = this->layout();
		new_layout.reverse();
		return {new_layout, types::base_};
	}
public:
	constexpr basic_const_array reversed() const&{return reversed_aux();}
	constexpr basic_array       reversed()      &{return reversed_aux();}
	constexpr basic_array       reversed()     &&{return reversed_aux();}

	friend constexpr basic_const_array reversed(basic_array const& s){return           s .reversed();}
	friend constexpr basic_array       reversed(basic_array      & s){return           s .reversed();}
	friend constexpr basic_array       reversed(basic_array     && s){return std::move(s).reversed();}

	friend constexpr decltype(auto)   rotated(basic_array const& s){return s.  rotated();}
	friend constexpr decltype(auto) unrotated(basic_array const& s){return s.unrotated();}

	constexpr decltype(auto)   rotated(dimensionality_type = 1) &     {return operator()();}
	constexpr decltype(auto)   rotated(dimensionality_type = 1) &&    {return operator()();}
	constexpr decltype(auto)   rotated(dimensionality_type = 1) const&{return operator()();}

	constexpr decltype(auto) unrotated(dimensionality_type = 1) &     {return operator()();}
	constexpr decltype(auto) unrotated(dimensionality_type = 1) &&    {return operator()();}
	constexpr decltype(auto) unrotated(dimensionality_type = 1) const&{return operator()();}

	constexpr decltype(auto) operator<<(dimensionality_type i) const{return   rotated(i);}
	constexpr decltype(auto) operator>>(dimensionality_type i) const{return unrotated(i);}

	using       iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_ptr      >;//, typename types::reference>;
	using const_iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_const_ptr>;
	using reverse_iterator = std::reverse_iterator<iterator>;

private:
	constexpr       iterator begin_aux() const{return {this->base_                 , this->stride_};}
	constexpr       iterator end_aux  () const{return {this->base_ + types::nelems_, this->stride_};}
public:
	constexpr const_iterator begin()const&{return begin_aux();}
	constexpr       iterator begin()     &{return begin_aux();}
	constexpr       iterator begin()    &&{return begin_aux();}

	constexpr const_iterator end  ()const&{return end_aux();}
	constexpr       iterator end  ()     &{return end_aux();}
	constexpr       iterator end  ()    &&{return end_aux();}

	friend constexpr const_iterator begin(basic_array const& s){return           s .begin();}
	friend constexpr       iterator begin(basic_array      & s){return           s .begin();}
	friend constexpr       iterator begin(basic_array     && s){return std::move(s).begin();}

	friend constexpr const_iterator end  (basic_array const& s){return           s .end();}
	friend constexpr       iterator end  (basic_array      & s){return           s .end();}
	friend constexpr       iterator end  (basic_array     && s){return std::move(s).end();}

	constexpr const_iterator cbegin() const{return begin();}
	constexpr const_iterator cend  () const{return end()  ;}

	friend constexpr auto cbegin(basic_array const& s){return s.cbegin();}
	friend constexpr auto cend  (basic_array const& s){return s.cend()  ;}

	template<class TT, class... As>//, DELETE((not std::is_assignable<typename basic_array::reference, typename basic_array<TT, 1, As...>::reference>{}))>
//	constexpr 
	auto operator=(basic_array<TT, 1, As...> const& other)&&
	->decltype(adl_copy(other.begin(), other.end(), std::declval<iterator>()), std::declval<basic_array&&>()){assert(this->extensions() == other.extensions());
		MULTI_MARK_SCOPE(std::string{"multi::operator= D=1 from "}+typeid(TT).name()+" to "+typeid(T).name() );
		return adl_copy(other.begin(), other.end(), this->begin()                                 ), std::move(*this);             }

	template<class TT, class... As>//, DELETE((not std::is_assignable<typename basic_array::reference, typename basic_array<TT, 1, As...>::reference>{}))>
	basic_array&  operator=(basic_array<TT, 1, As...> const& other)&{assert(this->extensions() == other.extensions());
		adl_copy(other.begin(), other.end(), this->begin());
		return *this;
	}

	template<class It> constexpr auto assign(It f)&& //	->decltype(adl::copy_n(f, this->size(), begin(std::move(*this))), void()){
	->decltype(adl_copy_n(f, this->size(), std::declval<iterator>()), void()){
		return adl_copy_n(f, this->size(), std::move(*this).begin()), void();}

	template<typename Array>//, typename = std::enable_if_t<not std::is_base_of<basic_array, Array>{}> >
	constexpr bool operator==(Array const& o) const&{ // TODO assert extensions are equal?
		return (this->extension()==extension(o)) and adl_equal(this->begin(), this->end(), adl_begin(o));
	}
	
	constexpr bool operator<(basic_array const& o) const&{return lexicographical_compare(*this, o);}//operator< <basic_array const&>(o);}
	template<class Array> constexpr void swap(Array&& o)&&{assert(this->extension() == o.extension());
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> constexpr void swap(A&& o)&{return swap(std::forward<A>(o));}
	friend constexpr void swap(basic_array&& a, basic_array&& b){std::move(a).swap(std::move(b));}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend constexpr void swap(basic_array&& s, A&& a){s.swap(a);}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend constexpr void swap(A&& a, basic_array&& s){s.swap(a);}
private:
	template<class A1, class A2>
	static constexpr auto lexicographical_compare(A1 const& a1, A2 const& a2){
	//	using multi::extension;
		if(extension(a1).first() > extension(a2).first()) return true;
		if(extension(a1).first() < extension(a2).first()) return false;
		return adl_lexicographical_compare(adl_begin(a1), adl_end(a1), adl_begin(a2), adl_end(a2));
	}
public:
	template<class O>
	constexpr bool operator<(O const& o) const{return lexicographical_compare(*this, o);}
	template<class O>
	constexpr bool operator>(O const& o) const{return lexicographical_compare(o, *this);}
public:
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	constexpr basic_array<T2, 1, P2> static_array_cast() const{//(basic_array&& o){  // name taken from std::static_pointer_cast
		return {this->layout(), static_cast<P2>(this->base())};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>, class... Args>
	constexpr basic_array<T2, 1, P2> static_array_cast(Args&&... args) const{//(basic_array&& o){  // name taken from std::static_pointer_cast
		return {this->layout(), P2{this->base(), std::forward<Args>(args)...}};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 std::decay_t<Element>::*
	>
	constexpr basic_array<T2, 1, P2> member_cast(PM pm) const{
		static_assert(sizeof(T)%sizeof(T2) == 0, 
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");
#if defined(__GNUC__) and (not defined(__INTEL_COMPILER))
		auto&& r1 = (*((typename basic_array::element_type*)(basic_array::base_))).*pm;//->*pm;
		auto p1 = &r1; auto p2 = (P2)p1;
		return {this->layout().scale(sizeof(T)/sizeof(T2)), p2};
#else
		return {this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))}; // this crashes the gcc compiler
#endif
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	basic_array<std::decay_t<T2>, 1, P2> reinterpret_array_cast() const&{
		static_assert( sizeof(T)%sizeof(T2)== 0, "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
//			this->layout().scale(sizeof(T)/sizeof(T2));
		static_assert( sizeof(P2) == sizeof(typename basic_array::element_ptr), "reinterpret on equal size?");
		auto const thisbase = this->base();
		P2 new_base; std::memcpy((void*)&new_base, (void const*)&thisbase, sizeof(P2)); //reinterpret_cast<P2 const&>(thisbase) // TODO find a better way, fancy pointers wouldn't need reinterpret_cast
		return {this->layout().scale(sizeof(T)/sizeof(T2)), new_base};
	}
	
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2 const> >
	constexpr basic_array<std::decay_t<T2>, 2, P2> reinterpret_array_cast(size_type n) const&{
		static_assert( sizeof(T)%sizeof(T2)== 0, 
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
	//	assert( sizeof(T )%(sizeof(T2)*n)== 0 );
		auto thisbase = this->base();		
		return basic_array<std::decay_t<T2>, 2, P2>{
			layout_t<2>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}, 
			static_cast<P2>(static_cast<void*>(thisbase))
		}.rotated();
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, 2, P2> reinterpret_array_cast(size_type n)&{
		static_assert( sizeof(T)%sizeof(T2)== 0, 
			"error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
	//	assert( sizeof(T )%(sizeof(T2)*n)== 0 );
		auto thisbase = this->base();		
		return basic_array<std::decay_t<T2>, 2, P2>{
			layout_t<2>{this->layout().scale(sizeof(T)/sizeof(T2)), 1, 0, n}, 
			static_cast<P2>(static_cast<void*>(thisbase))
		}.rotated();
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	constexpr basic_array<std::decay_t<T2>, 2, P2> reinterpret_array_cast(size_type n)&&{return this->reinterpret_array_cast(n);}

	template<class TT = typename basic_array::element_type>
	constexpr decltype(auto) fill(TT const& value = TT{})&
//	->decltype(adl_fill_n(this->begin(), typename basic_array::size_type{}, value), *this){
	{	return adl_fill_n(this->begin(), this->size(), value), *this;}
	template<class TT = typename basic_array::element_type>
	constexpr decltype(auto) fill(TT const& value = TT{})&&
//	->decltype(std::move(this->fill(value))){
	{	return std::move(this->fill(value));}
};

template<class T2, class P2, class Array, class... Args>
constexpr decltype(auto) static_array_cast(Array&& a, Args&&... args){
	return a.template static_array_cast<T2, P2>(std::forward<Args>(args)...);
}

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref : 
//TODO	multi::partially_ordered2<array_ref<T, D, ElementPtr>, void>,
	basic_array<T, D, ElementPtr>
{
protected:
	constexpr array_ref() noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{}, nullptr}{}
//#if __cplusplus >= 201703L and not defined(__INTEL_COMPILER)
//protected: 
//	[[deprecated("references are not copyable, use &&")]]
//	array_ref(array_ref const&) = default; // don't try to use `auto` for references, use `auto&&` or explicit value type
//#else
//public:
//	array_ref(array_ref const&) = default;
//#endif
//#if defined(__NVCC__)
//	array_ref(array_ref const&) = default;
//	array_ref(array_ref&&) = default;
//#endif
protected:
	[[deprecated("references are not copyable, use auto&&")]]
	array_ref(array_ref const&) = default; // don't try to use `auto` for references, use `auto&&` or explicit value type
public:
	array_ref(array_ref&&) = default; // this needs to be public in c++14
public:
	template<class OtherPtr, class=std::enable_if_t<not std::is_same<OtherPtr, ElementPtr>{}>>
	constexpr array_ref(array_ref<T, D, OtherPtr>&& other)
		: basic_array<T, D, ElementPtr>{other.layout(), ElementPtr{other.base()}}{}
	constexpr explicit array_ref(typename array_ref::element_ptr p, typename array_ref::extensions_type e = {}) noexcept // TODO eliminate this ctor
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p}{}

	constexpr array_ref(typename array_ref::extensions_type e, typename array_ref::element_ptr p) noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p}{}

	template<class TT, std::size_t N> // cppcheck-suppress noExplicitConstructor ; because a reference to c-array can be represented as an array_ref
	constexpr array_ref(TT(&t)[N]) : array_ref((typename array_ref::element_ptr)&t, extensions(t)){}

	using basic_array<T, D, ElementPtr>::operator=;
	using basic_array<T, D, ElementPtr>::operator==;
private:
	template<class It> constexpr auto copy_elements(It first){
		return adl_copy_n(first, array_ref::num_elements(), array_ref::data_elements());
	}
	template<class It> constexpr auto equal_elements(It first) const{
		return adl_equal(first, first + this->num_elements(), this->data_elements());
	}
	template<class TT, std::size_t N> using const_carr = TT const[N];
	template<class TT, std::size_t N> using carr       = TT      [N];
public:
	template<class TT, std::size_t N, std::enable_if_t<std::is_same<typename array_ref::element_type, std::decay_t<std::remove_all_extents_t<const_carr<TT, N>>>>{}, int> =0>
	constexpr operator const_carr<TT, N>&() const&{assert(extensions(*(const_carr<TT, N>*)this)==this->extensions());
		return *reinterpret_cast<const_carr<TT, N>*>(this->base_);
	}
	template<class TT, std::size_t N, std::enable_if_t<std::is_same<typename array_ref::element_type, std::decay_t<std::remove_all_extents_t<carr<TT, N>>> >{}, int> =0>
	constexpr operator carr<TT, N>&()&{assert(extensions(*(carr<TT, N>*)this)==this->extensions());
		return *reinterpret_cast<carr<TT, N>*>(this->base_);
	}
	constexpr typename array_ref::element_ptr data_elements() const&{return array_ref::base_;}
	constexpr array_ref&& operator=(array_ref const& o) &&{assert(this->num_elements()==o.num_elements());
		return array_ref::copy_elements(o.data_elements()), std::move(*this);
	}
	template<typename TT, dimensionality_type DD = D, class... As>
//	constexpr 
	array_ref& operator=(array_ref<TT, DD, As...> const& o)&{assert(this->extensions() == o.extensions());
		MULTI_MARK_SCOPE(std::string{"multi::operator= D="}+std::to_string(D)+" from "+typeid(TT).name()+" to "+typeid(T).name() );
		adl_copy_n(o.data(), o.num_elements(), this->data());
		return *this;
	}
	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr array_ref&& operator=(array_ref<TT, DD, As...> const& o)&&{return std::move(operator=(o));}

	using  elements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_ptr      >;
	using celements_type = array_ref<typename array_ref::element_type, 1, typename array_ref::element_const_ptr>;

private:
	constexpr elements_type elements_() const{return elements_type{data_elements(), this->num_elements()};}
public:
	constexpr  elements_type elements()         &     {return elements_();}
	constexpr  elements_type elements()         &&    {return elements_();}
	constexpr celements_type elements()         const&{return elements_();}

	friend constexpr elements_type elements(array_ref &      self){return           self . elements();}	
	friend constexpr elements_type elements(array_ref &&     self){return std::move(self). elements();}
	friend constexpr celements_type elements(array_ref const& self){return           self . elements();}

	       constexpr celements_type celements()         const&   {return {array_ref::data(), array_ref::num_elements()};}
	friend constexpr celements_type celements(array_ref const& s){return s.celements();}
	
	template<typename TT, dimensionality_type DD = D, class... As>
	constexpr bool operator==(array_ref<TT, DD, As...>&& o) const&{
		if( this->extensions() != o.extensions() ) return false; // TODO, or assert?
		return equal_elements(std::move(o).data_elements());
	}
	       constexpr typename array_ref::element_ptr data_elements()        &&   {return array_ref::base_;}
	friend constexpr typename array_ref::element_ptr data_elements(array_ref&& s){return std::move(s).data_elements();}

	       constexpr typename array_ref::element_ptr data()         const&   {return array_ref::base_;} 
	friend constexpr typename array_ref::element_ptr data(array_ref const& s){return s.data();}

	constexpr typename array_ref::decay_type const& operator*() const&{return static_cast<typename array_ref::decay_type const&>(*this);}
//	constexpr typename array_ref::decay_type const& operator*() const&{return *this;}
	
	constexpr typename array_ref::decay_type const& decay() const&{
		return static_cast<typename array_ref::decay_type const&>(*this);
	}
	friend constexpr typename array_ref::decay_type const& decay(array_ref const& s){return s.decay();}

	template<class Archive>
	auto serialize(Archive& ar, const unsigned int v){
	//	using boost::serialization::make_nvp;
//		if(this->num_elements() < (2<<8) ) 
			basic_array<T, D, ElementPtr>::serialize(ar, v);
//		else{
		//	using boost::serialization::make_binary_object;
		//	using boost::serialization::make_array;
//			if(std::is_trivially_copy_assignable<typename array_ref::element>{})
//				ar & multi::archive_traits<Archive>::make_nvp("binary_data", multi::archive_traits<Archive>::make_binary_object(this->data(), sizeof(typename array_ref::element)*this->num_elements())); //#include<boost/serialization/binary_object.hpp>
//			else ar & multi::archive_traits<Archive>::make_nvp("data", multi::archive_traits<Archive>::make_array(this->data(), this->num_elements()));
//		}
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
struct array_ptr : basic_array_ptr<basic_array<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>{
	using basic_ptr = basic_array_ptr<basic_array<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;
//	using basic_ptr = basic_array_ptr<array_ref<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;
//	using basic_ptr::basic_ptr;//array_ptr<array_ref<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>::basic_array_ptr;
public:
	constexpr array_ptr(Ptr p, index_extensions<D> x) : basic_ptr(p, multi::layout_t<D>{x}){}
	// cppcheck-suppress noExplicitConstructor ; because array_ptr can represent a null
	constexpr array_ptr(std::nullptr_t) : basic_ptr(nullptr, multi::layout_t<D>{}){}
	template<class TT, std::size_t N> // cppcheck-suppress noExplicitConstructor ; because array_ptr can represent a pointer to a c-array
	constexpr array_ptr(TT(*t)[N]) : basic_ptr(data_elements(*t), layout(*t)){}
	constexpr array_ref<T, D, Ptr> operator*() const{
		return array_ref<T, D, Ptr>{this->base(), (*this)->extensions()};//multi::layout_t<D>{x}};
	}
};

template<class T, typename Ptr>
class array_ptr<T, 0, Ptr> : multi::array_ref<T, 0, Ptr>{// Ref_;
public:
//	array_ptr(array_ptr&&) : Ref_{
	constexpr explicit array_ptr(Ptr p, typename multi::array_ref<T, 0, Ptr>::extensions_type x = {}) : multi::array_ref<T, 0, Ptr>(p, x){}
//	operator bool() const{return Ref_.base();}
	constexpr explicit operator Ptr () const{return this->base();}
	friend constexpr bool operator==(array_ptr const& self, array_ptr const& other){return self.base() == other.base();}
	friend constexpr bool operator!=(array_ptr const& self, array_ptr const& other){return self.base() != other.base();}
	constexpr multi::array_ref<T, 0, Ptr>& operator* () const{return const_cast<array_ptr&>(*this);}//               Ref_ ;}
	constexpr multi::array_ref<T, 0, Ptr>* operator->() const{return const_cast<array_ptr*>(this);}//std::addressof(Ref_);}
};

template<class TT, std::size_t N>
// auto operator&(TT(&t)[N]){ // c++ cannot overload & for primitive types
constexpr auto addressof(TT(&t)[N]){
	return array_ptr<
		std::decay_t<std::remove_all_extents_t<TT[N]>>, std::rank<TT[N]>{}, std::remove_all_extents_t<TT[N]>*
	>(&t);
}

template<class T, dimensionality_type D, typename Ptr = T*>
using array_cptr = array_ptr<T, D, 	typename std::pointer_traits<Ptr>::template rebind<T const>>;

template<dimensionality_type D, class P>
constexpr auto make_array_ref(P p, index_extensions<D> x){
	return array_ref<typename std::iterator_traits<P>::value_type, D, P>(p, x);
}

template<class P> auto make_array_ref(P p, index_extensions<1> x){return make_array_ref<1>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<2> x){return make_array_ref<2>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<3> x){return make_array_ref<3>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<4> x){return make_array_ref<4>(p, x);}
template<class P> auto make_array_ref(P p, index_extensions<5> x){return make_array_ref<5>(p, x);}

//In ICC you need to specify the dimensionality in make_array_ref<D>
//#if defined(__INTEL_COMPILER)
//template<dimensionality_type D, class P> 
//auto make_array_ref(P p, std::initializer_list<index_extension> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
//template<dimensionality_type D, class P> 
//auto make_array_ref(P p, std::initializer_list<index> il){return make_array_ref(p, detail::to_tuple<D, index_extension>(il));}
//#endif

#if defined(__cpp_deduction_guides)

template<class It, typename V = typename std::iterator_traits<It>::value_type> // pointer_traits doesn't have ::value_type
array_ptr(It, index_extensions<0> = {})->array_ptr<V, 0, It>;

template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<1>)->array_ptr<V, 1, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<2>)->array_ptr<V, 2, It>;
template<class It, typename V = typename std::iterator_traits<It>::value_type>
array_ptr(It, index_extensions<3>)->array_ptr<V, 3, It>;

template<class T, std::size_t N, typename V = typename std::remove_all_extents<T[N]>::type, std::size_t D = std::rank<T[N]>{}>
array_ptr(T(*)[N])->array_ptr<V, D>;

//#if not defined(__clang__)
//template<class It, dimensionality_type D, typename V = typename std::iterator_traits<It>::value_type>
//array_ref(It, index_extensions<D>)->array_ref<V, D, It>;
//#else
template<class It> array_ref(It, index_extensions<1>)->array_ref<typename std::iterator_traits<It>::value_type, 1, It>;
template<class It> array_ref(It, index_extensions<2>)->array_ref<typename std::iterator_traits<It>::value_type, 2, It>;
template<class It> array_ref(It, index_extensions<3>)->array_ref<typename std::iterator_traits<It>::value_type, 3, It>;
template<class It> array_ref(It, index_extensions<4>)->array_ref<typename std::iterator_traits<It>::value_type, 4, It>;
template<class It> array_ref(It, index_extensions<5>)->array_ref<typename std::iterator_traits<It>::value_type, 5, It>;
//#endif

template<class It, class Tuple> array_ref(It, Tuple)->array_ref<typename std::iterator_traits<It>::value_type, std::tuple_size<Tuple>::value, It>;
#endif

#if 1
template<class T, std::size_t N>
constexpr auto rotated(const T(&t)[N]) noexcept{ // TODO move to utility
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(t))>(
		base(t), extensions(t)
	).rotated();
}
template<class T, std::size_t N>
constexpr auto rotated(T(&t)[N]) noexcept{
	return multi::array_ref<std::remove_all_extents<T[N]>, std::rank<T[N]>{}, decltype(base(t))>(
		base(t), extensions(t)
	).rotated();
}
#endif

template<class TD, class Ptr> struct Array_aux;
template<class T, std::size_t D, class Ptr> struct Array_aux<T   [D], Ptr>{using type = array    <T, D, Ptr>  ;};
template<class T, std::size_t D, class Ptr> struct Array_aux<T(&)[D], Ptr>{using type = array_ref<T, D, Ptr>&&;};
template<class T, std::size_t D, class Ptr> struct Array_aux<T(*)[D], Ptr>{using type = array_ptr<T, D, Ptr>  ;};

template<class TD, class Second = 
	std::conditional_t<
		std::is_reference<TD>{} or std::is_pointer<TD>{}, 
		std::add_pointer_t<std::remove_all_extents_t<std::remove_reference_t<std::remove_pointer_t<TD>>>>,
		std::allocator<std::remove_all_extents_t<TD>>
	>
> using Array = typename Array_aux<TD, Second>::type;

template<class RandomAccessIterator, dimensionality_type D>
constexpr
multi::array_ptr<typename std::iterator_traits<RandomAccessIterator>::value_type, D, RandomAccessIterator>
operator/(RandomAccessIterator data, multi::iextensions<D> x){return {data, x};}

template<class T, dimensionality_type D, class... Ts>
constexpr std::true_type  is_basic_array_aux(basic_array<T, D, Ts...> const&);
constexpr std::false_type is_basic_array_aux(...);

template<class A> struct is_basic_array: decltype(is_basic_array_aux(std::declval<A>())){};

template<class In, class T, dimensionality_type N, class TP, class=std::enable_if_t<(N>1)>, class=decltype(adl_begin(*In{}), adl_end(*In{}))>
constexpr auto uninitialized_copy
// require N>1 (this is important because it forces calling placement new on the pointer
(In first, In last, multi::array_iterator<T, N, TP> dest){
	using std::begin; using std::end;
	while(first!=last){
		adl_uninitialized_copy(adl_begin(*first), adl_end(*first), adl_begin(*dest));
		++first;
		++dest;
	}
	return dest;
}

}}

namespace boost{
namespace multi{

// begin and end for forwarding reference are needed in this namespace 
// to overwrite the behavior of std::begin and std::end 
// which take rvalue-references as const-references.

template<class T> auto begin(T&& t)
->decltype(std::forward<T>(t).begin()){
	return std::forward<T>(t).begin();}

template<class T> auto end(T&& t)
->decltype(std::forward<T>(t).end()){
	return std::forward<T>(t).end();}

}}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include<cassert>

namespace multi = boost::multi;

int main(){

	{
		double a[4][5] = {{1.,2.},{2.,3.}};
		multi::array_ref<double, 2> A(&a[0][0], {4, 5});
		multi::array_ref<double, 2, double const*> B(&a[0][0], {4, 5});
		multi::array_ref<double const, 2> C(&a[0][0], {4, 5});
		multi::array_cref<double, 2> D(&a[0][0], {4, 5});
		A[1][1] = 2.;

		double d[4][5] = {{1.,2.},{2.,3.}};

		auto dd = (double const(&)[4][5])(d);
		assert( &(dd[1][2]) == &(d[1][2]) );
		assert(( & A[1].static_array_cast<double, double const*>()[1] == &A[1][1] ));
		assert(( &multi::static_array_cast<double, double const*>(A[1])[1] == &A[1][1] ));
	}
	{
		double const d2D[4][5] = {{1.,2.},{2.,3.}};
		multi::array_ref<double, 2, const double*> d2Rce(&d2D[0][0], {4, 5});
		assert( &d2Rce[2][3] == &d2D[2][3] );
		assert( d2Rce.size() == 4 );
		assert( num_elements(d2Rce) == 20 );
	}
	{
		std::string const dc3D[4][2][3] = {
			{{"A0a", "A0b", "A0c"}, {"A1a", "A1b", "A1c"}},
			{{"B0a", "B0b", "B0c"}, {"B1a", "B1b", "B1c"}},
			{{"C0a", "C0b", "C0c"}, {"C1a", "C1b", "C1c"}}, 
			{{"D0a", "D0b", "D0c"}, {"D1a", "D1b", "D1c"}}, 
		};
		multi::array_cref<std::string, 3> A(&dc3D[0][0][0], {4, 2, 3});
		assert( num_elements(A) == 24 and A[2][1][1] == "C1b" );
		auto const& A2 = A.sliced(0, 3).rotated()[1].sliced(0, 2).unrotated();
		assert( multi::rank<std::decay_t<decltype(A2)>>{} == 2 and num_elements(A2) == 6 );
		assert( std::get<0>(sizes(A2)) == 3 and std::get<1>(sizes(A2)) == 2 );
		
		auto const& A3 = A({0, 3}, 1, {0, 2});
		assert( multi::rank<std::decay_t<decltype(A3)>>{} == 2 and num_elements(A3) == 6 );
	}
}

#endif
#endif

