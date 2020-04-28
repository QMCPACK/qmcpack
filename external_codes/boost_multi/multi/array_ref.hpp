#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
for a in ./tests/*.cpp; do echo $a; sh $a || break; echo "\n"; done; exit;*/
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2018-2020

#if (defined(__clang__) and defined(__CUDA__)) or defined(__NVCC__)
#define BOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK // see comments https://www.boost.org/doc/libs/1_72_0/boost/utility/result_of.hpp
#endif

#ifndef BOOST_MULTI_ARRAY_REF_HPP
#define BOOST_MULTI_ARRAY_REF_HPP

#include "utility.hpp"

#include "./detail/layout.hpp"
#include "./detail/types.hpp"     // dimensionality_type
#include "./detail/operators.hpp" // random_iterable
#include "./detail/memory.hpp"    // pointer_traits
//#include "./utility/const_iterator.hpp"

#include "./config/NODISCARD.hpp"

//#include<iostream> // debug
#include<boost/pointer_cast.hpp>

#include<algorithm> // copy_n
#include<cstring> // for memset in reinterpret_cast

#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD
#endif

namespace boost{
namespace serialization{
//	template<class Archive> struct archive_traits;
//	template<class> struct nvp;
//	template<class T> const nvp<T> make_nvp(char const* name, T& t);
//	template<class T> class array_wrapper;
//	template<class T, class S> const array_wrapper<T> make_array(T* t, S s);
//	template<class T> 
//	class binary_object;
//	inline auto make_binary_object(const void * t, std::size_t size);
}}

#if 1
namespace boost{
namespace serialization{
	template<class Archive> struct archive_traits;
	template<class> struct nvp;
	template<class T> const nvp<T> make_nvp(char const* name, T& t);
	template<class T> class array_wrapper;
	template<class T, class S> const array_wrapper<T> make_array(T* t, S s);
//	template<class T> 
	class binary_object;
	inline const binary_object make_binary_object(const void * t, std::size_t size);
}}
#endif


namespace boost{
namespace multi{

template<class T> T& modify(T const& t){return const_cast<T&>(t);}

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>> 
struct basic_array;

template<typename T, dimensionality_type D, class A = std::allocator<T>> struct array;

template<typename T, dimensionality_type D, typename ElementPtr = T*, class Layout = layout_t<D>>
struct array_types : Layout{
	using element = T;
	using element_type = element; // this follows more closely https://en.cppreference.com/w/cpp/memory/pointer_traits
	constexpr static dimensionality_type dimensionality = D;
	using element_ptr = ElementPtr;
	using element_const_ptr = typename std::pointer_traits<ElementPtr>::template rebind<element_type const>; //multi::const_iterator<ElementPtr>; 
//	using element_const_ptr = typename multi::iterator_traits<ElementPtr>::rebind_const;
	using layout_t = Layout;
	using value_type = typename std::conditional<
		(dimensionality>1),
		array<element, dimensionality-1, typename pointer_traits<element_ptr>::default_allocator_type>, 
		typename std::conditional<
			dimensionality == 1,
			element,
			element
		>::type
	>::type;
	using decay_type = array<element, dimensionality, typename pointer_traits<element_ptr>::default_allocator_type>;
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

	constexpr element_ptr base() const{return base_;} //	element_const_ptr cbase() const{return base();}
	friend element_ptr base(array_types const& s){return s.base();}
	constexpr layout_t const& layout() const{return *this;}
	friend layout_t const& layout(array_types const& s){return s.layout();}
	element_ptr            origin() const{return base_+Layout::origin();} //	element_const_ptr     corigin() const{return origin();}
	friend decltype(auto)  origin(array_types const& s){return s.origin();} //	friend decltype(auto) corigin(array_types const& s){return s.corigin();}
protected:
	using derived = basic_array<T, D, ElementPtr, Layout>;
	element_ptr base_;
	array_types() = delete;
	constexpr array_types(std::nullptr_t np) : Layout{}, base_{np}{}
	array_types(array_types const&) = default;
public:
	constexpr array_types(layout_t l, element_ptr data): 
		Layout{l}, 
		base_{data}
	{}
//	template<class T2, class P2, class Array> friend decltype(auto) static_array_cast(Array&&);
public://TODO find why this needs to be public and not protected or friend
	template<class ArrayTypes, typename = std::enable_if_t<not std::is_base_of<array_types, std::decay_t<ArrayTypes>>{}>
		, typename = decltype(element_ptr{std::declval<ArrayTypes const&>().base_})
	>
	array_types(ArrayTypes const& a) : Layout{a}, base_{a.base_}{}
	template<typename ElementPtr2, 
		typename = decltype(Layout{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().layout()}),
		typename = decltype(element_ptr{std::declval<array_types<T, D, ElementPtr2, Layout> const&>().base_})
	>
	array_types(array_types<T, D, ElementPtr2, Layout> const& other) : Layout{other.layout()}, base_{other.base_}{}
	template<class T2, dimensionality_type D2, class E2, class L2> friend struct array_types;
};

template<class Ref, class Layout>
struct basic_array_ptr : 
	private Ref,
	boost::multi::iterator_facade<
		basic_array_ptr<Ref, Layout>, void, std::random_access_iterator_tag, 
		Ref const&, typename Layout::difference_type
	>,
	boost::multi::totally_ordered2<basic_array_ptr<Ref, Layout>, void>
{
	using pointer = Ref const*;
	using element_type = typename Ref::decay_type;
	using difference_type = typename Layout::difference_type;

	using value_type = element_type;
	using reference = Ref const&;
	using iterator_category = std::random_access_iterator_tag;

	basic_array_ptr(std::nullptr_t p = nullptr) : Ref{p}{}
	template<class, class> friend struct basic_array_ptr;
	basic_array_ptr(typename Ref::element_ptr p, layout_t<Ref::dimensionality-1> l) HD : Ref{l, p}{}
	basic_array_ptr(typename Ref::element_ptr p, index_extensions<Ref::dimensionality> e) : Ref{p, e}{}

	template<class Array, typename = decltype(typename Ref::element_ptr{typename Array::element_ptr{}})> 
	basic_array_ptr(Array const& o) : Ref{o->layout(), o->base()}{}//, stride_{o.stride_}{}
	basic_array_ptr(basic_array_ptr const& o) : Ref{static_cast<Layout const&>(o), o.base_}{}//, stride_{o.stride_}{}
	basic_array_ptr& operator=(basic_array_ptr const& other){
		this->base_ = other.base_;
		static_cast<Layout&>(*this) = other;
		return *this;
	}
	explicit operator bool() const{return this->base_;}
	constexpr Ref  operator* () const{return *this;}
	constexpr Ref* operator->() const{return  const_cast<basic_array_ptr*>(this);}
	constexpr Ref* operator->(){return  this;}
	constexpr Ref        operator[](difference_type n) const{return *(*this + n);}
//	template<class O> bool operator==(O const& o) const{return equal(o);}
	bool operator<(basic_array_ptr const& o) const{return distance_to(o) > 0;}
	constexpr basic_array_ptr(typename Ref::element_ptr p, Layout l) : Ref{l, p}{}
	template<typename T, dimensionality_type D, typename ElementPtr, class LLayout>
	friend struct basic_array;
	auto base() const{return this->base_;}
	friend auto base(basic_array_ptr const& self){return self.base();}
	using Ref::base_;
	using Ref::layout;
	bool operator==(basic_array_ptr const& o) const{return base_==o.base_ and layout()==o.layout();}
	template<class O> constexpr bool operator==(O const& o) const{return base()==o->base() and layout() == o->layout();}
	template<class O> constexpr bool operator!=(O const& o) const{return not ((*this)==o);}
	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr bool operator==(O const& o, basic_array_ptr const& s){return s==o;}
	template<class O, std::enable_if_t<not std::is_base_of<basic_array_ptr, O>{}, int> =0> friend constexpr bool operator!=(O const& o, basic_array_ptr const& s){return not(o==s);}
protected:
//	bool equal(basic_array_ptr const& o) const{return base_==o.base_ and layout()==o.layout();}
	void increment(){base_ += Ref::nelems();}
	void decrement(){base_ -= Ref::nelems();}
	void advance(difference_type n) HD{base_ += Ref::nelems()*n;}
	difference_type distance_to(basic_array_ptr const& other) const{
		assert( Ref::nelems() == other.Ref::nelems() and Ref::nelems() != 0 );
		assert( (other.base_ - base_)%Ref::nelems() == 0); 
		assert( layout() == other.layout() );
		return (other.base_ - base_)/Ref::nelems();
	}
public:
	basic_array_ptr& operator+=(difference_type n) HD{advance(n); return *this;}
};

template<class Element, dimensionality_type D, typename Ptr, class Ref 
#if 1
= 
	typename std::conditional<
			D != 1,
			basic_array<Element, D-1, 
				typename std::conditional<
					std::is_same<typename std::pointer_traits<Ptr>::element_type, void>{}, 
					typename std::pointer_traits<Ptr>::template rebind<Element>,
					Ptr
				>::type
			>,
			typename std::iterator_traits<Ptr>::reference
	//		typename std::iterator_traits<
	//				typename std::conditional<
	//					std::is_same<typename std::pointer_traits<Ptr>::element_type, void>{}, 
	//					typename std::pointer_traits<Ptr>::template rebind<Element>,
	//					Ptr
	//				>::type
	//		>::reference
		>::type
#endif
>
struct array_iterator;

template<class Element, dimensionality_type D, typename Ptr, class Ref>
struct array_iterator : 
	boost::multi::iterator_facade<
		array_iterator<Element, D, Ptr, Ref>, void, std::random_access_iterator_tag, 
		Ref const&, typename layout_t<D-1>::difference_type
	>,
	multi::decrementable<array_iterator<Element, D, Ptr, Ref>>,
	multi::incrementable<array_iterator<Element, D, Ptr, Ref>>,
	multi::affine<array_iterator<Element, D, Ptr, Ref>, multi::difference_type>,
	multi::totally_ordered2<array_iterator<Element, D, Ptr, Ref>, void>
{
	using difference_type = typename layout_t<D>::difference_type;
	using value_type = typename Ref::decay_type;
	using pointer = Ref*;
	using reference = Ref&&;//Ref const&;
//	using element_type = typename Ref::value_type;
	using iterator_category = std::random_access_iterator_tag;

	using rank = std::integral_constant<dimensionality_type, D>;

	using element = typename Ref::element;
	using element_ptr = typename Ref::element_ptr;
	array_iterator(std::nullptr_t p = nullptr) : ptr_{p}, stride_{1}{}//Ref{p}{}
	template<class, dimensionality_type, class, class> friend struct array_iterator;
	template<class Other, typename = decltype(typename Ref::types::element_ptr{typename Other::element_ptr{}})> 
	array_iterator(Other const& o) : /*Ref{o.layout(), o.base()},*/ ptr_{o.ptr_.base_, o.ptr_.layout()}, stride_{o.stride_}{}
	array_iterator(array_iterator const&) = default;
	array_iterator& operator=(array_iterator const& other){
		ptr_ = other.ptr_;
		stride_ = other.stride_;
		return *this;
	}
	explicit operator bool() const{return static_cast<bool>(ptr_.base_);}
	constexpr Ref operator*() const{/*assert(*this);*/ return {*ptr_};}//return *this;}
	constexpr decltype(auto) operator->() const{/*assert(*this);*/ return ptr_;}//return this;}
	constexpr Ref operator[](difference_type n) const{return *(*this + n);}
	template<class O> bool operator==(O const& o) const{return equal(o);}
	bool operator<(array_iterator const& o) const{return distance_to(o) > 0;}
	array_iterator(typename Ref::element_ptr p, layout_t<D-1> l, index stride) : /*Ref{l, p},*/
		ptr_{p, l}, 
		stride_{stride}
	{}
	template<typename T, dimensionality_type DD, typename ElementPtr, class LLayout>
	friend struct basic_array;
	auto base() const{return ptr_.base_;}//this->base_;}
	friend auto base(array_iterator const& self){return self.base();}
	auto stride() const{return stride_;}
	friend index stride(array_iterator const& s){return s.stride();}
private:
	basic_array_ptr<Ref, layout_t<D-1>> ptr_;
	index stride_ = {1}; // nice non-zero default
	bool equal(array_iterator const& o) const{return ptr_==o.ptr_ and stride_==o.stride_;}//base_==o.base_ && stride_==o.stride_ && ptr_.layout()==o.ptr_.layout();}
	void increment(){ptr_.base_ += stride_;}
	void decrement(){ptr_.base_ -= stride_;}
	void advance(difference_type n) HD{ptr_.base_ += stride_*n;}
	difference_type distance_to(array_iterator const& other) const{
		assert( stride_ == other.stride_);
		assert( stride_ != 0 );
	//	assert( this->stride()==stride(other) and this->stride() );// and (base(other.ptr_) - base(this->ptr_))%stride_ == 0
	//	assert( stride_ == other.stride_ and stride_ != 0 and (other.ptr_.base_-ptr_.base_)%stride_ == 0 and ptr_.layout() == other.ptr_.layout() );
	//	assert( stride_ == other.stride_ and stride_ != 0 and (other.base_ - base_)%stride_ == 0 and layout() == other.layout() );
		return (other.ptr_.base_ - ptr_.base_)/stride_;
	}
//	friend class boost::iterator_core_access;
public:
	array_iterator& operator++(){increment(); return *this;}
	array_iterator& operator--(){decrement(); return *this;}
	constexpr bool operator==(array_iterator const& o) const{return equal(o);}
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
	biiterator(It me, std::ptrdiff_t pos, std::ptrdiff_t stride) HD : me_{me}, pos_{pos}, stride_{stride}{}
	decltype(auto) operator++(){
		++pos_;
		if(pos_==stride_){
			++me_;
			pos_ = 0;
		}
		return *this;
	}
	bool operator==(biiterator const& o) const{return me_==o.me_ and pos_==o.pos_;}
	biiterator& operator+=(multi::difference_type n){me_ += n/stride_; pos_ += n%stride_; return *this;}
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

template<typename T, dimensionality_type D, typename ElementPtr, class Layout /*= layout_t<D>*/ >
struct basic_array : 
	multi::partially_ordered2<basic_array<T, D, ElementPtr, Layout>, void>,
//	multi::random_iterable<basic_array<T, D, ElementPtr, Layout>>,
	array_types<T, D, ElementPtr, Layout>
{
	using types = array_types<T, D, ElementPtr, Layout>;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr >;
	friend struct basic_array<typename types::element, typename Layout::rank{} + 1, typename types::element_ptr&>;
	using types::layout;
	constexpr auto layout() const /*HD*/{return array_types<T, D, ElementPtr, Layout>::layout();}
	using basic_const_array = basic_array<T, D, 
		typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>,
	//	typename multi::iterator_traits<ElementPtr>::rebind_const, 
		Layout
	>;
protected:
	using types::types;
	template<typename, dimensionality_type, class Alloc> friend struct static_array;
	template<class, class> friend struct basic_array_ptr;
#if __cplusplus >= 201703L
protected: basic_array(basic_array&&) = default; // if you need to generate a copy you can't use `auto` here, use `decay`.
#else
public   : basic_array(basic_array&&) = default; // in C++ < 17 this is necessary to return references from functions
#endif
public:
	basic_array(basic_array const&) = default;
	friend constexpr auto dimensionality(basic_array const& self){return self.dimensionality;}
	using typename types::reference;

	auto get_allocator() const{
		using multi::get_allocator;
		return get_allocator(this->base());
	}
	friend auto get_allocator(basic_array const& self){return self.get_allocator();}
//	using decay_type = array<typename types::element, D, decltype(default_allocator_of(std::declval<ElementPtr>()))>;
	template<class P>
	static decltype(auto) get_allocator_(P const& p){
		using multi::default_allocator_of;
		return default_allocator_of(p);
	}
	using decay_type = array<typename types::element_type, D, decltype(get_allocator_(std::declval<ElementPtr>()))>;
//	decay_type 
//	auto
//	static decay_type remake(std::initializer_list<typename basic_array::value_type> il){return decay_type(il);}
//	template<class... As> static auto remake(As&&... as) -> decay_type{return decay_type(std::forward<As>(as)...);}
	template<class Archive>
	auto serialize(Archive& ar, unsigned int /*file version*/){
		using boost::serialization::make_nvp;
		std::for_each(std::move(*this).begin(), std::move(*this).end(), [&](auto&& e){ar & make_nvp("item", e);});
	}

	decay_type decay() const{
		decay_type ret = std::move(modify(*this));
		return ret;
	}
	friend decay_type decay(basic_array const& self){return self.decay();}
	NODISCARD("because an identity decay was created")
	friend auto operator+(basic_array const& self){return self.decay();}

	constexpr typename types::const_reference operator[](index i) const&{ assert( this->extension().contains(i) );
		typename types::element_const_ptr new_base = typename types::element_ptr(this->base()) + std::ptrdiff_t{Layout::operator()(i)};
		return typename types::const_reference(this->layout().sub_, new_base);
	}
	constexpr typename types::reference       operator[](index i) &&{ assert( this->extension().contains(i) );
		typename types::element_ptr new_base = typename types::element_ptr(this->base()) + std::ptrdiff_t{Layout::operator()(i)};
		return typename types::reference(this->layout().sub_, new_base);
	}
	constexpr typename types::reference       operator[](index i) &{ assert( this->extension().contains(i) );
		typename types::element_ptr new_base = typename types::element_ptr(this->base()) + std::ptrdiff_t{Layout::operator()(i)};
		return typename types::reference(this->layout().sub_, new_base);
	}
	template<class Tp = std::array<index, static_cast<std::size_t>(D)>, typename = std::enable_if_t<(std::tuple_size<std::decay_t<Tp>>{}>1)> >
	auto operator[](Tp&& t) const
	->decltype(operator[](std::get<0>(t))[detail::tuple_tail(t)]){
		return operator[](std::get<0>(t))[detail::tuple_tail(t)];}
	template<class Tp, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tp>>::value==1> >
	auto operator[](Tp&& t) const
	->decltype(operator[](std::get<0>(t))){
		return operator[](std::get<0>(t));}
	template<class Tp = std::tuple<>, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tp>>::value==0> >
	decltype(auto) operator[](Tp&&) const{return *this;}
	using typename types::index;
	basic_const_array sliced(index first, index last) const&{
		typename types::layout_t new_layout = *this;
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};
	}
	basic_array sliced(index first, index last) &{
		typename types::layout_t new_layout = *this;
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};
	}
	basic_array sliced(index first, index last) &&{return sliced(first, last);}

	basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = *this; 
		new_layout.stride_*=s;
		return {new_layout, types::base_};
	}
	basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}
	using index_range = typename basic_array::index_range;
	decltype(auto) range(index_range ir) &     {return sliced(ir.front(), ir.front() + ir.size());}
	decltype(auto) range(index_range ir) &&    {return range(ir);}
	decltype(auto) range(index_range ir) const&{return sliced(ir.front(), ir.front() + ir.size());}

	auto range(typename types::index_range const& ir, dimensionality_type n) const{
		return rotated(n).range(ir).rotated(-n);
	}
	decltype(auto) flattened()&&{
		multi::biiterator<std::decay_t<decltype(std::move(*this).begin())>> biit{std::move(*this).begin(), 0, size(*(std::move(*this).begin()))};
		return basic_array<typename std::iterator_traits<decltype(biit)>::value_type, 1, decltype(biit)>{
			multi::layout_t<1>(1, 0, this->size()*size(*(std::move(*this).begin()))),
			biit
		};
	}
	friend decltype(auto) flattened(basic_array&& self){return std::move(self).flattened();}
	bool is_flattable() const{return this->stride() == this->layout().sub_.nelems_;}
	auto flatted() const{
		assert(is_flattable() && "flatted doesn't work for all layouts!");//this->nelems());
		multi::layout_t<D-1> new_layout{this->layout().sub_};
		new_layout.nelems_*=this->size();
		return basic_array<T, D-1, ElementPtr>{new_layout, types::base_};
	}
	friend auto flatted(basic_array const& self){return self.flatted();}
	template<typename Size>
	auto partitioned(Size const& s) const{
		assert(s!=0);
		assert(this->layout().nelems_%s==0);
		multi::layout_t<D+1> new_layout{this->layout(), this->layout().nelems_/s, 0, this->layout().nelems_};
		new_layout.sub_.nelems_/=s;
		return basic_array<T, D+1, ElementPtr>{new_layout, types::base_};
	}
	decltype(auto) transposed() const&{
		typename types::layout_t new_layout = *this;
		new_layout.transpose();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	friend decltype(auto) transposed(basic_array const& self){return self.transposed();}
	friend decltype(auto) operator~(basic_array const& self){return self.transposed();}

	
	decltype(auto) rotated()&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	decltype(auto) rotated()&&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		return basic_array{new_layout, types::base_};
	}
	decltype(auto) rotated() const&{
		typename types::layout_t new_layout = *this; new_layout.rotate();
		return basic_const_array{new_layout, types::base_};
	}
	friend decltype(auto) rotated(basic_array const&  self){return self.rotated();}
	friend decltype(auto) rotated(basic_array      && self){return std::move(self).rotated();}
	friend decltype(auto) rotated(basic_array      &  self){return self.rotated();}

	auto unrotated() &{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	auto unrotated() &&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_array<T, D, ElementPtr>{new_layout, types::base_};
	}
	auto unrotated() const&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate();
		return basic_const_array{new_layout, types::base_};
	}
	friend auto unrotated(basic_array const& self){return self.unrotated();}

	basic_array rotated(dimensionality_type i) &{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}
	basic_array rotated(dimensionality_type i) &&{return rotated(i);}
	basic_const_array rotated(dimensionality_type i) const&{
		typename types::layout_t new_layout = *this; 
		new_layout.rotate(i);
		return {new_layout, types::base_};
	}

	basic_array unrotated(dimensionality_type i) &{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}
	basic_array unrotated(dimensionality_type i) &&{return unrotated(i);}
	basic_const_array unrotated(dimensionality_type i) const&{
		typename types::layout_t new_layout = *this; 
		new_layout.unrotate(i);
		return {new_layout, types::base_};
	}

	decltype(auto) operator<<(dimensionality_type i) &{return rotated(i);}
	decltype(auto) operator>>(dimensionality_type i) &{return unrotated(i);}
	decltype(auto) operator<<(dimensionality_type i) &&{return std::move(*this).rotated(i);}
	decltype(auto) operator>>(dimensionality_type i) &&{return std::move(*this).unrotated(i);}
	decltype(auto) operator<<(dimensionality_type i) const&{return rotated(i);}
	decltype(auto) operator>>(dimensionality_type i) const&{return unrotated(i);}

	basic_array       operator()() &     {return *this;}
	basic_array       operator()() &&    {return operator()();}
	basic_const_array operator()() const&{return {this->layout(), this->base()};}

//	decltype(auto) operator()(index_range a) &     {return range(a);}
//	decltype(auto) operator()(index_range a) &&    {return std::move(*this).range(a);}
//	decltype(auto) operator()(index_range a) const&{return range(a);}
public:
	template<typename, dimensionality_type, typename, class> friend struct basic_array;
	basic_array       paren() &     {return *this;}
	basic_array       paren() &&    {return std::move(*this).operator()();}
	basic_const_array paren() const&{return {this->layout(), this->base()};}

	template<class... As> auto paren(index_range a, As... as) &     {return                  range(a).rotated().paren(as...).unrotated();}
	template<class... As> auto paren(index_range a, As... as) &&    {return std::move(*this).range(a).rotated().paren(as...).unrotated();}
	template<class... As> auto paren(index_range a, As... as) const&{return                  range(a).rotated().paren(as...).unrotated();}

	template<class... As> decltype(auto) paren(intersecting_range<index> inr, As... as) &     {return paren(intersection(this->extension(), inr), as...);}
	template<class... As> decltype(auto) paren(intersecting_range<index> inr, As... as) &&    {return std::move(*this).paren(intersection(this->extension(), inr), as...);}
	template<class... As> decltype(auto) paren(intersecting_range<index> inr, As... as) const&{return paren(intersection(this->extension(), inr), as...);}

	template<class... As> decltype(auto) paren(index i, As... as) &     {return                  operator[](i).paren(as...);}
	template<class... As> decltype(auto) paren(index i, As... as) &&    {return std::move(*this).operator[](i).paren(as...);}
	template<class... As> decltype(auto) paren(index i, As... as) const&{return                  operator[](i).paren(as...);}
public:
// default parameter here helps to determine irange type for {first, last} syntax
//	template<class... As>                                                                             decltype(auto) operator()(As... as                            )const&{return paren(as...                );}
//	template<class... As>                                                                             decltype(auto) operator()(As... as                            )    &&{return std::move(*this).paren(as...                );}
//	template<class... As>                                                                             decltype(auto) operator()(As... as                            )     &{return paren(as...                );}
//	template<class A1 = irange, class... As>                                                          decltype(auto) operator()(A1 a1, As... as                     )const&{return paren(a1, as...            );}
//	template<class A1 = irange, class A2 = irange, class... As>                                       decltype(auto) operator()(A1 a1, A2 a2, As... as              )const&{return paren(a1, a2, as...        );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class... As>                    decltype(auto) operator()(A1 a1, A2 a2, A3 a3, As... as       )const&{return paren(a1, a2, a3, as...    );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> decltype(auto) operator()(A1 a1, A2 a2, A3 a3, A4 a4, As... as)const&{return paren(a1, a2, a3, a4, as...);}
//	template<class... As>                                                                             decltype(auto) operator()(As... as                            )     &{return paren(as...                );}
//	template<class A1 = irange, class... As>                                                          decltype(auto) operator()(A1 a1, As... as                     )     &{return paren(a1, as...            );}
//	template<class A1 = irange, class A2 = irange, class... As>                                       decltype(auto) operator()(A1 a1, A2 a2, As... as              )     &{return paren(a1, a2, as...        );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class... As>                    decltype(auto) operator()(A1 a1, A2 a2, A3 a3, As... as       )     &{return paren(a1, a2, a3, as...    );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> decltype(auto) operator()(A1 a1, A2 a2, A3 a3, A4 a4, As... as)     &{return paren(a1, a2, a3, a4, as...);}
//	template<class... As>                                                                             decltype(auto) operator()(As... as                            )    &&{return std::move(*this).paren(as...                );}
//	template<class A1 = irange, class... As>                                                          decltype(auto) operator()(A1 a1, As... as                     )    &&{return std::move(*this).paren(a1, as...            );}
//	template<class A1 = irange, class A2 = irange, class... As>                                       decltype(auto) operator()(A1 a1, A2 a2, As... as              )    &&{return std::move(*this).paren(a1, a2, as...        );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class... As>                    decltype(auto) operator()(A1 a1, A2 a2, A3 a3, As... as       )    &&{return std::move(*this).paren(a1, a2, a3, as...    );}
//	template<class A1 = irange, class A2 = irange, class A3 = irange, class A4 = irange, class... As> decltype(auto) operator()(A1 a1, A2 a2, A3 a3, A4 a4, As... as)    &&{return std::move(*this).paren(a1, a2, a3, a4, as...);}
//	template<class A1 = irange, class... As>                    decltype(auto) operator()(A1 a1, irange a2, As... as       )const&{return paren(a1, a2, as...        );}
//	template<class A1 = irange, class... As>                    decltype(auto) operator()(A1 a1, irange a2, As... as       )     &{return paren(a1, a2, as...        );}
//	template<class A1 = irange, class... As>                    decltype(auto) operator()(A1 a1, irange a2, As... as       )    &&{return paren(a1, a2, as...        );}
//	template<class... As>                                       decltype(auto) operator()(irange r1, As... as              )const&{return paren(r1, as...        );}
//	template<class... As>                                       decltype(auto) operator()(irange r1, As... as              )     &{return paren(r1, as...        );}
//	template<class... As>                                       decltype(auto) operator()(irange r1, As... as              )    &&{return std::move(*this).paren(r1, as...        );}
/*
	template<class... As>                                       decltype(auto) operator()(index i1, irange r2, As... as              )const&{return paren(i1, r2, as...        );}
	template<class... As>                                       decltype(auto) operator()(index i1, irange r2, As... as              )     &{return paren(i1, r2, as...        );}
	template<class... As>                                       decltype(auto) operator()(index i1, irange r2, As... as              )    &&{return std::move(*this).paren(i1, r2, as...        );}
	template<class... As>                                       decltype(auto) operator()(irange r1, index i2, As... as              )const&{return paren(r1, i2, as...        );}
	template<class... As>                                       decltype(auto) operator()(irange r1, index i2, As... as              )     &{return paren(r1, i2, as...        );}
	template<class... As>                                       decltype(auto) operator()(irange r1, index i2, As... as              )    &&{return std::move(*this).paren(r1, i2, as...        );}
*/
#if 0
	auto operator()(index i) const&->decltype(operator[](i)){return operator[](i);}
	auto operator()(index i) &     ->decltype(operator[](i)){return operator[](i);}
	auto operator()(index i) &&    ->decltype(std::move(*this).operator[](i)){return std::move(*this).operator[](i);}

	template<class... As> auto operator()(irange r1)const&{return paren(r1);}
	template<class... As> auto operator()(irange r1)     &{return paren(r1);}
	template<class... As> auto operator()(irange r1)    &&{return std::move(*this).paren(r1);}

	auto operator()(irange r1, irange r2) const&{return paren(r1, r2);}
	auto operator()(irange r1, irange r2)      &{return paren(r1, r2);}
	auto operator()(irange r1, irange r2)     &&{return std::move(*this).paren(r1, r2);}

	decltype(auto) operator()(index r1, irange b2) const&{return paren(r1, b2);}
	decltype(auto) operator()(index r1, irange b2)     &{return paren(r1, b2);}
	decltype(auto) operator()(index r1, irange b2)    &&{return std::move(*this).paren(r1, b2);}

	auto operator()(irange r1, index b2) const&{return paren(r1, b2);}
	auto operator()(irange r1, index b2)      &{return paren(r1, b2);}
	auto operator()(irange r1, index b2)     &&{return std::move(*this).paren(r1, b2);}

	decltype(auto) operator()(index r1, index b2) const&{return paren(r1, b2);}
	decltype(auto) operator()(index r1, index b2)     &{return paren(r1, b2);}
	decltype(auto) operator()(index r1, index b2)    &&{return std::move(*this).paren(r1, b2);}

//	decltype(auto) operator()(index r1, index r2, index r3, index r4)const&{return paren(r1, r2, r3, r4);}
//	decltype(auto) operator()(index r1, index r2, index r3, index r4)    &&{return paren(r1, r2, r3, r4);}
//	decltype(auto) operator()(index r1, index r2, index r3, index r4)     &{return paren(r1, r2, r3, r4);}
#endif
/*
	template<class... As> auto operator()(irange r1, irange r2, irange r3, irange r4, As... as              )const&{return paren(r1, r2, r3, r4, as...        );}
	template<class... As> auto operator()(irange r1, irange r2, irange r3, irange r4, As... as              )     &{return paren(r1, r2, r3, r4, as...        );}
	template<class... As> auto operator()(irange r1, irange r2, irange r3, irange r4, As... as              )    &&{return std::move(*this).paren(r1, r2, r3, r4, as...);}

	template<class B1, class... As> decltype(auto) operator()(B1 b1, irange r2, irange r3, irange r4, As... as)const&{return paren(b1, r2, r3, r4, as...        );}
	template<class B2, class... As> decltype(auto) operator()(irange r1, B2 b2, irange r3, irange r4, As... as)const&{return paren(r1, b2, r3, r4, as...        );}
	template<class B3, class... As> decltype(auto) operator()(irange r1, irange r2, B3 b3, irange r4, As... as)const&{return paren(r1, r2, b3, r4, as...        );}
	template<class B4, class... As> decltype(auto) operator()(irange r1, irange r2, irange r3, B4 b4, As... as)const&{return paren(r1, r2, r3, b4, as...        );}

	template<class B1, class B2, class... As> decltype(auto) operator()(B1 b1, B2 b2, irange r3, irange r4, As... as)const&{return paren(b1, b2, r3, r4, as...        );}
	template<class B2, class B3, class... As> decltype(auto) operator()(irange r1, B2 b2, B3 b3, irange r4, As... as)const&{return paren(r1, b2, b3, r4, as...        );}
	template<class B3, class B4, class... As> decltype(auto) operator()(irange r1, irange r2, B3 b3, B4 b4, As... as)const&{return paren(r1, r2, b3, b4, as...        );}
*/

// the default template parameters below help interpret for {first, last} simple syntax as iranges
	template<class B1 = irange>                                                                       decltype(auto) operator()(B1 b1)                                const&{return paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         const&{return paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  const&{return paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) const&{return paren(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       decltype(auto) operator()(B1 b1)                                &{return paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         &{return paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  &{return paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) &{return paren(b1, b2, b3, b4, as...);}

	template<class B1 = irange>                                                                       decltype(auto) operator()(B1 b1)                                &&{return std::move(*this).paren(b1);}
	template<class B1 = irange, class B2 = irange>                                                    decltype(auto) operator()(B1 b1, B2 b2)                         &&{return std::move(*this).paren(b1, b2);}
	template<class B1 = irange, class B2 = irange, class B3 = irange>                                 decltype(auto) operator()(B1 b1, B2 b2, B3 b3)                  &&{return std::move(*this).paren(b1, b2, b3);}
	template<class B1 = irange, class B2 = irange, class B3 = irange, class B4 = irange, class... As> decltype(auto) operator()(B1 b1, B2 b2, B3 b3, B4 b4, As... as) &&{return std::move(*this).paren(b1, b2, b3, b4, as...);}

private:
	using Layout::nelems_;
	using Layout::stride_;
	using Layout::sub_;
public:
	using iterator = array_iterator<typename types::element, D, typename types::element_ptr, typename types::reference>;
	using const_iterator = array_iterator<typename types::element, D, typename types::element_const_ptr, typename types::const_reference>;
private:
	template<class Iterator>
	struct basic_reverse_iterator : 
		std::reverse_iterator<Iterator>,
		boost::multi::totally_ordered2<basic_reverse_iterator<Iterator>, void>
	{
		template<class O, typename = decltype(std::reverse_iterator<Iterator>{base(std::declval<O const&>())})>
		basic_reverse_iterator(O const& o) : std::reverse_iterator<Iterator>{base(o)}{}
		basic_reverse_iterator() : std::reverse_iterator<Iterator>{}{}
		explicit basic_reverse_iterator(Iterator it) : std::reverse_iterator<Iterator>(std::prev(it)){}
		explicit operator Iterator() const{auto ret = this->base(); if(ret!=Iterator{}) return ++ret; else return Iterator{};}
		explicit operator bool() const{return bool(this->base());}
		bool operator==(basic_reverse_iterator const& other) const{return (this->base() == other.base());}
		typename Iterator::reference operator*() const{return this->current;}
		typename Iterator::pointer operator->() const{return &this->current;}
		typename Iterator::reference operator[](typename Iterator::difference_type n) const{return *(this->current - n);}
		bool operator<(basic_reverse_iterator const& o) const{return o.base() < this->base();}
	};
public:
	using reverse_iterator = basic_reverse_iterator<iterator>;
	using ptr = basic_array_ptr<basic_array, Layout>;
	ptr operator&() &&{return {this->base_, this->layout()};}
	constexpr iterator begin(dimensionality_type d) &&{
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(0       ), l.sub_, l.stride_};
	}
	constexpr iterator end(dimensionality_type d) &&{
		Layout l = static_cast<Layout const&>(*this); l.rotate(d);
		return {types::base_ + l(l.size()), l.sub_, l.stride_};
	}

	constexpr iterator  begin()&{return {types::base_          , sub_, stride_};}
	constexpr iterator  end  ()&{return {types::base_ + nelems_, sub_, stride_};}
	friend iterator begin(basic_array& self){return self.begin();}
	friend iterator end  (basic_array& self){return self.end()  ;}

	constexpr iterator  begin()&&{return begin();}
	constexpr iterator  end  ()&&{return end();}
	friend iterator begin(basic_array&& self){return std::move(self).begin();}
	friend iterator end  (basic_array&& self){return std::move(self).end()  ;}

	constexpr const_iterator  begin() const&{return {types::base_          , sub_, stride_};}
	constexpr const_iterator  end  () const&{return {types::base_ + nelems_, sub_, stride_};}
	friend const_iterator begin(basic_array const& self){return self.begin();}
	friend const_iterator end  (basic_array const& self){return self.end()  ;}

protected:
	template<class A> void intersection_assign_(A&& other)&{// using multi::extension
		for(auto i : intersection(types::extension(), multi::extension(other)))
			operator[](i).intersection_assign_(std::forward<A>(other)[i]);
	}
	template<class A> void intersection_assign_(A&& o)&&{intersection_assign_(std::forward<A>(o));}
public:
	template<class It> void assign(It first, It last) &&{assert( this->size() == std::distance(first, last) );
		adl_copy(first, last, std::move(*this).begin());
	}
	template<class It> void assign(It first, It last) & {assert( this->size() == std::distance(first, last) );
		adl_copy(first, last, std::move(*this).begin());
	}
	template<class It> It assign(It first)&&{
		return adl_copy_n(first, this->size(), std::move(*this).begin()), first + std::move(*this).size();
	}
//	template<class It> void assign(It first) &&
//	->decltype(adl::copy_n(first, this->size(), this->begin())){
//		return adl::copy_n(first, this->size(), this->begin()), ;}
	template<class Range> auto assign(Range&& r) &&
	->decltype(this->assign(adl_begin(std::forward<Range>(r)), adl_end(std::forward<Range>(r)))){
		return this->assign(adl_begin(std::forward<Range>(r)), adl_end(std::forward<Range>(r)));}

	template<class Range> auto assign(Range&& r) &
	->decltype(this->assign(adl_begin(r), adl_end(r))){
		return this->assign(adl_begin(r), adl_end(r));}

	void assign(std::initializer_list<typename basic_array::value_type> il) const{assert( il.size() == this->size() );
		assign(il.begin(), il.end());
	}


	template<class A>//, typename = std::enable_if_t<not std::is_same<basic_array, std::decay_t<A>>{}>>
	basic_array& operator=(A&& o)&{
	//	assert(extension(*this) == extension(o));
		assert(this->extension() == o.extension());
		this->assign(adl_begin(std::forward<A>(o)), adl_end(std::forward<A>(o)));
		return *this;
	}
	template<class A>//, typename = std::enable_if_t<not std::is_same<basic_array, std::decay_t<A>>{}>>
	basic_array&& operator=(A&& o)&&{
	//	assert(extension(*this) == extension(o));
	//	assert(this->extension() == o.extension());
	//	std::move(*this).assign(adl::begin(std::forward<A>(o)), adl::end(std::forward<A>(o)));
		return std::move(this->operator=(std::forward<A>(o)));
	}
	basic_array&& operator=(basic_array&& o) &&{
		assert( this->extension() == o.extension() );
		std::move(*this).assign(std::move(o).begin(), std::move(o).end() );
		return std::move(*this);
	}
	template<class Array> void swap(Array&& o) &&{assert( std::move(*this).extension() == std::forward<Array>(o).extension() );
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> void swap(A&& o) &{return swap(std::forward<A>(o));}

	friend void swap(basic_array&& a, basic_array&& b){std::move(a).swap(std::move(b));}
	template<class Array> void swap(basic_array const& s, Array&& a){s.swap(a);}
	template<class Array> void swap(Array&& a, basic_array const& s){s.swap(a);}
	template<class Array>//, std::enable_if_t<std::is_same<Array, basic_array>{}, int> =0> 
	bool operator==(Array const& o) const&{
		return (this->extension()==o.extension()) and adl_equal(this->begin(), this->end(), adl_begin(o));
	}
//	bool operator==(basic_array const& o) = delete;//&&{return operator==<basic_array>(o);}
	template<class It>
	bool equal(It begin) const&{
		return adl_equal(
			std::move(modify(*this)).begin(), 
			std::move(modify(*this)).end(),
			begin
		);
	}
private:
	friend bool lexicographical_compare(basic_array&& a1, basic_array&& a2){
		if(a1.extension().first() > a2.extension().first()) return true;
		if(a1.extension().first() < a2.extension().first()) return false;
		return adl_lexicographical_compare(
			std::move(a1).begin(), std::move(a1).end(), 
			std::move(a2).begin(), std::move(a2).end()
		);
	}
public:
	template<class O> bool operator<(O&& o)&&{return lexicographical_compare(std::move(*this), std::move(o));}
	template<class O> bool operator>(O&& o)&&{return lexicographical_compare(std::move(o), std::move(*this));}
public:
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	constexpr basic_array<T2, D, P2> static_array_cast() const{
		P2 p2{this->base_};
		return basic_array<T2, D, P2>{this->layout(), p2};
	}
//	template<class T2, class P2 = decltype(boost::static_pointer_cast<T2>(std::declval<typename basic_array::element_ptr>()))>
//	auto static_array_cast() const HD ->basic_array<T2, D, P2>{
//		return basic_array<T2, D, P2>{this->layout(), boost::static_pointer_cast<T2>(this->base_)};
//	}
//	template<class T2>
//	auto static_array_cast() const HD -> basic_array<T2, D, decltype(boost::static_pointer_cast<T2>(std::declval<typename basic_array::element_ptr>()))>{
//		return {this->layout(), boost::static_pointer_cast<T2>(this->base_)};
//	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>,
		class Element = typename basic_array::element,
		class PM = T2 Element::*
	>
	basic_array<T2, D, P2> member_cast(PM pm) const{
		static_assert(sizeof(T)%sizeof(T2) == 0, 
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");
	//	return {this->layout().scale(sizeof(T)/sizeof(T2)), &(this->base_->*pm)};
		return basic_array<T2, D, P2>{this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}
	template<class T2, class P2 = T2*>//typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2> >
	basic_array<std::decay_t<T2>, D, P2> reinterpret_array_cast() const{
		static_assert( sizeof(T)%sizeof(T2)== 0, "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
		auto thisbase = this->base();
		return {
			this->layout().scale(sizeof(T)/sizeof(T2)), 
			static_cast<P2>(static_cast<void*>(thisbase))
		};
	}
	template<class T2, class P2 = T2*>
	basic_array<std::decay_t<T2>, D, P2> const_array_cast()&&{
		return {this->layout(), const_cast<P2>(this->base())};
	}
//	template<class T2, class P2 = T2*> // TODO implement move pointer
//	basic_array<std::decay_t<T2>, D, P2> move_array_cast()&&{
//		return {this->layout(), multi::make_const_iterator(this->base())};
//	}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class To, class From, std::enable_if_t<std::is_convertible<From, To>{},int> =0>
To _implicit_cast(From&& f){return static_cast<To>(f);}

template<class To, class From, std::enable_if_t<std::is_constructible<To, From>{} and not std::is_convertible<From, To>{},int> =0>
To _explicit_cast(From&& f){return static_cast<To>(f);}

template<class Element, typename Ptr, typename Ref>
struct array_iterator<Element, 1, Ptr, Ref> : 
	boost::multi::iterator_facade<
		array_iterator<Element, 1, Ptr, Ref>, 
		Element, std::random_access_iterator_tag, 
		Ref, multi::difference_type
	>,
	multi::affine<array_iterator<Element, 1, Ptr, Ref>, multi::difference_type>,
	multi::decrementable<array_iterator<Element, 1, Ptr, Ref>>,
	multi::incrementable<array_iterator<Element, 1, Ptr, Ref>>,
	multi::totally_ordered2<array_iterator<Element, 1, Ptr, Ref>, void>
{
	using affine = multi::affine<array_iterator<Element, 1, Ptr, Ref>, multi::difference_type>;
	using difference_type = typename affine::difference_type;

	array_iterator() = default;
	array_iterator(array_iterator const& other) = default;
	template<class Other, typename = decltype(_implicit_cast<Ptr>(typename Other::pointer{}))> 
	constexpr array_iterator(Other const& o) : data_{o.data_}, stride_{o.stride_}{}
	template<class Other, typename = decltype(_explicit_cast<Ptr>(typename Other::pointer{}))> 
	explicit constexpr array_iterator(Other const& o, int = 0) : data_{o.data_}, stride_{o.stride_}{}

	template<class EE, dimensionality_type, class PP, class RR> friend struct array_iterator;
	constexpr array_iterator(std::nullptr_t nu) HD : data_{nu}, stride_{1}{}
	constexpr array_iterator(Ptr const& p) HD : data_{p}, stride_{1}{}
	template<class EElement, typename PPtr, typename RRef>
	constexpr array_iterator(array_iterator<EElement, 1, PPtr, RRef> other) : data_{other.data_}, stride_{other.stride_}{} 
	explicit constexpr operator bool() const{return static_cast<bool>(this->data_);}
	constexpr Ref operator[](typename array_iterator::difference_type n) const{return *((*this) + n);}
	constexpr Ptr operator->() const{return data_;}
	using element = Element;
	using element_ptr = Ptr;
	using pointer = element_ptr;
	using rank = std::integral_constant<dimensionality_type, 1>;
	bool operator<(array_iterator const& o) const{return distance_to(o) > 0;}
private:
	constexpr array_iterator(Ptr d, typename basic_array<Element, 1, Ptr>::index s) : data_{d}, stride_{s}{}
	friend struct basic_array<Element, 1, Ptr>;
	Ptr data_ = nullptr;
	multi::index stride_;
	Ref dereference() const{return *data_;}
//	bool equal(array_iterator const& o) const{
//		assert(stride_ == o.stride_);
//		return data_==o.data_;// and stride_==o.stride_;
//	}
//	void increment(){data_ += stride_;}
//	void decrement(){data_ -= stride_;}
//	void advance(typename array_iterator::difference_type n) HD{data_ += stride_*n;}
	constexpr difference_type distance_to(array_iterator const& other) const{
		assert(stride_==other.stride_ and (other.data_-data_)%stride_ == 0);
		return (other.data_ - data_)/stride_;
	}
	auto base() const{return data_;}
	friend auto base(array_iterator const& self){return self.base();}
public:
	auto data() const HD{return data_;}
	auto stride() const{return stride_;}
	friend auto stride(array_iterator const& self){return self.stride();}
	array_iterator& operator++(){data_+=stride_; /*increment()*/; return *this;}
	array_iterator& operator--(){data_-=stride_; /*decrement()*/; return *this;}
	bool operator==(array_iterator const& o) const{return data_== o.data_;/*return equal(o);*/}
	bool operator!=(array_iterator const& o) const{return data_!= o.data_;/*return equal(o);*/}
	Ref operator*() const{return *data_;}//dereference();}
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
	using element_ref = typename std::iterator_traits<typename basic_array::element_ptr>::reference;//decltype(*typename basic_array::element_ptr{});
	decltype(auto) operator=(typename basic_array::element_type const& e) const&{
		adl_copy_n(&e, 1, this->base_); return *this;
	}
	bool operator==(typename basic_array::element const& e) const&{
		return adl_equal(&e, &e + 1, this->base_);
	}
	bool operator!=(typename basic_array::element const& e) const&{return not((*this)==e);}
	bool operator==(basic_array const& o) const&{
		return adl_equal(o.base_, o.base_ + 1, this->base_);
	}
	bool operator!=(basic_array const& o) const&{return not operator==(o);}
//	bool operator==(typename basic_array::element_type const& e) const&{
//		using std::equal; return equal(&e, &e + 1, this->base_);
//	}
//	bool operator!=(typename basic_array::element_type const& e) const&{return not operator==(e);}
//	operator element_ref() const{return *(this->base_);}
//	template<class TT> operator TT(){return static_cast<TT>(element_ref());}
	typename basic_array::element_ptr operator&() const{return this->base_;}
	using decay_type = typename types::element;
//	basic_array&
	element_ref operator()() const&{return *(this->base_);}
	operator element_ref()&&{return *(this->base_);}
//	decltype(auto) operator()() &&{return std::move(*(this->base_));}
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int){
		using boost::serialization::make_nvp;
		ar & make_nvp("element",  *(this->base_));
	//	std::for_each(this->begin(), this->end(), [&](auto&& e){ar & make_nvp("item", e);});
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
	auto get_allocator() const{return default_allocator_of(basic_array::base());}
	friend auto get_allocator(basic_array const& self){return self.get_allocator();}
	using decay_type = array<typename types::element, dimensionality_type{1}, decltype(default_allocator_of(std::declval<ElementPtr>()))>;
	       decay_type decay(          )&&      {return decay_type{std::move(*this)};}
	friend decay_type decay(basic_array&& self){return std::move(self).decay();}
	using basic_const_array = basic_array<
		T, 1, 
		typename std::pointer_traits<ElementPtr>::template rebind<typename basic_array::element_type const>,
	//	typename multi::iterator_traits<ElementPtr>::rebind_const, 
		Layout
	>;
protected:
	template<class A>
	void intersection_assign_(A&& other)&{
		for(auto idx : intersection(types::extension(), extension(other)))
			operator[](idx) = std::forward<A>(other)[idx];
	}
	template<class A> void intersection_assign_(A&& o)&&{intersection_assign_(std::forward<A>(o));}
protected:
	template<class TT, dimensionality_type DD, typename EP, class LLayout> friend struct basic_array;
	template<class TT, dimensionality_type DD, class Alloc> friend struct static_array;
	basic_array(basic_array const&) = default;
	template<class T2, class P2, class TT, dimensionality_type DD, class PP>
	friend decltype(auto) static_array_cast(basic_array<TT, DD, PP> const&);
public:
//	template<class Archive>
//	void serialize(Archive& ar, unsigned int){
//		for(auto&& e : *this) ar & BOOST_SERIALIZATION_NVP(e);
//	}
#if __cplusplus >= 201703L
protected: basic_array(basic_array&&) = default; // if you need to generate a copy you can't use `auto` here, use `decay` or `aut&&`.
#else
public   : basic_array(basic_array&&) = default; // in C++ < 17 this is necessary to return references from functions
// in c++17 things changed and non-moveable non-copyable types can be returned from functions and captured by auto
#endif
protected:
	template<class, class> friend struct basic_array_ptr;
	template<class, dimensionality_type D, class, class>
	friend struct array_iterator;
public:
	friend constexpr auto dimensionality(basic_array const& self){return self.dimensionality;}
	template<class BasicArray, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<BasicArray>>{}>, typename = decltype(types(std::declval<BasicArray&&>()))> 
	constexpr basic_array(BasicArray&& other) : types{std::forward<BasicArray>(other)}{}
	basic_array_ptr<basic_array, Layout> operator&() const{
		return {this->base_, this->layout()};
	}
	void assign(std::initializer_list<typename basic_array::value_type> il) const{assert( il.size() == static_cast<std::size_t>(this->size()) );
		assign(il.begin(), il.end());
	}
	template<class It> 
	basic_array&  assign(It first, It last) &{assert( std::distance(first, last) == this->size() );
		return adl_copy(first, last, this->begin()), *this;
	}
	template<class It> 
	basic_array&& assign(It first, It last)&&{return std::move(assign(first, last));}
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int){
		using boost::serialization::make_nvp;
		std::for_each(std::move(*this).begin(), std::move(*this).end(),[&](auto&& e){ar&make_nvp("item",e);});
	}
	template<class A, 
		typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}>,
		typename = decltype(
			std::declval<typename basic_array::reference&&>() 
				= std::declval<typename multi::array_traits<typename std::remove_reference_t<A>>::reference&&>()
		)
	>
	basic_array&& operator=(A&& o)&&{
		assert(this->extension() == extension(o));
		return this->assign(adl_begin(std::forward<A>(o)), adl_end(std::forward<A>(o)));
	}
	basic_array&& operator=(basic_array const& o)&&{
		return std::move(*this).assign(std::move(modify(o)).begin()), std::move(*this);
	}
	basic_array&& operator=(basic_array&& o)&&{ assert(this->extension() == o.extension());
		return std::move(*this).assign(std::move(o).begin()), std::move(*this);
	}
	template<class TT, dimensionality_type DD, class... As>
	basic_array&& operator=(basic_array<TT, DD, As...>&& o)&&{assert(this->extension() == o.extension());
		std::move(*this).assign(std::move(o).begin(), std::move(o).end()); return std::move(*this);
	}
	typename basic_array::const_reference operator[](index i) const&{assert( this->extension().contains(i) );
		return *(this->base() + Layout::operator()(i)); // in C++17 this is allowed even with syntethic references
	}
	typename basic_array::reference operator[](index i)&{assert(this->extension().contains(i));
		return *(this->base() + Layout::operator()(i));
	}
	typename basic_array::reference operator[](index i)&&{return this->operator[](i);}

	template<class Tuple, typename = std::enable_if_t<(std::tuple_size<std::decay_t<Tuple>>{}>1) > >
	constexpr auto operator[](Tuple&& t) const
	->decltype(operator[](std::get<0>(t))[detail::tuple_tail(t)]){
		return operator[](std::get<0>(t))[detail::tuple_tail(t)];}
	template<class Tuple, typename = std::enable_if_t<std::tuple_size<std::decay_t<Tuple>>{}==1> >
	decltype(auto) operator[](Tuple&& t) const{return operator[](std::get<0>(t));}
	decltype(auto) operator[](std::tuple<>) const{return *this;}

	using typename types::index;
	basic_array sliced(typename types::index first, typename types::index last)&{
		typename types::layout_t new_layout = *this; 
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};		
	}
	basic_array sliced(typename types::index first, typename types::index last)&&{return sliced(first, last);}
	basic_const_array sliced(typename types::index first, typename types::index last) const&{
		typename types::layout_t new_layout = *this; 
		(new_layout.nelems_/=Layout::size())*=(last - first);
		return {new_layout, types::base_ + Layout::operator()(first)};		
	}
	basic_array strided(typename types::index s) const{
		typename types::layout_t new_layout = this->layout();
		new_layout.stride_*=s;
		return {new_layout, types::base_};//+ Layout::operator()(this->extension().front())};
	}
	basic_array sliced(typename types::index first, typename types::index last, typename types::index stride) const{
		return sliced(first, last).strided(stride);
	}
	auto range(index_range const& ir) &{return sliced(ir.front(), ir.last());}
	auto range(index_range const& ir) &&{return sliced(ir.front(), ir.last());}
	auto range(index_range const& ir) const&{return sliced(ir.front(), ir.last());}


	decltype(auto) operator()()&&{return std::move(*this);}

	auto operator()(index_range const& ir) &{return range(ir);}
	auto operator()(index_range const& ir) &&{return std::move(*this).range(ir);}
	auto operator()(index_range const& ir) const&{return range(ir);}

	decltype(auto) operator()(index i) &     {return operator[](i);}
	decltype(auto) operator()(index i) &&    {return std::move(*this).operator[](i);}
	decltype(auto) operator()(index i) const&{return operator[](i);}

	auto paren() &{return operator()();}
	auto paren() &&{return std::move(*this).operator()();}
	auto paren() const&{return operator()();}

	auto paren(index_range const& ir) &{return range(ir);}
	auto paren(index_range const& ir) &&{return std::move(*this).range(ir);}
	auto paren(index_range const& ir) const&{return range(ir);}

	decltype(auto) paren(index i) &     {return operator[](i);}
	decltype(auto) paren(index i) &&    {return std::move(*this).operator[](i);}
	decltype(auto) paren(index i) const&{return operator[](i);}

	template<typename Size>
	auto partitioned(Size const& s) const{
		assert( this->layout().nelems_%s==0 );
		multi::layout_t<2> new_layout{this->layout(), this->layout().nelems_/s, 0, this->layout().nelems_};
		new_layout.sub_.nelems_/=s;
		return basic_array<T, 2, ElementPtr>{new_layout, types::base_};
	}
	friend decltype(auto) rotated(basic_array const& self){return self.rotated();}
	friend decltype(auto) unrotated(basic_array const& self){return self.unrotated();}
//	friend decltype(auto) transposed(basic_array const& self){return self.transposed();}
//	friend decltype(auto) operator~(basic_array const& self){return transposed(self);}

	decltype(auto) rotated(dimensionality_type = 1) &     {return operator()();}
	decltype(auto) rotated(dimensionality_type = 1) &&    {return std::move(*this).operator()();}
	decltype(auto) rotated(dimensionality_type = 1) const&{return operator()();}

	decltype(auto) unrotated(dimensionality_type = 1) &     {return operator()();}
	decltype(auto) unrotated(dimensionality_type = 1) &&    {return std::move(*this).operator()();}
	decltype(auto) unrotated(dimensionality_type = 1) const&{return operator()();}

	decltype(auto) operator<<(dimensionality_type i) const{return rotated(i);}
	decltype(auto) operator>>(dimensionality_type i) const{return unrotated(i);}

	using iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_ptr, typename types::reference>;
	using const_iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_const_ptr>;
	using reverse_iterator = std::reverse_iterator<iterator>;

//	using const_iterator = typename multi::array_iterator<typename types::element, 1, typename types::element_const_ptr, typename types::const_reference>;
//	using const_iterator = multi::const_iterator<iterator>;

	constexpr iterator begin()&           {return {this->base_, this->stride_};}
	constexpr iterator begin()&&          {return begin();}
	constexpr const_iterator begin()const&{return iterator{this->base_, this->stride_};}

	constexpr iterator end  ()&{return {basic_array::base_ + types::nelems_, basic_array::stride_};}
	constexpr iterator end  ()&&{return end();}
	constexpr const_iterator end  ()const&{return iterator{basic_array::base_ + types::nelems_, basic_array::stride_};}

	friend iterator begin(basic_array      & s){return s.begin();}
	friend iterator begin(basic_array     && s){return std::move(s).begin();}
	friend const_iterator begin(basic_array const& s){return s.begin();}

	friend iterator end(basic_array      & s){return s.end();}
	friend iterator end(basic_array     && s){return std::move(s).end();}
	friend const_iterator end(basic_array const& s){return s.end();}

//	constexpr const_iterator begin() const&{return {basic_array::base_                 , basic_array::stride_};}
//	constexpr const_iterator end  () const&{return {basic_array::base_ + types::nelems_, basic_array::stride_};}

	template<class It> auto assign(It f)&& //	->decltype(adl::copy_n(f, this->size(), begin(std::move(*this))), void()){
	->decltype(adl_copy_n(f, this->size(), std::declval<iterator>()), void()){
		return adl_copy_n(f, this->size(), std::move(*this).begin()), void();}

	template<typename Array>//, typename = std::enable_if_t<not std::is_base_of<basic_array, Array>{}> >
	constexpr bool operator==(Array const& o) const&{ // TODO assert extensions are equal?
		return (this->extension()==extension(o)) and adl_equal(this->begin(), this->end(), adl_begin(o));
	}
//	constexpr bool operator==(basic_array const& other) const&{
//		return (this->extension()==extension(std::move(modify(other))))
//			and adl_equal(std::move(*this).begin(), std::move(*this).end(), std::move(modify(other)).begin());
//		return this->operator==<basic_array>(std::move(other));
//	}
//	constexpr bool operator==(basic_array&& other)&&{
//		return (basic_array::extension()==extension(std::move(other))) 
//			and adl_equal(std::move(*this).begin(), std::move(*this).end(), adl::begin(std::move(other)));
//		return this->operator==<basic_array>(std::move(other));
//	}
	bool operator<(basic_array const& o) const&{return lexicographical_compare(*this, o);}//operator< <basic_array const&>(o);}
	template<class Array> void swap(Array&& o)&&{{using multi::extension; assert(this->extension() == extension(o));}
		adl_swap_ranges(this->begin(), this->end(), adl_begin(std::forward<Array>(o)));
	}
	template<class A> void swap(A&& o)&{return swap(std::forward<A>(o));}
	friend void swap(basic_array&& a, basic_array&& b){std::move(a).swap(std::move(b));}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend void swap(basic_array&& s, A&& a){s.swap(a);}
	template<class A, typename = std::enable_if_t<not std::is_base_of<basic_array, std::decay_t<A>>{}> > friend void swap(A&& a, basic_array&& s){s.swap(a);}
private:
	template<class A1, class A2>
	static auto lexicographical_compare(A1 const& a1, A2 const& a2){
		using multi::extension;
		if(extension(a1).first() > extension(a2).first()) return true;
		if(extension(a1).first() < extension(a2).first()) return false;
		return adl_lexicographical_compare(adl_begin(a1), adl_end(a1), adl_begin(a2), adl_end(a2));
	}
public:
	template<class O>
	bool operator<(O const& o) const{return lexicographical_compare(*this, o);}
	template<class O>
	bool operator>(O const& o) const{return lexicographical_compare(o, *this);}
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
		class PM = T2 Element::*
	>
	basic_array<T2, 1, P2> member_cast(PM pm) const{
		static_assert(sizeof(T)%sizeof(T2) == 0, 
			"array_member_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom alignas structures (to the interesting member(s) sizes) or custom pointers to allow reintrepreation of array elements");
		return {this->layout().scale(sizeof(T)/sizeof(T2)), static_cast<P2>(&(this->base_->*pm))};
	}
	template<class T2, class P2 = typename std::pointer_traits<typename basic_array::element_ptr>::template rebind<T2>>
	basic_array<T2, 1, P2> reinterpret_array_cast() const{
		static_assert( sizeof(T)%sizeof(T2)== 0, "error: reinterpret_array_cast is limited to integral stride values, therefore the element target size must be multiple of the source element size. Use custom pointers to allow reintrepreation of array elements in other cases" );
//			this->layout().scale(sizeof(T)/sizeof(T2));
		static_assert( sizeof(P2) == sizeof(typename basic_array::element_ptr), "reinterpret on equal size?");
		auto const thisbase = this->base();
		P2 new_base; std::memcpy(&new_base, &thisbase, sizeof(P2)); //reinterpret_cast<P2 const&>(thisbase) // TODO find a better way, fancy pointers wouldn't need reinterpret_cast
		return {this->layout().scale(sizeof(T)/sizeof(T2)), new_base};
	}
};

template<class T2, class P2, class Array, class... Args>
constexpr decltype(auto) static_array_cast(Array&& a, Args&&... args){return a.template static_array_cast<T2, P2>(std::forward<Args>(args)...);}

template<class T2, class Array, class P2 = typename std::pointer_traits<typename std::decay<Array>::type::element_ptr>::template rebind<T2> , class... Args>
constexpr decltype(auto) static_array_cast(Array&& a, Args&&... args){return a.template static_array_cast<T2, P2>(std::forward<Args>(args)...);}

template<
	class T2, class Array,
	class P2 = typename std::pointer_traits<typename std::decay<Array>::type::element_ptr>::template rebind<T2>
>
decltype(auto) reinterpret_array_cast(Array&& a){
	return a.template reinterpret_array_cast<T2, P2>();
}

template<class T2, class Array, class P2 = typename std::pointer_traits<typename std::decay_t<Array>::element_ptr>::template rebind<T2>,
	class PM = T2 std::decay_t<Array>::element::*
>
decltype(auto) member_array_cast(Array&& a, PM pm){return a.template member_cast<T2, P2>(pm);}

template<
	class T2, class Array, class P2 = typename std::pointer_traits<typename std::decay_t<Array>::element_ptr>::template rebind<T2>,
	class OtherT2, class OtherElement, std::enable_if_t<not std::is_same<T2, OtherElement>{}, int> =0
>
decltype(auto) member_array_cast(Array&& a, OtherT2 OtherElement::* pm){
	static_assert(sizeof(OtherElement)==sizeof(typename std::decay_t<Array>::element_type),
		"member cast does not implicitly reinterprets between types of different size");
	static_assert(sizeof(OtherT2)==sizeof(T2), 
		"member cast does not implicitly reinterprets between types of different size");
	return reinterpret_array_cast<OtherElement>(std::forward<Array>(a)).template member_cast<T2>(pm);
}

template<typename T, dimensionality_type D, typename ElementPtr = T*>
struct array_ref : 
//TODO	multi::partially_ordered2<array_ref<T, D, ElementPtr>, void>,
	basic_array<T, D, ElementPtr>
{
protected:
	constexpr array_ref() noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{}, nullptr}{}
public:
	[[deprecated("references are not copyable, use &&")]]
	array_ref(array_ref const&) = default; // don't try to use `auto` for references, use `auto&&` or explicit value type
	constexpr array_ref(typename array_ref::element_ptr p, typename array_ref::extensions_type e = {}) noexcept
		: basic_array<T, D, ElementPtr>{typename array_ref::types::layout_t{e}, p}{}

	template<class TT, std::size_t N>
	constexpr array_ref(TT(&t)[N]) : basic_array<T, D, ElementPtr>{data_elements(t), extensions(t)}{}

//	constexpr array_ref(typename array_ref::element_ptr p, std::initializer_list<index_extension> il) noexcept
//		: array_ref(p, detail::to_tuple<D, index_extension>(il)){}
//	template<class Extension>//, typename = decltype(array_ref(std::array<Extension, D>{}, allocator_type{}, std::make_index_sequence<D>{}))>
//	constexpr array_ref(typename array_ref::element_ptr p, std::array<Extension, D> const& x) 
//		: array_ref(p, x, std::make_index_sequence<D>{}){}
//	using basic_array<T, D, ElementPtr>::operator[];
	using basic_array<T, D, ElementPtr>::operator=;
	using basic_array<T, D, ElementPtr>::operator==;
//	using basic_array<T, D, ElementPtr>::operator<;
//	using basic_array<T, D, ElementPtr>::operator>;
//	template<class ArrayRef> explicit array_ref(ArrayRef&& a) : array_ref(a.data(), extensions(a)){}
private:
	template<class It> auto copy_elements(It first){
		return adl_copy_n(first, array_ref::num_elements(), array_ref::data_elements());
	}
	template<class It> auto equal_elements(It first) const{
		return adl_equal(first, first + this->num_elements(), this->data_elements());
	}
	template<class TT, std::size_t N> using carr = TT[N];
public:
	template<class TT, std::size_t N>//, std::enable_if<std::is_convertible<typename array_ref::element_ptr, std::remove_all_extents_t<carr<TT, N>> >{}, int> =0> 
	operator carr<TT, N>&() const&{
		assert(extensions(*(carr<TT, N>*)this)==this->extensions());
		return *reinterpret_cast<carr<TT, N>*>(this->base_);
	}
	typename array_ref::element_ptr data_elements() const&{return array_ref::base_;}
	array_ref&& operator=(array_ref const& o) &&{assert(this->num_elements()==o.num_elements());
		return array_ref::copy_elements(o.data_elements()), std::move(*this);
	}
	template<typename TT, dimensionality_type DD = D, class... As>
	array_ref&& operator=(array_ref<TT, DD, As...> const& o)&&{assert(this->extensions() == o.extensions());
		return adl_copy_n(o.data(), o.num_elements(), this->data()), std::move(*this);
	}
	template<typename TT, dimensionality_type DD = D, class... As>
	bool operator==(array_ref<TT, DD, As...>&& o) const&{
		if( this->extensions() != o.extensions() ) return false; // TODO, or assert?
		return equal_elements(std::move(o).data_elements());
	}
	typename array_ref::element_ptr data_elements()&&{return array_ref::base_;}
	friend typename array_ref::element_ptr data_elements(array_ref&& s){return std::move(s).data_elements();}

	typename array_ref::element_ptr data() const& HD{return array_ref::base_;} 
	friend typename array_ref::element_ptr data(array_ref const& self){return self.data();}
	friend decltype(auto) operator*(array_ref const& self){
		return static_cast<typename array_ref::decay_type const&>(self);
	}
	
//	explicit 
//	operator typename array_ref::decay_type const&() const&{
//		return static_cast<typename array_ref::decay_type const&>(*this);		
// 	}
//	operator typename array_ref::decay_type() &&{
//		return static_cast<typename array_ref::decay_type&&>(std::move(*this));		
//	}
//	operator typename array_ref::decay_type&&() & = delete;
	typename array_ref::decay_type const& decay() const&{
		return static_cast<typename array_ref::decay_type const&>(*this);
	}
	friend typename array_ref::decay_type const& decay(array_ref const& s){return s.decay();}
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int v){
		using boost::serialization::make_nvp;
		if(this->num_elements() < (2<<8) ) std::move(*this).basic_array<T, D, ElementPtr>::serialize(ar, v);
		else{
			using boost::serialization::make_binary_object;
			using boost::serialization::make_array;
			if(std::is_trivially_copy_assignable<typename array_ref::element>{})
				ar & make_nvp("binary_data", make_binary_object(this->data(), sizeof(typename array_ref::element)*this->num_elements())); //#include<boost/serialization/binary_object.hpp>
			else ar & make_nvp("data", make_array(this->data(), this->num_elements()));
		}
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

//template<class T, dimensionality_type D, typename Ptr = T*>
//using array_ptr = basic_array_ptr<array_ref<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;

template<class T, dimensionality_type D, typename Ptr = T*>
struct array_ptr : basic_array_ptr<basic_array<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>{
	using basic_ptr = basic_array_ptr<basic_array<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;
//	using basic_ptr = basic_array_ptr<array_ref<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>;
//	using basic_ptr::basic_ptr;//array_ptr<array_ref<T, D, Ptr>, typename array_ref<T, D, Ptr>::layout_t>::basic_array_ptr;
public:
	array_ptr(Ptr p, index_extensions<D> x) : basic_ptr(p, multi::layout_t<D>{x}){}
	template<class TT, std::size_t N>
	constexpr array_ptr(TT(*t)[N]) : basic_ptr(data_elements(*t), layout(*t)){}
	array_ref<T, D, Ptr> operator*() const{
		return {this->base(), (*this)->extensions()};//multi::layout_t<D>{x}};
	}
};

template<class TT, std::size_t N>
// auto operator&(TT(&t)[N]){ // c++ cannot overload & for primitive types
auto addressof(TT(&t)[N]){
	return array_ptr<
		std::decay_t<std::remove_all_extents_t<TT[N]>>, std::rank<TT[N]>{}, std::remove_all_extents_t<TT[N]>*
	>(&t);
}

template<class T, dimensionality_type D, typename Ptr = T*>
using array_cptr = array_ptr<T, D, 	typename std::pointer_traits<Ptr>::template rebind<T const>>;

template<dimensionality_type D, class P>
array_ref<typename std::iterator_traits<P>::value_type, D, P> 
make_array_ref(P p, index_extensions<D> x){return {p, x};}

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

#if __cpp_deduction_guides
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
constexpr auto rotated(const T(&t)[N]) noexcept{
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

template<class RandomAccessIterator, dimensionality_type D>
multi::array_ptr<typename std::iterator_traits<RandomAccessIterator>::value_type, D, RandomAccessIterator>
operator/(RandomAccessIterator data, multi::iextensions<D> x){return {data, x};}

}}

#undef HD

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


#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MULTI_ARRAY_REF

#include<cassert>
#include<numeric> // iota
#include<iostream>
#include<vector>

using std::cout;
namespace multi = boost::multi;


int main(){

	{
		double a[4][5] = {{1.,2.},{2.,3.}};
		multi::array_ref<double, 2> A(&a[0][0], {4, 5});
		multi::array_ref<double, 2, double const*> B(&a[0][0], {4, 5});
		multi::array_ref<double const, 2> C(&a[0][0], {4, 5});
		multi::array_cref<double, 2> D(&a[0][0], {4, 5});
		A[1][1] = 2.;
//		A[1].cast<double const*>()[1] = 2.;
		double d[4][5] = {{1.,2.},{2.,3.}};
	//	typedef d45 = double const[4][5];
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
		{auto x = extensions(A2);
		for(auto i : std::get<0>(x) ){
			for(auto j : std::get<1>(x) ) cout<< A2[i][j] <<' ';
			cout<<'\n';
		}}
		auto const& A3 = A({0, 3}, 1, {0, 2});
		assert( multi::rank<std::decay_t<decltype(A3)>>{} == 2 and num_elements(A3) == 6 );
		{
			auto x = extensions(A3);
			for(auto i : std::get<0>(x)){
				for(auto j : std::get<1>(x)) cout<< A3[i][j] <<' ';
				cout<<'\n';
			}
		}
	}
	return 0;
}
#undef HD

#endif
#endif

