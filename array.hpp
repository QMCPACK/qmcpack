#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
echo $X
[ ! -d build.$X ] && (mkdir build.$X && cd build.$X && cmake ..)
cd build.$X && make -j && ctest -j --output-on-failure
exit
#endif
// $CXXX $CXXFLAGS $0 -o $0.$X&&$0.$X&&rm $0.$X;exit
//  © Alfredo A. Correa 2018-2019

#ifndef BOOST_MULTI_ARRAY_HPP 
#define BOOST_MULTI_ARRAY_HPP

#include "./array_ref.hpp"
#include "./config/NO_UNIQUE_ADDRESS.hpp"

#include "./memory/allocator.hpp"
#include "./detail/memory.hpp"
#include "./detail/adl.hpp"

#include<memory>

namespace boost{
namespace multi{

template<class Allocator> struct array_allocator{
	using allocator_type = Allocator;
protected:
	MULTI_NO_UNIQUE_ADDRESS allocator_type alloc_;
	allocator_type& alloc(){return alloc_;}
	array_allocator(allocator_type const& a = {}) : alloc_{a}{}
	typename std::allocator_traits<allocator_type>::pointer 
	allocate(typename std::allocator_traits<allocator_type>::size_type n){
		return n?std::allocator_traits<allocator_type>::allocate(alloc_, n):nullptr;
	}
	typename std::allocator_traits<allocator_type>::pointer 
	allocate(typename std::allocator_traits<allocator_type>::size_type n, typename std::allocator_traits<allocator_type>::const_void_pointer hint){
		return n?std::allocator_traits<allocator_type>::allocate(alloc_, n, hint):nullptr;
	}
	auto uninitialized_fill_n(typename std::allocator_traits<allocator_type>::pointer base, typename std::allocator_traits<allocator_type>::size_type num_elements, typename std::allocator_traits<allocator_type>::value_type e){
		return adl_alloc_uninitialized_fill_n(alloc_, base, num_elements, e);
	}
	template<typename It> 
	auto uninitialized_copy_n(It first, size_type n, typename std::allocator_traits<allocator_type>::pointer data){
		return adl_alloc_uninitialized_copy_n(alloc_, first, n, data);
	}
	template<typename It> 
	auto destroy_n(It first, size_type n){
		return adl_alloc_destroy_n(this->alloc(), first, n);
	}
public:
	allocator_type get_allocator() const{return alloc_;}
	friend allocator_type get_allocator(array_allocator const& s){return s.get_allocator();}
};

template<class T, class Ptr = T*> struct move_ptr : std::move_iterator<Ptr>{
	using std::move_iterator<Ptr>::move_iterator;
	explicit operator Ptr() const{return std::move_iterator<Ptr>::base();}
};

// static_array is not a value type because it doesn't define assignment for static_arrays of different extensions
template<class T, dimensionality_type D, class Alloc = std::allocator<T>>
struct static_array : 
	protected array_allocator<Alloc>,
	public array_ref<T, D, typename std::allocator_traits<typename array_allocator<Alloc>::allocator_type>::pointer>
{
private:
	using array_alloc = array_allocator<Alloc>;
public:	
	static_assert( std::is_same<typename std::allocator_traits<Alloc>::value_type, typename static_array::element>{}, 
		"allocator value type must match array value type");
	static_assert( std::is_same<typename std::allocator_traits<Alloc>::pointer, typename static_array::element_ptr>{}, 
		"allocator pointer type must match array pointer type");
	using array_alloc::get_allocator;
	using typename array_allocator<Alloc>::allocator_type;
//	using allocator_type = typename static_array::allocator_type;
	using decay_type = array<T, D, Alloc>;
protected:
	using alloc_traits = typename std::allocator_traits<typename static_array::allocator_type>;
	using ref = array_ref<T, D, typename std::allocator_traits<typename std::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;
	auto uninitialized_value_construct(){
		return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}
	auto uninitialized_default_construct(){
	//	return std::uninitialized_default_construct_n(this->base_, this->num_elements());
		return adl_alloc_uninitialized_default_construct_n(static_array::alloc(), this->base_, this->num_elements());
	}
	template<typename It> auto uninitialized_copy_elements(It first){
		return array_alloc::uninitialized_copy_n(first, this->num_elements(), this->data_elements());
	}
	void destroy_aux(std::false_type){array_alloc::destroy_n(this->data_elements(), this->num_elements());}
	void destroy_aux(std::true_type ){}
	void destroy(){destroy_aux(std::is_trivially_destructible<typename static_array::element>{});}
	void allocate(){this->base_ = array_alloc::allocate(static_array::num_elements());}
public:
	using value_type = typename std::conditional<
		(static_array::dimensionality>1),
		array<typename static_array::element, static_array::dimensionality-1, allocator_type>, 
		typename std::conditional<
			static_array::dimensionality == 1,
			typename static_array::element,
			typename static_array::element // TODO or void?
		>::type
	>::type;

	using typename ref::size_type;
	using typename ref::difference_type;
	explicit static_array(typename static_array::allocator_type const& a) : array_alloc{a}{}
protected:
	static_array(static_array&& other, typename static_array::allocator_type const& a) noexcept     //6b
	:	array_alloc{a},
		ref{other.base_, other.extensions()}
	{
		other.ref::layout_t::operator=({});
		other.base_ = nullptr;
	}
public:
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(
		basic_array<typename static_array::element, static_array::dimensionality, multi::move_ptr<typename static_array::element, typename static_array::element_ptr>>&& other, 
		typename static_array::allocator_type const& a = {}
	) noexcept : 
		array_alloc{a},
		ref{
			other.layout()==typename static_array::layout_t(other.extensions())?
				other.base_.base():
				array_alloc::allocate(other.num_elements())
			,
			other.extensions()
		}
	{
		if(other.base_.base() != static_array::base_)
			recursive<D>::alloc_uninitialized_copy(static_array::alloc(), 
				other.template static_array_cast<typename static_array::element, typename static_array::element_ptr>().begin(), 
				other.template static_array_cast<typename static_array::element, typename static_array::element_ptr>().end()  , 
				this->begin()
			);
	}
//	template<class Array>//, std::enable_if_t<std::is_same<Array, basic_array>{}, int> =0> 
//	auto operator==(Array&& o) const&
//	->decltype(std::move(modify(*this)).ref::operator==(std::forward<Array>(o))){
//		return std::move(modify(*this)).ref::operator==(std::forward<Array>(o));}

//	auto operator==(static_array const& o) const&{return std::move(modify(*this)).ref::operator==(std::move(modify(o)));}

	template<class TT, class... Args> bool operator==(basic_array<TT, D, Args...> const& other) const{return ref::operator==(other);}
	template<class TT, class... Args> bool operator!=(basic_array<TT, D, Args...> const& other) const{return ref::operator!=(other);}

	template<class It, class=typename std::iterator_traits<std::decay_t<It>>::difference_type>//edecltype(std::distance(std::declval<It>(), std::declval<It>()), *std::declval<It>())>      
	// analogous to std::vector::vector (5) https://en.cppreference.com/w/cpp/container/vector/vector
	static_array(It first, It last, typename static_array::allocator_type const& a = {}) : 
		array_alloc{a},
		ref{
			array_alloc::allocate(typename static_array::layout_t{index_extension(adl_distance(first, last))*multi::extensions(*first)}.num_elements()), 
			index_extension(adl_distance(first, last))*multi::extensions(*first)
		}
	{
//		recursive<D>::alloc_uninitialized_copy(static_array::alloc(), first, last, this->begin());
		adl_alloc_uninitialized_copy(static_array::alloc(), first, last, ref::begin());
//		adl::uninitialized_copy(first, last, ref::begin());
	}

	template<
		class Range, class=std::enable_if_t<not std::is_base_of<static_array, std::decay_t<Range>>{}>, 
		class=decltype(/*static_array*/(std::declval<Range&&>().begin(), std::declval<Range&&>().end())), // instantiation of static_array here gives a compiler error in 11.0
		class=std::enable_if_t<not is_basic_array<Range&&>{}>// TODO add is_assignable<value_type> check
	> 
	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions
	static_array(Range&& rng) : static_array(std::forward<Range>(rng).begin(), std::forward<Range>(rng).end()){}

	template<class TT> 
	auto uninitialized_fill_elements(TT const& value){
		return array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), value);
	}

	template<class TT, class... As> 
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(array_ref<TT, D, As...> const& other, typename static_array::allocator_type const& a = {}) :
		array_alloc{a},
		ref{array_alloc::allocate(other.num_elements()), other.extensions()}
	{
		adl_alloc_uninitialized_copy_n(static_array::alloc(), other.data_elements(), other.num_elements(), this->data_elements());
	}

	static_array(typename static_array::extensions_type x, typename static_array::element const& e, typename static_array::allocator_type const& a) : //2
		array_alloc{a}, 
		ref(array_alloc::allocate(typename static_array::layout_t{x}.num_elements()), x)
	{
		array_alloc::uninitialized_fill_n(this->data_elements(), this->num_elements(), e);
	}
	template<class Element, std::enable_if_t<std::is_convertible<Element, typename static_array::element>{} and D==0, int> = 0>
	explicit static_array(Element const& e, typename static_array::allocator_type const& a) :
		static_array(typename static_array::extensions_type{}, e, a){}

	static_array(typename static_array::extensions_type x, typename static_array::element const& e) : //2
		array_alloc{}, ref(array_alloc::allocate(typename static_array::layout_t{x}.num_elements()), x)
	{
		array_alloc::uninitialized_fill_n(this->base(), this->num_elements(), e);
	}
	explicit static_array(typename static_array::extensions_type x, typename std::allocator_traits<Alloc>::const_void_pointer hint) :
		array_alloc{}, ref(array_alloc::allocate(typename static_array::layout_t{x}.num_elements(), hint), x)
	{}
//	template<class Elem, typename = std::enable_if_t<std::is_convertible<Elem, typename static_array::element>{} and D==0>>
//	static_array(Elem const& e)  //2
//	:	static_array(multi::iextensions<D>{}, e){}

//	explicit static_array(typename static_array::index n, typename static_array::value_type const& v, typename static_array::allocator_type const& a = {})
//	: 	static_array(typename static_array::index_extension(n), v, a){}
	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>>
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v, typename static_array::allocator_type const& a = {}) //3
	= delete;
//	: static_array(e*extensions(v), a) {
//		adl::fill(this->begin(), this->end(), v); // TODO this should be alloc_unintialized_fill
//	}
//	template<class Allocator, typename = std::enable_if_t<std::is_same<Allocator, allocator_type>{}> >
//	explicit 
// analgous to std::vector::vector ((4)) https://en.cppreference.com/w/cpp/container/vector/vector
	explicit static_array(typename static_array::extensions_type x, typename static_array::allocator_type const& a = typename static_array::allocator_type{}) :
		array_alloc{a}, ref{array_alloc::allocate(typename static_array::layout_t{x}.num_elements()), x}
	{
		if(not std::is_trivially_default_constructible<typename static_array::element_type>{})
			uninitialized_default_construct();
	}
	template<class TT, class... Args, 
		class=std::enable_if_t<std::is_assignable<typename ref::element_ref, typename multi::basic_array<TT, D, Args...>::element>{}>,
		class=decltype(adl_copy(std::declval<multi::basic_array<TT, D, Args...> const&>().begin(), std::declval<multi::basic_array<TT, D, Args...> const&>().end(), std::declval<typename static_array::iterator>()))
	>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(multi::basic_array<TT, D, Args...> const& o, typename static_array::allocator_type const& a = {})
	: static_array(o.extensions(), a) // TODO: should be uninitialized_copy
	{
		adl_copy(o.begin(), o.end(), this->begin()); // TODO: should be uninitialized_copy, and recursive
	}

	template<class TT, class... Args>
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(array_ref<TT, D, Args...>&& o)
	: array_alloc{}, ref{array_alloc::allocate(o.num_elements()), o.extensions()}{
		static_array::uninitialized_copy_elements(std::move(o).data_elements());
	}
	static_array(static_array const& o, typename static_array::allocator_type const& a) //5b
	: array_alloc{a}, ref{static_array::allocate(num_elements(o)), extensions(o)}{
		uninitialized_copy_from(data(o));
	}
	static_array(static_array const& o)                                  //5b
	: array_alloc{o.get_allocator()}, ref{array_alloc::allocate(num_elements(o)), extensions(o)}{
		uninitialized_copy_elements(o.data_elements());
	}
//	TODO static_array(static_array&& o)                                  //5b'
//	: array_alloc{o.get_allocator()}, ref{array_alloc::allocate(num_elements(o)), extensions(o)}{
//		array_alloc::uninitialized_move_elements(data_elements(o));
//	}
	// cppcheck-suppress noExplicitConstructor ; to allow assigment-like construction of nested arrays
	constexpr static_array(std::initializer_list<typename static_array<T, D>::value_type> mil) : 
		static_array(static_array<T, D>(mil.begin(), mil.end()))
	{}
	static_array(
		std::initializer_list<typename static_array<T, D>::value_type> mil, 
		typename static_array::allocator_type const& a
	) : static_array(static_array<T, D>(mil.begin(), mil.end()), a){}
	template<class TT, std::size_t N> 
	// cppcheck-suppress noExplicitConstructor ; to allow assigment-like construction from c-arrays
	constexpr static_array(TT(&array)[N]) : static_array(std::begin(array), std::end(array)){}
	template<class It> static auto distance(It a, It b){using std::distance; return distance(a, b);}
protected:
	void deallocate(){
		if(this->num_elements()) alloc_traits::deallocate(this->alloc(), this->base_, static_cast<typename alloc_traits::size_type>(this->num_elements()));
	}
	void clear() noexcept{
		this->destroy();
		deallocate();
		layout_t<D>::operator=({});
	}
	template<class... Ts>
	constexpr static_array&& reindex(Ts... a)&&{
		static_array::layout_t::reindex(a...);
		return std::move(*this);
	}
	template<class... Ts>
	constexpr static_array& reindex(Ts... a)&{
		static_array::layout_t::reindex(a...);
		return *this;
	}
public:
	static_array() = default;
	~static_array() noexcept{destroy(); deallocate();}
	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element const>;
	using element_move_ptr  = std::move_iterator<typename static_array::element_ptr>;
	using reference = typename std::conditional<
		(static_array::dimensionality > 1), 
		basic_array<typename static_array::element, static_array::dimensionality-1, typename static_array::element_ptr>, 
		typename std::conditional<
			static_array::dimensionality == 1,
			typename std::iterator_traits<typename static_array::element_ptr>::reference,
			void
		>::type
	//	typename pointer_traits<typename static_array::element_ptr>::element_type&
	>::type;
	using const_reference = typename std::conditional<
		(static_array::dimensionality > 1), 
		basic_array<typename static_array::element, static_array::dimensionality-1, typename static_array::element_const_ptr>, // TODO should be const_reference, but doesn't work witn rangev3
		typename std::conditional<
			static_array::dimensionality == 1,
			decltype(*std::declval<typename static_array::element_const_ptr>()),
		//	typename std::iterator_traits<typename static_array::element_const_ptr>::reference,
			void
		>::type
	//	typename pointer_traits<typename static_array::element_ptr>::element_type const&
	>::type;
	using       iterator = multi::array_iterator<T, static_array::dimensionality, typename static_array::element_ptr      >;//, reference>;
	using const_iterator = multi::array_iterator<T, static_array::dimensionality, typename static_array::element_const_ptr>;//, const_reference>;
//	reference
//	      reference operator[](index i)     &&{return std::move(*this).ref::operator[](i);}
//	      reference operator[](index i)      &{return ref::operator[](i);}
//	const_reference operator[](index i) const&{return ref::operator[](i);}
//	typename static_array::allocator_type get_allocator() const{return static_cast<typename static_array::allocator_type const&>(*this);}
	friend typename static_array::allocator_type get_allocator(static_array const& self){return self.get_allocator();}

//	[[deprecated("use ::data_elements()")]] typename static_array::element_ptr data()       {return ref::data_elements();}
//	[[deprecated("use ::data_elements()")]] constexpr auto data() const{return typename static_array::element_const_ptr{ref::data_elements()};}
//#ifndef __NVCC__ // deprecated friend doesn't work in nvcc
//	[[deprecated("use data_elements()")]] 
//#else
//	__attribute__((deprecated))
//#endif
//	friend typename static_array::element_ptr       data(static_array&       s){return s.data_elements();}
//#ifndef __NVCC__ // deprecated friend doesn't work in nvcc
//	[[deprecated("use data_elements()")]] 
//#else
//	__attribute__((deprecated))
//#endif
//	friend typename static_array::element_const_ptr data(static_array const& s){return s.data_elements();}

	element_const_ptr                   data_elements() const&{return this->base_;}
	typename static_array::element_ptr  data_elements()      &{return this->base_;}
	static_array::element_move_ptr      data_elements()     &&{return std::make_move_iterator(this->base_);}

	friend auto data_elements(static_array const& s){return           s .data_elements();}
	friend auto data_elements(static_array      & s){return           s .data_elements();}
	friend auto data_elements(static_array     && s){return std::move(s).data_elements();}

	constexpr typename static_array::element_ptr       base()      {return ref::base();}
	constexpr typename static_array::element_const_ptr base() const{return typename static_array::element_const_ptr{ref::base()};}
	friend typename static_array::element_ptr       base(static_array&       s){return s.base();}
	friend typename static_array::element_const_ptr base(static_array const& s){return s.base();}

	typename static_array::element_ptr       origin()      {return ref::origin();}
	typename static_array::element_const_ptr origin() const{return ref::origin();}
	friend typename static_array::element_ptr       origin(static_array&       s){return s.origin();}
	friend typename static_array::element_const_ptr origin(static_array const& s){return s.origin();}

//	template<class... Args> decltype(auto) operator()(Args const&... args)&{return ref::operator()(args...);}
//	template<class... Args> decltype(auto) operator()(Args const&... args) const&{return ref::operator()(args...);}
//	using ref::operator();

//	basic_array<T, D, typename static_array::element_ptr> 
//	decltype(auto) operator()() &{
//		this->template static_array_cast<typename static_array::element_type>();
	//	return static_array_cast<typename static_array::element_type>(*this);//(std::forward<Ts>(t)...);
	//	return ref::operator()();
	//	return *this;
//	}
//	basic_array<T, D, typename static_array::element_const_ptr> operator()() const&{
//		this->template static_array_cast<typename static_array::element_type>();
	//	return static_array_cast<typename static_array::element_type const>(*this);
	//	return basic_array<T, D, typename static_array::element_const_ptr>{this->layout(), this->base_};
//	}
//	template<class... Ts> decltype(auto) operator()(Ts&&... t) &     {assert(0); return ref::operator()(std::forward<Ts>(t)...);}
//	template<class... Ts> decltype(auto) operator()(Ts&&... t) &&    {return std::move(*this).ref::operator()(std::forward<Ts>(t)...);}
//	template<class... Ts> decltype(auto) operator()(Ts&&... t) const&{return ref::operator()(std::forward<Ts>(t)...);}

//	template<class... Ts> decltype(auto) operator()(Ts&&... t) const{return static_array_cast<typename static_array::element_type const>(*this)(std::forward<Ts>(t)...);}

#if 0
	template<class... As> decltype(auto) paren(index a, As... as) &     {return ref::paren(a, std::forward<As>(as)...);}
	template<class... As> decltype(auto) paren(index a, As... as) && = delete;//{return ref::operator()(a, std::forward<As>(as)...);}
	template<class... As> decltype(auto) paren(index a, As... as) const&{return ref::paren(a, std::forward<As>(as)...);}

	template<class... As> decltype(auto) paren(index_range a, As... as) &{return ref::paren(a, std::forward<As>(as)...);}
	template<class... As> decltype(auto) paren(index_range a, As... as) && = delete;//{return ref::operator()(a, std::forward<As>(as)...);}
	template<class... As> decltype(auto) paren(index_range a, As... as) const&{return ref::paren(a, std::forward<As>(as)...);}

	decltype(auto) operator()(index i) & {return operator[](i);}
	decltype(auto) operator()(index i) && {return std::move(*this).operator[](i);}
	decltype(auto) operator()(index i) const& {return operator[](i);}
//#define SARRAY1(A1) auto operator()(A1 a1) const{return operator()<>(a1);}
#define SARRAY2(A1, A2)	\
	auto operator()(A1 a1, A2 a2) const& JUSTRETURN(paren<A2>(a1, a2)) \
	auto operator()(A1 a1, A2 a2) && = delete;/*     JUSTRETURN(std::move(*this).static_array::template paren<A2>(a1, a2))*/ \
	auto operator()(A1 a1, A2 a2) &      JUSTRETURN(paren<A2>(a1, a2))  
	SARRAY2(index, index ); SARRAY2(irange, index );
	SARRAY2(index, irange); SARRAY2(irange, irange);
#undef SARRAY2
#if 0
#define SARRAY3(A1, A2, A3) auto operator()(A1 a1, A2 a2, A3 a3) const{return operator()<A2, A3>(a1, a2, a3);} auto operator()(A1 a1, A2 a2, A3 a3){return operator()<A2, A3>(a1, a2, a3);}
	SARRAY3(index, index , index ); SARRAY3(irange, index , index );
	SARRAY3(index, index , irange); SARRAY3(irange, index , irange);
	SARRAY3(index, irange, index ); SARRAY3(irange, irange, index );
	SARRAY3(index, irange, irange); SARRAY3(irange, irange, irange);
#undef SARRAY3
#define SARRAY4(A1, A2, A3, A4) auto operator()(A1 a1, A2 a2, A3 a3, A4 a4) const{return operator()<A2, A3, A4>(a1, a2, a3, a4);} auto operator()(A1 a1, A2 a2, A3 a3, A4 a4) {return operator()<A2, A3, A4>(a1, a2, a3, a4);}
	SARRAY4(index, index, index , index ); SARRAY4(index, irange, index , index );
	SARRAY4(index, index, index , irange); SARRAY4(index, irange, index , irange);
	SARRAY4(index, index, irange, index ); SARRAY4(index, irange, irange, index );
	SARRAY4(index, index, irange, irange); SARRAY4(index, irange, irange, irange);
	SARRAY4(irange, index, index , index ); SARRAY4(irange, irange, index , index );
	SARRAY4(irange, index, index , irange); SARRAY4(irange, irange, index , irange);
	SARRAY4(irange, index, irange, index ); SARRAY4(irange, irange, irange, index );
	SARRAY4(irange, index, irange, irange); SARRAY4(irange, irange, irange, irange);
#undef SARRAY4
#endif
#endif
//	using const_reverse_iterator = basic_reverse_iterator<const_iterator>;
	constexpr auto rotated(dimensionality_type d = 1) const&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, D, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto rotated(dimensionality_type d = 1)&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}
	constexpr auto rotated(dimensionality_type d = 1)&&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}
//	friend decltype(auto) rotated(static_array const& self){return self.rotated();}
//	template<class Array, typename = std::enable_if_t<std::is_same<static_array, std::decay_t<Array>>{}> > 
	friend constexpr decltype(auto) rotated(static_array&       s){return s.rotated();}
	friend constexpr decltype(auto) rotated(static_array const& s){return s.rotated();}

	constexpr auto unrotated(dimensionality_type d = 1) const&{
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, D, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto unrotated(dimensionality_type d = 1)&{
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, D, typename static_array::element_ptr>{new_layout, this->base_};
	}
	friend constexpr decltype(auto) unrotated(static_array& self){return self.unrotated();}
	friend constexpr decltype(auto) unrotated(static_array const& self){return self.unrotated();}

	constexpr decltype(auto) operator<<(dimensionality_type d)      {return   rotated(d);}
	constexpr decltype(auto) operator>>(dimensionality_type d)      {return unrotated(d);}
	constexpr decltype(auto) operator<<(dimensionality_type d) const{return   rotated(d);}
	constexpr decltype(auto) operator>>(dimensionality_type d) const{return unrotated(d);}

	static_array& operator=(static_array const& other) &{
		assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}
	template<class TT, class... As>
	static_array& operator=(static_array<TT, static_array::dimensionality, As...> const& other)&{assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}
//	template<class... As>
//	static_array operator=(static_array<static_array::value_type, static_array::dimensionality, As...> const& o){assert( extensions(o) == static_array::extensions() );
//		return adl::copy_elements(o.data_elements()), *this;
//	}
	constexpr operator basic_array<typename static_array::value_type, static_array::dimensionality, typename static_array::element_const_ptr, typename static_array::layout_t>()&{
		return this->template static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
//		return static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}
};

template<class T, class Alloc>
struct static_array<T, dimensionality_type{0}, Alloc> : 
	protected array_allocator<Alloc>,
	public array_ref<T, 0, typename std::allocator_traits<typename array_allocator<Alloc>::allocator_type>::pointer>
{
private:
	using array_alloc = array_allocator<Alloc>;
public:	
	static_assert( std::is_same<typename std::allocator_traits<Alloc>::value_type, typename static_array::element>{}, 
		"allocator value type must match array value type");
	using array_alloc::get_allocator;
	using allocator_type = typename static_array::allocator_type;
protected:
	using alloc_traits = typename std::allocator_traits<allocator_type>;
	using ref = array_ref<T, 0, typename std::allocator_traits<typename std::allocator_traits<Alloc>::template rebind_alloc<T>>::pointer>;
	auto uninitialized_value_construct(){return adl_alloc_uninitialized_value_construct_n(static_array::alloc(), this->base_, this->num_elements());}
	template<typename It> auto uninitialized_copy(It first){return adl_alloc_uninitialized_copy_n(this->alloc(), first, this->num_elements(), this->data_elements());}
	template<typename It>
	auto uninitialized_move(It first){
		return adl_alloc_uninitialized_move_n(this->alloc(), first, this->num_elements(), this->data_elements());
	}
	void destroy(){array_alloc::destroy_n(this->data_elements(), this->num_elements());}
public:
	using typename ref::value_type;
	using typename ref::size_type;
	using typename ref::difference_type;
	constexpr explicit static_array(allocator_type const& a) : array_alloc{a}{}
protected:
	constexpr static_array(static_array&& other, allocator_type const& a)                           //6b
	:	array_alloc{a},
		ref{other.base_, other.extensions()}
	{
		other.ref::layout_t::operator=({});
	}
public:
	using ref::operator==;
	using ref::operator!=;
	
	template<class TT, 
		std::enable_if_t<std::is_constructible<T, TT>{}, int>* = nullptr, 
		class = decltype(adl_uninitialized_copy_n(&std::declval<TT&>(), 1, std::declval<typename static_array::element_ptr&>()))
	>
	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions
	static_array(TT&& tt) : ref(static_array::allocate(typename static_array::layout_t{}.num_elements()), {}){
		adl_uninitialized_copy_n(&tt, 1, this->base());
	}
//	template<class Range0, 
//		std::enable_if_t<std::is_constructible<typename static_array::element, typename Range0::element>{}, int>* = nullptr, 
//		class = decltype(adl_uninitialized_copy_n(&std::declval<Range0&>(), 1, std::declval<typename static_array::element_ptr&>()))
//	>
////	template<class Range0, class = decltype(adl_uninitialized_copy_n(&std::declval<Range0&>(), 1, std::declval<typename static_array::element_ptr&>()))>
//	// cppcheck-suppress noExplicitConstructor ; because I want to use equal for lazy assigments form range-expressions
//	static_array(Range0&& r) : ref(static_array::allocate(typename static_array::layout_t{}.num_elements()), {}){
//		adl_uninitialized_copy_n(&r, 1, this->base());
//	}
	static_array(typename static_array::extensions_type x, typename static_array::element const& e, allocator_type const& a) : //2
		array_alloc{a}, 
		ref(static_array::allocate(typename static_array::layout_t{x}.num_elements()), x)
	{
		uninitialized_fill(e);
	}
	static_array(typename static_array::element_type const& e, allocator_type const& a)
		: static_array(typename static_array::extensions_type{}, e, a){}
	auto uninitialized_fill(typename static_array::element const& e){array_alloc::uninitialized_fill_n(this->base_, this->num_elements(), e);}
	static_array(typename static_array::extensions_type const& x, typename static_array::element const& e)  //2
		: array_alloc{}, ref(static_array::allocate(typename static_array::layout_t{x}.num_elements()), x)
	{
		uninitialized_fill(e);
	}

	static_array() : static_array(multi::iextensions<0>{}){}

	explicit static_array(typename static_array::element const& e)  //2
		: static_array(multi::iextensions<0>{}, e)
	{}

	template<class ValueType, typename = std::enable_if_t<std::is_same<ValueType, typename static_array::value_type>{}>> 
	explicit static_array(typename static_array::index_extension const& e, ValueType const& v, allocator_type const& a = {}) //3
		: static_array(e*extensions(v), a)
	{
		using std::fill; fill(this->begin(), this->end(), v);
	}

	explicit static_array(typename static_array::extensions_type const& x, allocator_type const& a = {}) //3
	: array_alloc{a}, ref{static_array::allocate(typename static_array::layout_t{x}.num_elements()), x}{
		if(not std::is_trivially_default_constructible<typename static_array::element>{}) uninitialized_value_construct();
	}
	template<class TT, class... Args> 
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(multi::basic_array<TT, 0, Args...> const& other, allocator_type const& a = {})
		: array_alloc{a}, ref(static_array::allocate(other.num_elements()), extensions(other))
	{
		using std::copy; copy(other.begin(), other.end(), this->begin());
	}
	template<class TT, class... Args> 
	// cppcheck-suppress noExplicitConstructor ; because argument can be well-represented
	static_array(array_ref<TT, 0, Args...> const& other)
	:	array_alloc{}, ref{static_array::allocate(other.num_elements()), extensions(other)}{
		uninitialized_copy_(other.data_elements());
	}
	static_array(static_array const& other, allocator_type const& a)                      //5b
	:	array_alloc{a}, ref{static_array::allocate(other.num_elements()), extensions(other)}{
	//	assert(0);
		uninitialized_copy_(other.data_elements());
	}
	static_array(static_array const& o) :                                  //5b
		array_alloc{o.get_allocator()}, 
		ref{static_array::allocate(o.num_elements()), o.extensions()}
	{
		uninitialized_copy(o.data_elements());
	}
	static_array(static_array&& o) :       // it is private because it is a valid operation for derived classes //5b
		array_alloc{o.get_allocator()}, 
		ref{static_array::allocate(o.num_elements()), o.extensions()}
	{
		uninitialized_move(o.data_elements()); // TODO: uninitialized_move?
	}
	template<class It> static auto distance(It a, It b){using std::distance; return distance(a, b);}
protected:
	void deallocate(){ // TODO move this to array_allocator
		if(this->num_elements()) alloc_traits::deallocate(this->alloc(), this->base_, static_cast<typename alloc_traits::size_type>(this->num_elements()));
	}
	void clear() noexcept{
		this->destroy();
		deallocate();
		layout_t<0>::operator=({});
	}
public:
//	static_array() = default;
	~static_array() noexcept{
		this->destroy();
		deallocate();
	}
	using element_const_ptr = typename std::pointer_traits<typename static_array::element_ptr>::template rebind<typename static_array::element const>;
	friend allocator_type get_allocator(static_array const& self){return self.get_allocator();}

//	[[deprecated("use data_elements() instead of data()")]]
//	constexpr typename static_array::element_ptr       data()      {return ref::data_elements();}
//	[[deprecated("use data_elements() instead of data()")]]
//	constexpr auto data() const{return typename static_array::element_const_ptr{ref::data_elements()};}

	// TODO find how to use `deprecated` with nvcc
//	friend constexpr typename static_array::element_ptr       data(static_array&       s)
//	{return s.data_elements();}
//	friend constexpr typename static_array::element_const_ptr data(static_array const& s)
//	{return s.data_elements();}

	       constexpr typename static_array::element_ptr       base()                 &   {return ref::base();}
	       constexpr typename static_array::element_const_ptr base()            const&   {return ref::base();}
	friend constexpr typename static_array::element_ptr       base(static_array&       s){return s.base();}
	friend constexpr typename static_array::element_const_ptr base(static_array const& s){return s.base();}

	constexpr typename static_array::element_ptr       origin()      {return ref::origin();}
	constexpr typename static_array::element_const_ptr origin() const{return ref::origin();}
	friend constexpr typename static_array::element_ptr       origin(static_array&       s){return s.origin();}
	friend constexpr typename static_array::element_const_ptr origin(static_array const& s){return s.origin();}

//	template<class... Args> decltype(auto) operator()(Args const&... args)&{return ref::operator()(args...);}
//	template<class... Args> decltype(auto) operator()(Args const&... args) const&{return ref::operator()(args...);}

//	using ref::operator typename static_array::element_ref;//{return *(this->base_);}
//	operator decltype(auto)(){return *(this->base_);}
//	operator decltype(auto)() const{return *(this->base_);}
//	operator typename static_array::element_type() const{return *(this->base_);}
//	using ref::operator();	
//	operator typename static_array::element_ref() const& = delete;//{return *(this->base_);}
	constexpr operator typename std::iterator_traits<typename static_array::element_const_ptr>::reference() const&{
		return *(this->base_);
	}
	constexpr operator typename std::add_rvalue_reference<typename std::iterator_traits<typename static_array::element_ptr>::reference>::type()&&{
		return *(this->base_);
	}
	constexpr operator typename std::iterator_traits<typename static_array::element_ptr>::reference()&{
		return *(this->base_);
	}
	constexpr explicit operator typename std::iterator_traits<typename static_array::element_const_ptr>::value_type(){
		return *(this->base_);
	}
//	basic_array<T, D, typename static_array::element_ptr> 
//	decltype(auto) operator()()&{
//		return ref::operator()();
	//	return *this;
//	}
//	basic_array<T, 0, typename static_array::element_const_ptr> operator()() const&{
//		return basic_array<T, 0, typename static_array::element_const_ptr>{this->layout(), this->base_};
//	}

//	typename std::iterator_traits<typename static_array::element_const_ptr>::reference operator()() const&{
//		return *(this->base_);
//	}


//	using const_reverse_iterator = basic_reverse_iterator<const_iterator>;
	constexpr auto rotated(dimensionality_type d = 1) const&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto rotated(dimensionality_type d = 1)&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}
	constexpr auto rotated(dimensionality_type d = 1)&&{
		typename static_array::layout_t new_layout = *this;
		new_layout.rotate(d);
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}
//	friend decltype(auto) rotated(static_array const& self){return self.rotated();}
//	template<class Array, typename = std::enable_if_t<std::is_same<static_array, std::decay_t<Array>>{}> > 
	friend constexpr decltype(auto) rotated(static_array&       s){return s.rotated();}
	friend constexpr decltype(auto) rotated(static_array const& s){return s.rotated();}

	constexpr auto unrotated(dimensionality_type d = 1) const&{
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, 0, typename static_array::element_const_ptr>{new_layout, this->base_};
	}
	constexpr auto unrotated(dimensionality_type d = 1)&{
		typename static_array::layout_t new_layout = *this;
		new_layout.unrotate(d);
		return basic_array<T, 0, typename static_array::element_ptr>{new_layout, this->base_};
	}
	friend constexpr decltype(auto) unrotated(static_array& self){return self.unrotated();}
	friend constexpr decltype(auto) unrotated(static_array const& self){return self.unrotated();}

	constexpr decltype(auto) operator<<(dimensionality_type d){return rotated(d);}
	constexpr decltype(auto) operator>>(dimensionality_type d){return unrotated(d);}
	constexpr decltype(auto) operator<<(dimensionality_type d) const{return rotated(d);}
	constexpr decltype(auto) operator>>(dimensionality_type d) const{return unrotated(d);}

	static_array& operator=(static_array const& other){assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}
	template<class TT, class... As>
	constexpr static_array& operator=(static_array<TT, static_array::dimensionality, As...> const& other)&{assert( extensions(other) == static_array::extensions() );
		adl_copy_n(other.data_elements(), other.num_elements(), this->data_elements());
		return *this;
	}

	constexpr operator basic_array<typename static_array::value_type, static_array::dimensionality, typename static_array::element_const_ptr, typename static_array::layout_t>()&{
		return this->template static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>();
	//	return static_array_cast<typename static_array::value_type, typename static_array::element_const_ptr>(*this);
	}
	
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int version){
//		auto extensions = this->extensions();
//		using boost::serialization::make_nvp;
//		ar & make_nvp("extensions", extensions);
//		if(extensions != this->extensions()){clear(); this->reextent(extensions);}
//		assert(extensions == this->extensions());
		std::move(*this).ref::serialize(ar, version);
	}
};

template<typename T, class Alloc>
struct array<T, dimensionality_type{0}, Alloc>
	: static_array<T, 0, Alloc>
{
	using static_ = static_array<T, 0, Alloc>;
	using static_::static_;
	void reextent(typename array::extensions_type const&){}
};

template<class T, dimensionality_type D, class Alloc>
struct array : static_array<T, D, Alloc>,
	boost::multi::random_iterable<array<T, D, Alloc> >
{
	using static_ = static_array<T, D, Alloc>;
	static_assert(std::is_same<typename array::alloc_traits::value_type, T>{} or std::is_same<typename array::alloc_traits::value_type, void>{}, "!");
public:
//	array_ptr<T, D, typename array::element_const_ptr> operator&() const&{return {this->base(), this->extensions()};}
//	array_ptr<T, D, typename array::element_ptr> operator&() &{return {this->base(), this->extensions()};}
//	array_ptr<T, D, typename array::element_ptr> operator&() && = delete;
	template<class Archive>
	auto serialize(Archive& ar, const unsigned int version){
		auto extensions = this->extensions();
		ar & multi::archive_traits<Archive>::make_nvp("extensions", extensions);
		if(extensions != this->extensions()){clear(); this->reextent(extensions);}
		assert(extensions == this->extensions());
		static_::serialize(ar, version);
	}
	using static_::static_;
	using typename static_::value_type;
	array() = default;
	array(array const&) = default;
	array& reshape(typename array::extensions_type x) &{
		typename array::layout_t new_layout{x};
		assert( new_layout.num_elements() == this->num_elements() );
		static_cast<typename array::layout_t&>(*this)=new_layout;
		return *this;
	}
	array& clear() noexcept{
		static_::clear();
		return *this;
	}
	friend array& clear(array& self) noexcept{return self.clear();}

	friend auto data_elements(array const& self){return self.data_elements();}
	friend auto data_elements(array      & self){return self.data_elements();}
	friend auto data_elements(array     && self){return std::move(self).data_elements();}
	
	basic_array<typename array::element, array::dimensionality, multi::move_ptr<typename array::element> >
	move() &{
		basic_array<typename array::element, array::dimensionality, multi::move_ptr<typename array::element> >
		ret = multi::static_array_cast<typename array::element, multi::move_ptr<typename array::element>>(*this);
		layout_t<array::dimensionality>::operator=({});
		return ret;
	}
	friend 	
	basic_array<typename array::element, array::dimensionality, multi::move_ptr<typename array::element> >
	move(array& self){return self.move();}

//	explicit	
//	array(array const& other)                                              // 5a
//	:	allocator_type{other}, ref{allocate(other.num_elements()), extensions(other)}{
//		uninitialized_copy_(other.data());
//	}
//	explicit
//s	using static_::static_array;
//	template<class Array> array(Array const& arr) : static_{arr}{}
//	using static_::static_;
//	template<class... As>
//	array(typename array::extensions_type x, As&&... as) : static_{x, std::forward<As>(as)...}{} //2
//	array(array const& other) : static_{static_cast<static_ const&>(other)}{}
	array(array&& o, typename array::allocator_type const& a) noexcept : static_{std::move(o), a}{}
	array(array&& o) noexcept : array{std::move(o), o.get_allocator()}{}
	friend typename array::allocator_type get_allocator(array const& self){return self.get_allocator();}
#if 0
	template<class A//, typename = std::enable_if_t<not std::is_base_of<array, std::decay_t<A>>{}>,
//		typename = std::enable_if_t<not std::is_convertible<std::decay_t<A>, typename array::element_type>{}>
	>
	array& operator=(A&& a){
		auto ext = extensions(a);
		if(ext==array::extensions()){
		//	const_cast<array const&>
			std::move(*this).static_::ref::operator=(std::forward<A>(a));
		}else{
			this->clear(); //	this->ref::layout_t::operator=(layout_t<D>{extensions(a)}); //			this->base_ = allocate(this->num_elements());
			this->base_ = this->allocate(static_cast<typename array::alloc_traits::size_type>(this->static_::ref::layout_t::operator=(layout_t<D>{extensions(a)}).num_elements()));
			using std::begin; using std::end;
		//	alloc_uninitialized_copy(this->alloc(), begin(std::forward<A>(a)), end(std::forward<A>(a)), array::begin()); //	recursive_uninitialized_copy<D>(alloc(), begin(std::forward<A>(a)), end(std::forward<A>(a)), array::begin());
			adl::alloc_uninitialized_copy(this->alloc(), begin(std::forward<A>(a)), end(std::forward<A>(a)), array::begin()); //	recursive_uninitialized_copy<D>(alloc(), begin(std::forward<A>(a)), end(std::forward<A>(a)), array::begin());
		}
		return *this;
	}
#endif
#if 0
	template<class... As>
	array& operator=(array<array::value_type, array::dimensionality, As...> const& o){
		if(o.extensions()==array::extensions()) return static_::operator=(o);
		return operator=(array{o});
	/*
		else{
			array::destroy();
			array::deallocate();
			this->static_::ref::layout_t::operator=(layout_t<D>{extensions(o)});
			array::allocate();//array::num_elements());
			array::uninitialized_copy_elements(data_elements(o));//adl::alloc_uninitialized_copy_n(data_elements(o), num_elements(o), array::data_elements());
		}*/
	//	return *this;
	}
	template<class OtherT, class... As>
	auto operator=(array<OtherT, array::dimensionality, As...> const& o) &{
		if(extensions(o)==array::extensions()){
			static_::operator=(o);
		}else{
			array::destroy(); array::deallocate();
			this->static_::ref::layout_t::operator=(layout_t<D>{extensions(o)});
			array::allocate(); array::uninitialized_copy_elements(data_elements(o));
		}
		return *this;
	}
//	array& operator=(array other) noexcept{swap(other);} //TODO consider this
#endif
	void swap(array& other) noexcept{
		using std::swap;
		swap(this->alloc(), other.alloc());
		swap(this->base_, other.base_);
		swap(
			static_cast<typename array::layout_t&>(*this), 
			static_cast<typename array::layout_t&>(other)
		);
	}
#ifndef NOEXCEPT_ASSIGNMENT
	array& operator=(array&& other) noexcept{
		using std::exchange;
		clear();
		this->base_ = exchange(other.base_, nullptr);
		this->alloc() = std::move(other.alloc());
		static_cast<typename array::layout_t&>(*this) = exchange(static_cast<typename array::layout_t&>(other), {});
		return *this;
	}
	array& operator=(array const& o){
		if(array::extensions() == o.extensions()) static_::operator=(o);
		else operator=(array{o});
		return *this;
	}
#else
	array& operator=(array o) noexcept{return swap(o), *this;}
#endif
	template<class Range, class=std::enable_if_t<not std::is_base_of<array, std::decay_t<Range>>{}> >
	auto operator=(Range&& o) // check that LHS is not read-only
	->decltype(                                         static_::operator=(o)                     , std::declval<array&>()){
		return ((array::extensions() == o.extensions())?static_::operator=(o):operator=(array(o))), *this                 ;}

	array& operator=(basic_array<T, D, multi::move_ptr<typename array::element, typename array::element_ptr>>& other){
		if(other.layout() != this->layout()) return array::operator=(other.template static_array_cast<typename array::element, typename array::element_ptr>());
		if(this->base_ != other.base_) other.base_ = nullptr;
		return *this;
	}
	friend void swap(array& a, array& b){a.swap(b);}
	void assign(typename array::extensions_type x, typename array::element const& e){
		if(array::extensions()==x){
			fill_n(this->base_, this->num_elements(), e);
		}else{
			this->clear();
			this->layout_t<D>::operator=(layout_t<D>{x});
			this->base_ = this->allocate();
			uninitialized_fill_n(e);
		//	recursive_uninitialized_fill<dimensionality>(alloc(), begin(), end(), e);
		}
	}
//	template<class It, class Size> It assign_n(It first, Size n){
//		if(n == array::size() and multi::extensions(*first) == multi::extensions(*array::begin())){
//			return static_::ref::assign(first);
//		}
//		this->
//	}
	template<class It>
	array& assign(It first, It last){using std::next; using std::all_of;
	//	auto const s = adl::distance(first, last);
		if(adl_distance(first, last) == array::size()){// and multi::extensions(*first) == multi::extensions(*array::begin())){
			static_::ref::assign(first);
		}else{
			this->operator=(array(first, last));
		/*
			this->clear();
			this->layout_t<D>::operator=(layout_t<D>{std::tuple_cat(std::make_tuple(index_extension{array::extension().front(), array::extension().front() + distance(first, last)}), multi::extensions(*first))});
			using std::next;
			using std::all_of;
			if(first!=last) assert( all_of(next(first), last, [x=multi::extensions(*first)](auto& e){return extensions(e)==x;}) );
			this->base_ = this->allocate();
			multi::uninitialized_copy<D>(first, last, array::begin());
		*/
		}
		return *this;
	}
	void assign(std::initializer_list<typename array::value_type> il){assign(il.begin(), il.end());}
	template<class Range> auto assign(Range&& r) &
	->decltype(assign(adl_begin(r), adl_end(r))){
		return assign(adl_begin(r), adl_end(r));}
	array& operator=(std::initializer_list<typename array::value_type> il){assign(il.begin(), il.end()); return *this;}

//	template<class TT, class... Args> bool operator==(basic_array<TT, D, Args...> const& other) const{return static_::operator==(other);}
//	template<class TT, class... Args> bool operator!=(basic_array<TT, D, Args...> const& other) const{return not operator==(other);}

	void reextent(typename array::extensions_type const& x){
		if(x == this->extensions()) return;
		array tmp(x, this->get_allocator());
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
	}
	void reextent(typename array::extensions_type const& x, typename array::element const& e){
		if(x == this->extensions()) return;
		array tmp(x, e, this->get_allocator());
		auto const is = intersection(this->extensions(), x);
		tmp.apply(is) = this->apply(is);
		swap(tmp);
	}
	template<class... Ts> constexpr array&& reindex(Ts... a)&&{array::layout_t::reindex(a...); return std::move(*this);}
	template<class... Ts> constexpr array&  reindex(Ts... a)& {array::layout_t::reindex(a...); return           *this ;}
	~array() noexcept = default;
};

#if defined(__cpp_deduction_guides)
// clang cannot recognize templated-using, so don't replace IL<IL<T>> by IL2<T>, etc
//#ifndef __clang__
//template<class T, dimensionality_type D, class A=std::allocator<T>> static_array(multi::initializer_list_t<T, D>, A={})->static_array<T, D, A>;
//template<class T, dimensionality_type D, class A=std::allocator<T>> array(multi::initializer_list_t<T, D>, A={})->array<T, D, A>;
//#else
#define IL std::initializer_list
	template<class T, class A=std::allocator<T>> static_array(IL<T>                , A={})->static_array<T,1,A>; 
	template<class T, class A=std::allocator<T>> static_array(IL<IL<T>>            , A={})->static_array<T,2,A>;
	template<class T, class A=std::allocator<T>> static_array(IL<IL<IL<T>>>        , A={})->static_array<T,3,A>; 
	template<class T, class A=std::allocator<T>> static_array(IL<IL<IL<IL<T>>>>    , A={})->static_array<T,4,A>; 
	template<class T, class A=std::allocator<T>> static_array(IL<IL<IL<IL<IL<T>>>>>, A={})->static_array<T,5,A>;

//	template<class T> array(std::initializer_list<T>)->array<T, 1>; 

	template<class T> array(IL<T>                )->array<T,1>; 
	template<class T> array(IL<IL<T>>            )->array<T,2>;
	template<class T> array(IL<IL<IL<T>>>        )->array<T,3>; 
	template<class T> array(IL<IL<IL<IL<T>>>>    )->array<T,4>; 
	template<class T> array(IL<IL<IL<IL<IL<T>>>>>)->array<T,5>;

	template<class T, class A> array(IL<T>                , A)->array<T,1,A>; 
	template<class T, class A> array(IL<IL<T>>            , A)->array<T,2,A>;
	template<class T, class A> array(IL<IL<IL<T>>>        , A)->array<T,3,A>; 
	template<class T, class A> array(IL<IL<IL<IL<T>>>>    , A)->array<T,4,A>; 
	template<class T, class A> array(IL<IL<IL<IL<IL<T>>>>>, A)->array<T,5,A>;

#undef IL
//#endif

template<class T, class A=std::allocator<T>> array(T[]                  , A={})->array<T,1,A>;

//template<class Array, class E = typename multi::array_traits<Array>::element, class A=std::allocator<E>, class=std::enable_if_t<is_allocator<A>{}>> array(Array            , A={})->array<typename multi::array_traits<Array>::element, 1, A>;

template<dimensionality_type D, class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<D>, T)->array<T, D, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<0>, T)->array<T,0, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<1>, T)->array<T,1, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<2>, T)->array<T,2, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<3>, T)->array<T,3, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<4>, T)->array<T,4, std::allocator<T>>;
	template<class T, class=std::enable_if_t<not is_allocator<T>{}> > array(iextensions<5>, T)->array<T,5, std::allocator<T>>;

template<dimensionality_type D, class T, class A, typename = std::enable_if_t<std::is_same<typename std::allocator_traits<A>::value_type, T>{}> > static_array(iextensions<D>, T, A)->static_array<T, D, A>;
	template<class T, class A, typename = std::enable_if_t<std::is_same<typename std::allocator_traits<A>::value_type, T>{}>> static_array(T, A)->static_array<T, 0, A>;
//	template<class T, class A, typename = std::enable_if_t<not is_allocator<T>{}> > static_array(iextensions<1>, T, A = {})->static_array<T, 1, A>;

template<dimensionality_type D, class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<D>, A)->array<T, D, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<0>, A)->array<T, 0, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<1>, A)->array<T, 1, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<2>, A)->array<T, 2, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<3>, A)->array<T, 3, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<4>, A)->array<T, 4, A>;
	template<class A, class=std::enable_if_t<is_allocator<A>{}>, typename T = typename std::allocator_traits<A>::value_type> array(iextensions<5>, A)->array<T, 5, A>;

	template<class T> array(iextensions<0>, T)->array<T, 0>;
	template<class T> array(iextensions<1>, T)->array<T, 1>; template<class T> array(multi::size_type, T)->array<T, 1>;
	template<class T> array(iextensions<2>, T)->array<T, 2>;
	template<class T> array(iextensions<3>, T)->array<T, 3>;

template<class T, class MR, class A=memory::allocator<T, MR>> array(iextensions<1>, T, MR*)->array<T, 1, A>;
template<class T, class MR, class A=memory::allocator<T, MR>> array(iextensions<2>, T, MR*)->array<T, 2, A>;
template<class T, class MR, class A=memory::allocator<T, MR>> array(iextensions<3>, T, MR*)->array<T, 3, A>;
template<class T, class MR, class A=memory::allocator<T, MR>> array(iextensions<4>, T, MR*)->array<T, 4, A>;
template<class T, class MR, class A=memory::allocator<T, MR>> array(iextensions<5>, T, MR*)->array<T, 5, A>;

template<class MatrixRef, class DT = typename MatrixRef::decay_type, class T = typename DT::element, dimensionality_type D = DT::dimensionality, class Alloc = typename DT::allocator_type>
array(MatrixRef)->array<T, D, Alloc>;

template<typename T, dimensionality_type D, typename P> array(basic_array<T, D, P>)->array<T, D>;
#endif

template <class T, std::size_t N>
multi::array<typename std::remove_all_extents<T[N]>::type, std::rank<T[N]>{}> 
decay(const T(&t)[N]) noexcept{
	return multi::array_cref<typename std::remove_all_extents<T[N]>::type, std::rank<T[N]>{}>(data_elements(t), extensions(t));
}

template<class T, size_t N>
struct array_traits<T[N], void, void>{
	using reference = T&;
	using element = std::remove_all_extents_t<T[N]>;
	using decay_type = multi::array<T, 1>;
};

}}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#include<cassert>
#include<numeric> // iota
#include<iostream>
#include<algorithm>
#include<vector>

#include <random>
#include <boost/timer/timer.hpp>
#include<boost/multi_array.hpp>

using std::cout;
namespace multi = boost::multi;

#if 1
template<class Matrix, class Vector>
void solve(Matrix& m, Vector& y){
//	using std::size; // assert(size(m) == std::ptrdiff_t(size(y)));
	std::ptrdiff_t msize = size(m); 
	for(auto r = 0; r != msize; ++r){ //	auto mr = m[r]; //  auto const mrr = mr[r];// assert( mrr != 0 ); // m[r][r] = 1;
		auto mr = m[r];
		auto mrr = mr[r];
		for(auto c = r + 1; c != msize; ++c) mr[c] /= mrr;
		auto yr = (y[r] /= mrr);
		for(auto r2 = r + 1; r2 != msize; ++r2){ //	auto mr2 = m[r2]; //	auto const mr2r = mr2[r]; // m[r2][r] = 0;
			auto mr2 = m[r2];
			auto const& mr2r = mr2[r];
			auto const& mr_cr = m[r];
			for(auto c = r + 1; c != msize; ++c) mr2[c] -= mr2r*mr_cr[c];
			y[r2] -= mr2r*yr;
		}
	}
	for(auto r = msize - 1; r > 0; --r){ //	auto mtr = m.rotated(1)[r];
		auto const& yr = y[r];
		for(auto r2 = r-1; r2 >=0; --r2)
			y[r2] -= yr*m[r2][r];
	}
}
#endif

void f(boost::multi::array<double, 4> const& A){
	A[1][2];
	auto&& a = A[1][2]; (void)a; // careful, a is a reference here, don't use auto, 
	auto const& b = A[1][2]; (void)b; // use auto const& if possible
//	A[1][2][3][4] = 5; // fail, element is read-only
}

template<class C>
void set_99(C&& c){
	for(auto j : c.extension(0))
		for(auto k : c.extension(1))
				c[j][k] = 99.;
}

namespace multi = boost::multi;

template<class T> void fun(T const& t){
	std::cout << typeid(t).name() << std::endl;
}

template<class T> struct extension{};
template<class T> void gun(extension<T>){
	std::cout << typeid(T).name() << std::endl;
}

typedef double a1010[10][10];

struct A{
	double const* p;
	A(std::initializer_list<double> il){ p = &*(il.begin() + 1); };
};

template<class T> void what(T&&) = delete;

int main(){

{
	multi::array<double, 1> A1 = {1.,2.,3.}; 
	assert( A1.dimensionality==1 and A1.num_elements()==3 );

	multi::array<double, 2> A2 {
		 {1.,2.,3.},
		 {4.,5.,6.}
	};
	*A2.begin()->begin() = 99;
	assert(A2[0][0] == 99 );
}
{
	double A[2][3] = {{1.,2.,3.}, {4.,5.,6.}};
	using multi::decay;
	auto A_copy = decay(A);
}
	{
		multi::array<double, 2> A = {
			{1},
			{2},
			{3}
		};
		assert( size(A) == 3 );
		assert( size(rotated(A)) == 1 );
		assert( stride(A) == 1 );
		assert( stride(rotated(A)) == 1 );
		assert( A.extensions() );
	}
	{
		multi::array<double, 1> A = {1};
		assert( size(A) == 1 );
		assert( stride(A) == 1 );
	}
	{
		multi::array<double, 0> A = 3.;
	//	assert( stride(A) == 1 );
	}
	{
		double D3[3] = {0, 1, 2};
		multi::array<double, 1> A(D3);
		assert( A[1] == 1 );
	}
}
#endif
#endif

