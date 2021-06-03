#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0.$X -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_DETAIL_ADL_HPP
#define MULTI_DETAIL_ADL_HPP

#include<cstddef>     // std::size_t
#include<type_traits> // std::conditional_t
#include<utility>

#include<memory>    // uninitialized_copy, etc
#include<algorithm> // std::copy, std::copy_n, std::equal, etc
#include<iterator>  // begin, end

#include "../detail/memory.hpp"
#include "../config/MAYBE_UNUSED.hpp"

#if defined(__NVCC__)
#include<thrust/copy.h>
#include<thrust/equal.h>
#endif

namespace boost{namespace multi{
	template<std::size_t I> struct priority : std::conditional_t<I==0, std::true_type, struct priority<I-1>>{}; 
}}

#define DECLRETURN(ExpR) ->decltype(ExpR){return ExpR;}
#define JUSTRETURN(ExpR)                 {return ExpR;}

#define BOOST_MULTI_DEFINE_ADL(FuN) \
namespace boost{namespace multi{ \
namespace adl{ \
	namespace custom{template<class...> struct FuN##_t;} 	__attribute__((unused))  \
	static constexpr class FuN##_t{ \
		template<class... As> [[deprecated]] auto _(priority<0>,        As&&... as) const = delete; \
		template<class... As>          auto _(priority<1>,        As&&... as) const DECLRETURN(std::FuN(std::forward<As>(as)...)) \
		template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(     FuN(std::forward<As>(as)...)) \
		template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).FuN(std::forward<As>(as)...))     \
		template<class... As>          auto _(priority<4>,        As&&... as) const DECLRETURN(custom::FuN##_t<As&&...>::_(std::forward<As>(as)...)) \
	public: \
		template<class... As> auto operator()(As&&... as) const->decltype(_(priority<4>{}, std::forward<As>(as)...)){return _(priority<4>{}, std::forward<As>(as)...);} \
	} FuN; \
} \
}}

namespace boost{namespace multi{

constexpr class adl_copy_n_t{
	template<class... As>          constexpr auto _(priority<0>,        As&&... as) const DECLRETURN(              std::copy_n              (std::forward<As>(as)...)) // it is important to terminate with SFINAE
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(           thrust::copy_n              (std::forward<As>(as)...))
#endif
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   copy_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  copy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).copy_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_copy_n;

constexpr class adl_fill_n_t{
	template<         class... As> constexpr auto _(priority<0>,        As&&... as) const DECLRETURN(              std::fill_n              (std::forward<As>(as)...))
#if defined(__NVCC__)
	template<         class... As> constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(           thrust::fill_n              (std::forward<As>(as)...))
#endif
	template<         class... As> constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   fill_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  fill_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).fill_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_fill_n;


constexpr class adl_equal_fn__{
	template<         class...As> constexpr auto _(priority<1>,      As...as) const DECLRETURN(   std::equal(as...))
#ifdef THRUST_VERSION
//	template<         class...As> constexpr auto _(priority<2>,      As...as) const DECLRETURN(thrust::equal(as...))
#endif
#ifndef THRUST_VERSION
	template<         class...As> constexpr auto _(priority<3>,      As...as) const DECLRETURN(        equal(as...))
#endif
	template<class T, class...As> constexpr auto _(priority<4>, T t, As...as) const DECLRETURN(      t.equal(as...))
public:
	template<class...As> constexpr auto operator()(As...as) const DECLRETURN(_(priority<4>{}, as...))
#ifdef THRUST_VERSION
	template<class It, class...As, class=std::enable_if_t<(It::dimensionality > 1)> > 
	                     constexpr auto operator()(It begin, As... as) const DECLRETURN(_(priority<1>{}, begin, as...))
#endif
} adl_equal;

constexpr class adl_copy_fn__{
	template<class InputIt, class OutputIt, 
		class=std::enable_if_t<std::is_assignable<typename std::iterator_traits<OutputIt>::reference, typename std::iterator_traits<InputIt>::reference>{}>
	>
	                               constexpr auto _(priority<1>, InputIt first, InputIt last, OutputIt d_first)
	                                                                                const DECLRETURN(              std::copy(first, last, d_first))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(           thrust::copy(std::forward<As>(as)...))
#endif
	template<         class... As> constexpr auto _(priority<3>,        As&&... as) const DECLRETURN(                   copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<5>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).copy(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN( _(priority<5>{}, std::forward<As>(as)...) ) \
} adl_copy;

namespace adl{ \
	namespace custom{template<class...> struct fill_t;} __attribute__((unused))
	static constexpr class fill_t{ \
		template<class... As>          auto _(priority<1>,        As&&... as) const DECLRETURN(              std::fill              (std::forward<As>(as)...)) \
		template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(                   fill              (std::forward<As>(as)...)) \
		template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).fill              (std::forward<As>(as)...)) \
		template<class... As>          auto _(priority<4>,        As&&... as) const DECLRETURN(custom::           fill_t<As&&...>::_(std::forward<As>(as)...)) \
	public: \
		template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...)) \
	} fill; \
} \
}}

namespace boost{namespace multi{

template<class Alloc> 
struct alloc_construct_elem_t{
	Alloc* palloc_;
	template<class T> auto operator()(T&& p) const
	->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(p))){
		return std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(p));}
};

namespace xtd{

//template<class T>
//constexpr auto adl_to_address(const T& p) noexcept;

template<class T> // this one goes last!!!
constexpr auto to_address(const T& p) noexcept;

template<class T> 
constexpr auto _to_address(priority<0>, const T& p) noexcept
->decltype(to_address(p.operator->())){
	return to_address(p.operator->());}

template<class T>
constexpr auto _to_address(priority<1>, const T& p) noexcept
->decltype(std::pointer_traits<T>::to_address(p)){
	return std::pointer_traits<T>::to_address(p);}

template<class T, std::enable_if_t<std::is_pointer<T>{}, int> = 0>
constexpr T _to_address(priority<2>, T const& p) noexcept{
    static_assert(!std::is_function<T>{}, "!");
    return p;
}

template<class T> // this one goes last!!!
constexpr auto to_address(const T& p) noexcept
->decltype(_to_address(priority<2>{}, p))
{
	return _to_address(priority<2>{}, p);}

template<class Alloc, class ForwardIt, class Size, typename Value = typename std::iterator_traits<ForwardIt>::value_type>
ForwardIt alloc_uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size n)
//->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), Value()), first)>
{
	ForwardIt current = first;
	try{
		for (; n > 0 ; ++current, --n) std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), Value()); // !!!!!!!!!!!!!! if you are using std::complex type consider making complex default constructible (e.g. by type traits)
		//	::new (static_cast<void*>(std::addressof(*current))) Value();
		return current;
	}catch(...){
		for(; current != first; ++first) std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*first));
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T = typename std::iterator_traits<ForwardIt>::value_type>
auto alloc_uninitialized_default_construct_n(Alloc& alloc, ForwardIt first, Size n)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first)), first)>
{
	ForwardIt curr = first;
	if(std::is_trivially_default_constructible<T>{}) std::advance(curr, n);
	else{
		using _ = std::allocator_traits<Alloc>;
		try{
			for(;n > 0; ++curr, --n) _::construct(alloc, std::addressof(*curr));
		}catch(...){
			for(;curr!=first; ++first) _::destroy(alloc, std::addressof(*first));
			throw;
		}
	}
	return curr;
}

template<class ForwardIt, class Size>
ForwardIt uninitialized_default_construct_n( ForwardIt first, Size n ){
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    ForwardIt current = first;
    try {
        for (; n > 0 ; (void) ++current, --n) {
            ::new (static_cast<void*>(std::addressof(*current))) T;
        }
        return current;
    }  catch (...) {assert(0);
//        std::destroy(first, current);
        throw;
    }
}

}

template<class Alloc> struct alloc_destroy_elem_t{
	Alloc* palloc_;
	template<class T> constexpr auto operator()(T&& p) const
//	->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::forward<T>(t)...)){
	{	return std::allocator_traits<Alloc>::destroy(*palloc_, std::addressof(p));}
};

template<class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto destroy_n(BidirIt first, Size n)
->std::decay_t<decltype(std::addressof(*(first-1)), first)>{
//	first += (n-1); // nullptr case gives UB here
//	for (; n > 0; --first, --n)
//		std::allocator_traits<Alloc>::destroy(a, std::addressof(*first));
	first += n;
	for (; n != 0; --first, --n)
		std::addressof(*(first-1))->~T();
	return first;
}

template<class Alloc, class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto alloc_destroy_n(Alloc& a, BidirIt first, Size n)
->std::decay_t<decltype(std::addressof(*(first-1)), first)>{
//	first += (n-1); // nullptr case gives UB here
//	for (; n > 0; --first, --n)
//		std::allocator_traits<Alloc>::destroy(a, std::addressof(*first));
	first += n;
	for (; n != 0; --first, --n)
		std::allocator_traits<Alloc>::destroy(a, std::addressof(*(first-1)));
	return first;
}

constexpr class adl_uninitialized_copy_fn__{
	template<class InIt, class FwdIt, class=decltype(std::addressof(*FwdIt{}))> // sfinae friendy std::uninitialized_copy
	                               constexpr auto _(priority<1>, InIt f, InIt l, FwdIt d) const DECLRETURN(              std::uninitialized_copy(f, l, d))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as)       const DECLRETURN(                   uninitialized_copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as)       const DECLRETURN(  std::decay_t<T>::uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as)       const DECLRETURN(std::forward<T>(t).uninitialized_copy(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_uninitialized_copy;

namespace xtd{

template<class InputIt, class Size, class ForwardIt, class Value = typename std::iterator_traits<ForwardIt>::value_type>
auto uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first)
->std::decay_t<decltype(::new (static_cast<void*>(std::addressof(*d_first))) Value(*first), d_first)>{
	ForwardIt current = d_first;
	try{
		for (; count > 0; ++first, (void) ++current, --count) {
			::new (static_cast<void*>(std::addressof(*current))) Value(*first);
		}
	}catch(...){
		for(; d_first != current; ++d_first) d_first->~Value();
		throw;
	}
	return current;
}

template<class InputIt, class Size, class ForwardIt, class Value = typename std::iterator_traits<ForwardIt>::value_type>
auto uninitialized_move_n(InputIt first, Size count, ForwardIt d_first)
->std::decay_t<decltype(::new (static_cast<void*>(std::addressof(*d_first))) Value(std::move(*first)), d_first)>{
	ForwardIt current = d_first;
	try{
		for (; count > 0; ++first, (void) ++current, --count) {
			::new (static_cast<void*>(std::addressof(*current))) Value(std::move(*first));
		}
	}catch(...){
		for(; d_first != current; ++d_first) d_first->~Value();
		throw;
	}
	return current;
}

}

constexpr class adl_uninitialized_copy_n_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::uninitialized_copy_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   uninitialized_copy_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  uninitialized_copy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_copy_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return _(priority<4>{}, std::forward<As>(as)...);}
} adl_uninitialized_copy_n;

MAYBE_UNUSED static constexpr class adl_uninitialized_move_n_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              xtd::uninitialized_move_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   uninitialized_move_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  uninitialized_move_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_move_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return _(priority<4>{}, std::forward<As>(as)...);}
} adl_uninitialized_move_n;

namespace xtd{

template<class T, class InputIt, class Size, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
constexpr auto alloc_uninitialized_copy_n(std::allocator<T>&, InputIt f, Size n, ForwardIt d){
	return adl_uninitialized_copy_n(f, n, d);}

template<class T, class InputIt, class Size, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
constexpr auto alloc_uninitialized_move_n(std::allocator<T>&, InputIt f, Size n, ForwardIt d){
	return adl_uninitialized_move_n(f, n, d);}

template<class Alloc, class InputIt, class Size, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
auto alloc_uninitialized_copy_n(Alloc& a, InputIt f, Size n, ForwardIt d)
//->std::decay_t<decltype(a.construct(std::addressof(*d), *f), d)>
{
	ForwardIt c = d;
	try{
		for(; n > 0; ++f, ++c, --n) std::allocator_traits<Alloc>::construct(a, std::addressof(*c), *f);
		return c;
	}catch(...){
		for(; d != c; ++d) std::allocator_traits<Alloc>::destroy(a, std::addressof(*d));
		throw;
	}
}

template<class Alloc, class InputIt, class Size, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
auto alloc_uninitialized_move_n(Alloc& a, InputIt f, Size n, ForwardIt d)
//->std::decay_t<decltype(a.construct(std::addressof(*d), *f), d)>
{
	ForwardIt c = d;
	try{
		for(; n > 0; ++f, ++c, --n) std::allocator_traits<Alloc>::construct(a, std::addressof(*c), std::move(*f));
		return c;
	}catch(...){
		for(; d != c; ++d) std::allocator_traits<Alloc>::destroy(a, std::addressof(*d));
		throw;
	}
}

template<class T, class InputIt, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
constexpr auto alloc_uninitialized_copy(std::allocator<T>&, InputIt f, InputIt l, ForwardIt d){
	return adl_uninitialized_copy(f, l, d);
}

template<class Alloc, class InputIt, class ForwardIt, class=decltype(std::addressof(*std::declval<ForwardIt>())), class=std::enable_if_t<std::is_constructible<typename std::iterator_traits<ForwardIt>::value_type, typename std::iterator_traits<InputIt>::reference>{}> >
auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt d_first)
//->std::decay_t<decltype(a.construct(std::addressof(*d_first), *first), d_first)> // problematic in clang-11 + gcc-9
{
	ForwardIt current = d_first;
	try{
		for(; first != last; ++first, (void)++current){
			std::allocator_traits<std::decay_t<Alloc>>::construct(a, std::addressof(*current), *first);
		//	a.construct(std::addressof(*current), *first);
		}
		return current;
	}catch(...){
		for(; d_first != current; ++d_first){
			std::allocator_traits<std::decay_t<Alloc>>::destroy(a, std::addressof(*d_first));
		//	a.destroy(std::addressof(*d_first));
		}
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T>
auto alloc_uninitialized_fill_n(Alloc& a, ForwardIt first, Size n, T const& v)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(a, std::addressof(*first), v), first)>
{
	ForwardIt current = first; // using std::to_address;
	try{
		for(; n > 0; ++current, --n) std::allocator_traits<Alloc>::construct(a, std::addressof(*current), v);
		return current;
	}catch(...){
		for(; first != current; ++first) std::allocator_traits<Alloc>::destroy(a, std::addressof(*first)); 
		throw;
	}
}
}

}}

namespace boost{
namespace multi{

constexpr class adl_distance_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::distance(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   distance(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::distance(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).distance(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_distance;

constexpr class adl_begin_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::begin(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   begin(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::begin(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).begin(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_begin;

constexpr class adl_end_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::end(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   end(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::end(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).end(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_end;

constexpr class adl_swap_ranges_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::swap_ranges(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   swap_ranges(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::swap_ranges(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).swap_ranges(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_swap_ranges;

constexpr class adl_lexicographical_compare_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              std::lexicographical_compare(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   lexicographical_compare(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::lexicographical_compare(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).lexicographical_compare(std::forward<As>(as)...))
public:
	template<class... As> constexpr	auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_lexicographical_compare;

constexpr class adl_alloc_uninitialized_value_construct_n_fn__{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_value_construct_n(std::forward<As>(as)...)) // TODO: use boost?
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   alloc_uninitialized_value_construct_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_value_construct_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_value_construct_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return (_(priority<4>{}, std::forward<As>(as)...));}
} adl_alloc_uninitialized_value_construct_n;

constexpr class adl_alloc_uninitialized_default_construct_n_t{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...)) // TODO: use boost?
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_default_construct_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return (_(priority<4>{}, std::forward<As>(as)...));}
} adl_alloc_uninitialized_default_construct_n;

constexpr class adl_uninitialized_default_construct_n_t{
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const{return            xtd::uninitialized_default_construct_n              (std::forward<As>(as)...);}
//		template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(                   uninitialized_default_construct_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_default_construct_n              (std::forward<As>(as)...))
public:
		template<class... As> constexpr auto operator()(As&&... as) const{return (_(priority<4>{}, std::forward<As>(as)...));}
} adl_uninitialized_default_construct_n;

constexpr class destroy_n_fn__ {
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(            multi::destroy_n              (std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   destroy_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::destroy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).destroy_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_destroy_n;

constexpr class alloc_destroy_n_fn__ {
	template<class T, class... As> constexpr auto _(priority<1>, T&&  , As&&... as) const DECLRETURN(                     adl_destroy_n              (std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(multi::            alloc_destroy_n              (std::forward<As>(as)...)) // TODO: use boost?
	template<class... As>          constexpr auto _(priority<3>,        As&&... as) const DECLRETURN(                   alloc_destroy_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::alloc_destroy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<5>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_destroy_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_alloc_destroy_n;

constexpr class adl_alloc_uninitialized_copy_fn__ {
	template<class A, class... As> constexpr auto _(priority<1>, A&&  , As&&... as) const DECLRETURN(                     adl_uninitialized_copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<2>, T&& t, As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(                   alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<5>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_copy(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_alloc_uninitialized_copy;

constexpr class alloc_uninitialized_copy_n_fn__ {
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const{return(                  xtd::alloc_uninitialized_copy_n(std::forward<As>(as)...));}
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   alloc_uninitialized_copy_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_copy_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return _(priority<3>{}, std::forward<As>(as)...);} \
} adl_alloc_uninitialized_copy_n;

constexpr class alloc_uninitialized_move_n_fn__ {
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const{return(                  xtd::alloc_uninitialized_move_n(std::forward<As>(as)...));}
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                   alloc_uninitialized_move_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_move_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const{return _(priority<3>{}, std::forward<As>(as)...);} \
} adl_alloc_uninitialized_move_n;

constexpr class uninitialized_fill_n_fn__ {
	template<class... As>          constexpr auto _(priority<1>,        As&&... as) const DECLRETURN(               std::uninitialized_fill_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(                    uninitialized_fill_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>, T&& t, As&&... as) const DECLRETURN( std::forward<T>(t).uninitialized_fill_n(std::forward<As>(as)...))
public:
	template<class T1, class... As> constexpr auto operator()(T1&& t1, As&&... as) const DECLRETURN(_(priority<3>{}, t1, std::forward<As>(as)...))
} adl_uninitialized_fill_n;

constexpr class alloc_uninitialized_fill_n_fn__ {
	template<class T, class... As> constexpr auto _(priority<1>, T&&  , As&&... as) const DECLRETURN(                      adl_uninitialized_fill_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>,        As&&... as) const DECLRETURN(               xtd::alloc_uninitialized_fill_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<3>,        As&&... as) const DECLRETURN(                    alloc_uninitialized_fill_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>, T&& t, As&&... as) const DECLRETURN( std::forward<T>(t).alloc_uninitialized_fill_n(std::forward<As>(as)...))
public:
	template<class T1, class... As> constexpr auto operator()(T1&& t1, As&&... as) const DECLRETURN(_(priority<4>{}, t1, std::forward<As>(as)...))
} adl_alloc_uninitialized_fill_n;

template<dimensionality_type N>
struct recursive{
	template<class Alloc, class InputIt, class ForwardIt>
	static constexpr auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		using std::begin; using std::end;
		while(first!=last){
			recursive<N-1>::alloc_uninitialized_copy(a, begin(*first), end(*first), begin(*dest));
			++first;
			++dest;
		}
		return dest;
	}
};

template<> struct recursive<1>{
	template<class Alloc, class InputIt, class ForwardIt>
	static auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		return adl_alloc_uninitialized_copy(a, first, last, dest);
	}
};

}}


//BOOST_MULTI_DEFINE_ADL(lexicographical_compare);
//BOOST_MULTI_DEFINE_ADL(swap_ranges);

//BOOST_MULTI_DEFINE_ADL(begin);
//BOOST_MULTI_DEFINE_ADL(end);

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi ADL"
#ifdef BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#else
#include<boost/test/included/unit_test.hpp>
#endif

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_detail_adl){

	std::vector<double> v{1., 2., 3.};
	std::vector<double> w(3);

	multi::adl_copy_n(v.data(), 3, w.data());

	BOOST_REQUIRE(
		multi::adl_equal(v.data(), v.data()+v.size(), w.data())
	);

}
#endif
#endif

