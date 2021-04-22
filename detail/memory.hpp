#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x&&$0x&&rm -f $0x; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "layout.hpp"
#include "../utility.hpp"
#include<memory>
#include<algorithm> // copy_n
#include "../detail/adl.hpp"

//#include<boost/core/alloc_construct.hpp>

namespace boost{
namespace multi{

namespace memory{

template<class Alloc>
struct allocator_traits : std::allocator_traits<Alloc>{
	template<class Ptr, class... Args>
	static auto construct(Alloc& a, Ptr p, Args&&... args)
	->decltype(a.construct(p, std::forward<Args>(args)...)){
		return a.construct(p, std::forward<Args>(args)...);}

	template<class Ptr> 
	static auto destroy(Alloc& a, Ptr p)
	->decltype(a.destroy(p)){
		return a.destroy(p);}
};

}

using memory::allocator_traits;
//using memory::to_address;
////////////////////////////////////////////////////////////////////////////////
//template<class Ptr> auto to_address(Ptr const& p) noexcept
//->decltype(p.operator->()){
//	return p.operator->();}

//template<class T> constexpr T* to_address(T* p) noexcept{return p;}

//https://en.cppreference.com/w/cpp/memory/destroy
template<class Alloc, class ForwardIt, typename = std::enable_if_t<!(has_rank<ForwardIt>{})>>//, typename = std::enable_if_t<typename ForwardIt::rank{} == 1> >//, typename AT = typename std::allocator_traits<Alloc> >
void destroy(Alloc& a, ForwardIt f, ForwardIt l){
	for(; f != l; ++f) allocator_traits<Alloc>::destroy(a, std::addressof(*f));
}

template<class Alloc, class ForwardIt, typename = std::enable_if_t<has_rank<ForwardIt>{} and typename ForwardIt::rank{} == 1>>//, typename = std::enable_if_t<typename ForwardIt::rank{} == 1> >//, typename AT = typename std::allocator_traits<Alloc> >
void destroy(Alloc& a, ForwardIt first, ForwardIt last, double* = 0){
	//	using multi::to_address;
	for(; first != last; ++first) a.destroy(to_address(first)); //	AT::destroy(a, to_address(first)); //	AT::destroy(a, addressof(*first)); // a.destroy(addressof(*first));
}
template<class Alloc, class ForwardIt, typename = std::enable_if_t<has_rank<ForwardIt>{} and typename ForwardIt::rank{} != 1>>//, typename AT = typename std::allocator_traits<Alloc> >
void destroy(Alloc& a, ForwardIt first, ForwardIt last, void* = 0){
	for(; first != last; ++first) destroy(a, begin(*first), end(*first));
}

template<class Alloc, class InputIt, class Size, class ForwardIt>//, typename AT = std::allocator_traits<Alloc> >
ForwardIt uninitialized_move_n(Alloc& a, InputIt f, Size n, ForwardIt d){
	ForwardIt c = d;
//	using std::addressof;
	try{
		for(; n > 0; ++f, ++c, --n)
		//	alloc_construct(a, to_address(c), std::move(*f));
			a.construct(std::addressof(*c), std::move(*f));
		//	AT::construct(a, to_address(c), *f);
		//	AT::construct(a, addressof(*c), *f);
		//	a.construct(addressof(*c), *f);
		return c;
	}catch(...){destroy(a, d, c); throw;}
}

template<class Alloc, class ForwardIt, class Size>//, class AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_default_construct_n(Alloc& a, ForwardIt first, Size n){
	ForwardIt current = first;
	try{
		for(; n > 0; ++current, --n)
		//	allocator_traits<Alloc>::construct(a, to_address(current));
			a.construct(to_address(current));
		return current;
	}catch(...){destroy(a, first, current); throw;}
}

template<class Alloc, class ForwardIt, class Size, 
	typename T = typename std::iterator_traits<ForwardIt>::value_type,
	typename = std::enable_if_t<not std::is_trivially_default_constructible<T>{}> 
>
ForwardIt uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n){
	ForwardIt current = first; // using std::addressof;
	try{
		for(; n > 0; ++current, --n)
			allocator_traits<Alloc>::construct(a, to_address(current), T());
		//	a.construct(to_address(current), T()); //	a.construct(std::pointer_traits<Ptr>::pointer_to(*current), T()); //	AT::construct(a, to_address(current), T()); //	AT::construct(a, addressof(*current), T()); //	a.construct(addressof(*current), T());
		return current;
	}catch(...){destroy(a, first, current); throw;}
}

template<class... Args> auto std_copy(Args&&... args){
	using std::copy;
	return copy(std::forward<Args>(args)...);
}
//using std::copy_n;
//template<class Alloc, class In, typename Size, class Fwd, class=std::enable_if_t<std::is_trivially_copyable<typename Fwd::element>{}>>
//auto alloc_uninitialized_copy_n(Alloc&, In first, Size count, Fwd d_first)
//->decltype(std_copy_n(first, count, d_first)){
//	assert(0);
//	return copy_n(first, count, d_first);}

//using std::copy;
//template<class Alloc, class In, class MIt, class=std::enable_if_t<std::is_trivially_copyable<typename MIt::element>{}>>
//auto alloc_uninitialized_copy(Alloc&, In f, In l, MIt d)
//->decltype(adl::copy(f, l, d)){
//	return adl::copy(f, l, d);}

/*
template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<has_rank<MIt>{}>, typename = std::enable_if_t<typename MIt::rank{}==1>, class=std::enable_if_t<!std::is_trivially_copyable<typename MIt::element>{}> >
MIt alloc_uninitialized_copy(Alloc& a, InputIt f, InputIt l, MIt const& d){
	MIt current = d; // using multi::to_address;
	try{
		for(; f != l; ++f, ++current) a.construct(to_address(current), *f);
		return current;
	}catch(...){destroy(a, d, current); throw;}
}*/

namespace xtd{
template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<!has_rank<MIt>{}> >
MIt alloc_uninitialized_copy(Alloc& a, InputIt f, InputIt l, MIt d, double* = 0){
	MIt current = d;
//	using multi::to_address;
	try{
		for(; f != l; ++f, ++current) a.construct(to_address(current), *f);
		return current;
	}catch(...){destroy(a, d, current); throw;}
}
}

// https://en.cppreference.com/w/cpp/memory/destroy_at
template<class Alloc, class T, typename AT = std::allocator_traits<Alloc> >
void destroy_at(Alloc& a, T* p){AT::destroy(a, p);}

//template<class Ptr> auto to_address(const Ptr& p) noexcept 
//->decltype(p.operator->());

//template<class Ptr> auto to_address(Ptr const& p) noexcept 
//->decltype(p.operator->());

//template<class Alloc, class ForwardIt>
//void destroy(Alloc& a, ForwardIt first, ForwardIt last);
/*template<class Alloc, class ForwardIt, class Size>
ForwardIt alloc_destroy_n(Alloc& a, ForwardIt first, Size n);
template<class Alloc, class InputIt, class Size, class ForwardIt>
ForwardIt alloc_uninitialized_copy_n(Alloc& a, InputIt f, Size n, ForwardIt d);
template<class Alloc, class InputIt, class Size, class ForwardIt>
ForwardIt alloc_uninitialized_move_n(Alloc& a, InputIt f, Size n, ForwardIt d);
template<class Alloc, class ForwardIt, class Size, class T>
ForwardIt alloc_uninitialized_fill_n(Alloc& a, ForwardIt first, Size n, const T& v);
template<class Alloc, class ForwardIt, class Size>
ForwardIt alloc_uninitialized_default_construct_n(Alloc& a, ForwardIt first, Size n);
*/
//template<class Alloc, class ForwardIt, class Size>
//ForwardIt uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n);

// https://en.cppreference.com/w/cpp/memory/destroy_n
template<class Alloc, class ForwardIt, class Size>//, typename AT = typename std::allocator_traits<Alloc> >
ForwardIt destroy_n(Alloc& a, ForwardIt first, Size n){
//	using std::addressof;
	for(; n > 0; ++first, --n)
		allocator_traits<Alloc>::destroy(a, to_address(first));
	//	a.destroy(to_address(first));
	//	AT::destroy(a, addressof(*first));
	//	a.destroy(addressof(*first));
	return first;
}

//template<class Alloc, class ForwardIt, class Size, 
//	typename T = typename std::iterator_traits<ForwardIt>::value_type,
//	typename = std::enable_if_t<std::is_trivially_default_constructible<T>{}> 
//>
//ForwardIt uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n, void* = 0){
//	return uninitialized_default_construct_n(a, to_address(first), n);
//}

template<class AA> class is_allocator{
	static std::false_type aux(...     );
	template<
		class A, 
		class P = typename A::pointer, class S = typename A::size_type,
		typename = decltype(
			std::declval<A const&>()==A{std::declval<A const&>()},
			std::declval<A&>().deallocate(P{std::declval<A&>().allocate(std::declval<S>())}, std::declval<S>())
		)
	>
	static std::true_type  aux(A const&);
public:
	static bool const value = decltype(aux(std::declval<AA>()))::value;
	constexpr operator bool() const{return value;}
};

}}

namespace boost{
namespace multi{

/*
template<class Alloc, class It1, class T, class Ptr, class Ref>
void uninitialized_copy(Alloc& a, It first, It last, multi::array_iterator<T, 1, Ptr, Ref> dest){
	while(first != last){
		using std::begin; using std::end;
		uninitialized_copy(a, begin(*first), end(*first), begin(*dest)); // to make it work with T[][]
		++first;
		++dest;
	}
	return dest;
}*/

#if 0
template<dimensionality_type N> struct recursive_uninitialized_copy_aux;

template<dimensionality_type N, class Alloc, class InputIt, class ForwardIt>
ForwardIt recursive_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
	return recursive_uninitialized_copy_aux<N>::call(a, first, last, dest);
}

template<dimensionality_type N>
struct recursive_uninitialized_copy_aux{
	template<class Alloc, class InputIt, class ForwardIt>
	static auto call(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			using std::begin; using std::end;
			recursive_uninitialized_copy<N-1>(a, begin(*first), end(*first), begin(*dest)); // to make it work with T[][]
			++first;
			++dest;
		}
		return dest;
	}
};

template<>
struct recursive_uninitialized_copy_aux<1>{
	template<class Alloc, class InputIt, class ForwardIt>
	static auto call(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		return uninitialized_copy(a, first, last, dest);
	}
};
#endif

template<dimensionality_type N, class InputIt, class ForwardIt>
auto uninitialized_copy(InputIt first, InputIt last, ForwardIt dest){
	while(first!=last){
		uninitialized_copy<N-1>(begin(*first), end(*first), begin(*dest));
		++first;
		++dest;
	}
	return dest;
}

template<dimensionality_type N> struct recursive_fill_aux;

template<dimensionality_type N, class Out, class T>
void recursive_fill(Out f, Out l, T const& value){
	return recursive_fill_aux<N>::call(f, l, value);
}

template<dimensionality_type N>
struct recursive_fill_aux{
	template<class Out, class T>
	static auto call(Out first, Out last, T const& value){
		using std::begin; using std::end;
		for(; first != last; ++first)
			recursive_fill<N-1>(begin(*first), end(*first), value); // (*first).begin() instead of first->begin() to make it work with T[][]
	}
};

template<> struct recursive_fill_aux<1>{
	template<class O, class T>  static auto call(O f, O l, T const& v){
		using std::fill; return fill(f, l, v);
	}
};

#if 0
template<dimensionality_type N> struct recursive_uninitialized_fill_aux;

template<dimensionality_type N, class Alloc, class Out, class T>
void recursive_uninitialized_fill(Alloc& a, Out f, Out l, T const& value){
	return recursive_uninitialized_fill_aux<N>::call(a, f, l, value);
}

template<dimensionality_type N>
struct uninitialized_fill_aux{
	template<class Alloc, class Out, class T>
	static auto call(Alloc& a, Out first, Out last, T const& v){
		using std::begin; using std::end;
		for(; first != last; ++first)
			recursive_uninitialized_fill<N-1>(a, begin(*first), end(*first), v); // (*first).begin() instead of first->begin() to make it work with T[][]
	}
};

template<> struct uninitialized_fill_aux<1>{template<class Alloc, class O, class T> 
	static auto call(Alloc& a, O f, O l, T const& v){return uninitialized_fill(a, f, l, v);}
};
#endif

}}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MULTI_DETAIL_MEMORY

#include<vector>

namespace multi = boost::multi;

template<class T> void what(T&&);

int main(){

	static_assert(multi::is_allocator<std::allocator<double>>{}, "!");
	static_assert(not multi::is_allocator<double>{}, "!");
	static_assert(not multi::is_allocator<std::vector<double>>{}, "!");

	{
		double* p = nullptr;
		auto a = multi::default_allocator_of(p);
		static_assert(std::is_same<decltype(a), std::allocator<double>>{}, "!");
	//	what(typename std::iterator_traits<double*>::value_type{});
	}
#if 0
	{
		std::vector<double>::iterator it;
		auto a = multi::default_allocator_of(it);
	//	what(typename std::iterator_traits<decltype(it)>::value_type{});
		static_assert(std::is_same<decltype(a), std::allocator<double>>{}, "!");
	}
#endif

}
#endif
#endif

