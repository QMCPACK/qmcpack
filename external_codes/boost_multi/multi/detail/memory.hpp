#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++2a -Wall -Wextra `#-Wfatal-errors` -D_TEST_BOOST_MULTI_DETAIL_MEMORY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "layout.hpp"
#include "../utility.hpp"
#include<memory>


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

template<class T, typename = decltype(std::pointer_traits<T>::to_address(std::declval<T const&>()))> 
                  auto use_address_aux(T const& p)->std::true_type;
template<class T> auto use_address_aux(...       )->std::false_type;

template<class T> struct use_address : decltype(use_address_aux<T>()){};
 
template<class T>
constexpr T* to_address(T* p) noexcept
{
    static_assert(!std::is_function<T>{}, "!");
    return p;
}

template<class T>
auto to_address_aux(const T& p, std::true_type) noexcept{
   return std::pointer_traits<T>::to_address(p);
}

template<class T>
auto to_address_aux(const T& p, std::false_type) noexcept
//->decltype(to_address(p.operator->()))
{
	return to_address(p.operator->());}
 
template<class T>
auto to_address(const T& p) noexcept{
	return to_address_aux(p, use_address<T>{});
//    if constexpr (use_address<T>::value) {
//       return std::pointer_traits<T>::to_address(p);
//    } else {
//       return memory::to_address(p.operator->());
//   }
}


}

using memory::allocator_traits;
using memory::to_address;
////////////////////////////////////////////////////////////////////////////////
//template<class Ptr> auto to_address(Ptr const& p) noexcept
//->decltype(p.operator->()){
//	return p.operator->();}

//template<class T> constexpr T* to_address(T* p) noexcept{return p;}

//https://en.cppreference.com/w/cpp/memory/destroy
template<class Alloc, class ForwardIt, typename = std::enable_if_t<!(has_rank<ForwardIt>{})>>//, typename = std::enable_if_t<typename ForwardIt::rank{} == 1> >//, typename AT = typename std::allocator_traits<Alloc> >
void destroy(Alloc& a, ForwardIt f, ForwardIt l){
	//	using multi::to_address;
	for(; f != l; ++f) allocator_traits<Alloc>::destroy(a, to_address(f));
	//	 a.destroy(to_address(first)); //	AT::destroy(a, to_address(first)); //	AT::destroy(a, addressof(*first)); // a.destroy(addressof(*first));
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
ForwardIt uninitialized_copy_n(Alloc& a, InputIt f, Size n, ForwardIt d){
	ForwardIt c = d;
//	using std::addressof;
	try{
		for(; n > 0; ++f, ++c, --n)
			a.construct(to_address(c), *f);
		//	AT::construct(a, to_address(c), *f);
		//	AT::construct(a, addressof(*c), *f);
		//	a.construct(addressof(*c), *f);
		return c;
	}catch(...){destroy(a, d, c); throw;}
}

template<class Alloc, class ForwardIt, class Size, class T>//, typename AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_fill_n(Alloc& a, ForwardIt first, Size n, const T& v){
	ForwardIt current = first; // using std::to_address;
	try{
		for(; n > 0; ++current, --n)
			allocator_traits<Alloc>::construct(a, to_address(current), v);
		//	a.construct(to_address(current), v); //	AT::construct(a, to_address(current), v); //	AT::construct(a, addressof(*current), v); //a.construct(addressof(*current), v);
		return current;
	}catch(...){destroy(a, first, current); throw;}
}

template<class Alloc, class ForwardIt, class Size>//, class AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_default_construct_n(Alloc& a, ForwardIt first, Size n){
	ForwardIt current = first;
	try{
		for(; n > 0; ++current, --n)
			allocator_traits<Alloc>::construct(a, to_address(current));
		//	a.construct(to_address(current));
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

template<class Alloc, class InputIt, class MMIt, 
	typename = std::enable_if_t<has_rank<MMIt>{}>,
	typename = std::enable_if_t<!std::is_trivially_copyable<typename MMIt::element>{}>, 
	typename = std::enable_if_t<typename MMIt::rank{} != 1>
>
MMIt uninitialized_copy(Alloc& a, InputIt f, InputIt l, MMIt d){
	MMIt c = d; using std::begin; using std::end;
	try{
		while(f!= l) uninitialized_copy(a, begin(*f), end(*f), begin(*c)), ++f, ++c; // to make it work with T[][]
		return c;
	}catch(...){destroy(a, d, c); throw;}
}

template<class Alloc, class In, class MIt, class=std::enable_if_t<std::is_trivially_copyable<typename MIt::element>{}>>
MIt uninitialized_copy(Alloc&, In f, In l, MIt d){using std::copy; return copy(f, l, d);}

template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<has_rank<MIt>{}>, typename = std::enable_if_t<typename MIt::rank{}==1>, class=std::enable_if_t<!std::is_trivially_copyable<typename MIt::element>{}> >
MIt uninitialized_copy(Alloc& a, InputIt f, InputIt l, MIt const& d){
	MIt current = d; // using multi::to_address;
	try{
		for(; f != l; ++f, ++current) a.construct(to_address(current), *f);
		return current;
	}catch(...){destroy(a, d, current); throw;}
}

template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<!has_rank<MIt>{}> >
MIt uninitialized_copy(Alloc& a, InputIt f, InputIt l, MIt d, double* = 0){
	MIt current = d;
//	using multi::to_address;
	try{
		for(; f != l; ++f, ++current) a.construct(to_address(current), *f);
		return current;
	}catch(...){destroy(a, d, current); throw;}
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
template<class Alloc, class ForwardIt, class Size>
ForwardIt destroy_n(Alloc& a, ForwardIt first, Size n);
template<class Alloc, class InputIt, class Size, class ForwardIt>
ForwardIt uninitialized_copy_n(Alloc& a, InputIt f, Size n, ForwardIt d);
template<class Alloc, class ForwardIt, class Size, class T>
ForwardIt uninitialized_fill_n(Alloc& a, ForwardIt first, Size n, const T& v);
template<class Alloc, class ForwardIt, class Size>
ForwardIt uninitialized_default_construct_n(Alloc& a, ForwardIt first, Size n);
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

template<class Alloc, class ForwardIt, class Size, 
	typename T = typename std::iterator_traits<ForwardIt>::value_type,
	typename = std::enable_if_t<std::is_trivially_default_constructible<T>{}> 
>
ForwardIt uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n, void* = 0){
	return uninitialized_default_construct_n(a, first, n);
}

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

template<class T, typename = decltype(std::declval<T const&>().default_allocator())>
std::true_type           has_default_allocator_aux(T const&);
std::false_type          has_default_allocator_aux(...);
template<class T> struct has_default_allocator : decltype(has_default_allocator_aux(std::declval<T>())){};

template<class Pointer, class T = void> struct pointer_traits;

template<class P>
struct pointer_traits<P, std::enable_if_t<has_default_allocator<P>{}> > : std::pointer_traits<P>{
	static auto default_allocator_of(typename pointer_traits::pointer const& p){return p.default_allocator();}
	using default_allocator_type = decltype(std::declval<P const&>().default_allocator());
};
template<class P>
struct pointer_traits<P, std::enable_if_t<not has_default_allocator<P>{}> > : std::pointer_traits<P>{
	using default_allocator_type = std::allocator<std::decay_t<typename pointer_traits::element_type> >;
	static default_allocator_type default_allocator_of(typename pointer_traits::pointer const&){return {};}
};

template<class P>
typename pointer_traits<P>::default_allocator_type 
default_allocator_of(P const& p){return pointer_traits<P>::default_allocator_of(p);}

}}

#if _TEST_BOOST_MULTI_DETAIL_MEMORY

namespace multi = boost::multi;

int main(){

	double* p;
	auto a = multi::default_allocator_of(p);
	static_assert(std::is_same<decltype(a), std::allocator<double>>{}, "!");

}
#endif
#endif

