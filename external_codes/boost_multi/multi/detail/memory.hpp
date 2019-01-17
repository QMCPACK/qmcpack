#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra `#-Wfatal-errors` -D_TEST_BOOST_MULTI_DETAIL_MEMORY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "layout.hpp"

#include<memory>

namespace boost{
namespace multi{

#if __cplusplus >= 201703L
using std::uninitialized_default_construct_n;
using std::uninitialized_value_construct_n;
using std::destroy_at;
using std::destroy;
using std::destroy_n;
#else
// https://en.cppreference.com/w/cpp/memory/destroy_at
template<class T> void destroy_at(T* p){p->~T();}
// https://en.cppreference.com/w/cpp/memory/destroy_n
template<class ForwardIt, class Size>
ForwardIt destroy_n(ForwardIt first, Size n){
	for(; n > 0; (void) ++first, --n) destroy_at(std::addressof(*first));
	return first;
}
//https://en.cppreference.com/w/cpp/memory/destroy
template<class ForwardIt>
void destroy(ForwardIt first, ForwardIt last){
	for(; first != last; ++first) destroy_at(std::addressof(*first));
}

template<class ForwardIt, class Size>
ForwardIt uninitialized_default_construct_n(ForwardIt first, Size n){
	using T = typename std::iterator_traits<ForwardIt>::value_type;
	ForwardIt current = first;
	try{
		for(; n > 0; (void) ++current, --n)
			::new (static_cast<void*>(std::addressof(*current))) T;
		return current;
	}catch(...){destroy(first, current); throw;}
}
template<class ForwardIt, class Size>
ForwardIt uninitialized_value_construct_n(ForwardIt first, Size n){
	using T = typename std::iterator_traits<ForwardIt>::value_type;
	ForwardIt current = first;
	try{
		for(; n > 0; (void) ++current, --n)
			::new (static_cast<void*>(std::addressof(*current))) T();
		return current;
    }catch(...){destroy(first, current); throw;}
}
#endif

template<dimensionality_type N> struct uninitialized_copy_aux;

template<dimensionality_type N, class InputIt, class ForwardIt>
ForwardIt uninitialized_copy(InputIt first, InputIt last, ForwardIt dest){
	return uninitialized_copy_aux<N>::call(first, last, dest);
}

template<dimensionality_type N>
struct uninitialized_copy_aux{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			using std::begin; using std::end;
			uninitialized_copy<N-1>(begin(*first), end(*first), begin(*dest)); // to make it work with T[][]
		//	uninitialized_copy<N-1>(first->begin(), first->end(), dest->begin());
			++first;
			++dest;
		}
		return dest;
	}
};

template<>
struct uninitialized_copy_aux<1>{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		using std::uninitialized_copy;
		return uninitialized_copy(first, last, dest);
	}
};

template<dimensionality_type N> struct fill_aux;

template<dimensionality_type N, class Out, class T>
void fill(Out f, Out l, T const& value){return fill_aux<N>::call(f, l, value);}

template<dimensionality_type N>
struct fill_aux{
	template<class Out, class T>
	static auto call(Out first, Out last, T const& value){
		using std::begin; using std::end;
		for(; first != last; ++first)
			fill<N-1>(begin(*first), end(*first), value); // (*first).begin() instead of first->begin() to make it work with T[][]
	}
};

template<> struct fill_aux<1>{
	template<class O, class T> 
	static auto call(O f, O l, T const& v){using std::fill; return fill(f, l, v);}
};

template<dimensionality_type N> struct uninitialized_fill_aux;

template<dimensionality_type N, class Out, class T>
void uninitialized_fill(Out f, Out l, T const& value){
	return uninitialized_fill_aux<N>::call(f, l, value);
}

template<dimensionality_type N>
struct uninitialized_fill_aux{
	template<class Out, class T>
	static auto call(Out first, Out last, T const& value){
		using std::begin; using std::end;
		for(; first != last; ++first)
			uninitialized_fill<N-1>(begin(*first), end(*first), value); // (*first).begin() instead of first->begin() to make it work with T[][]
	}
};

template<> struct uninitialized_fill_aux<1>{
	template<class O, class T> 
	static auto call(O f, O l, T const& v){using std::uninitialized_fill; return uninitialized_fill(f, l, v);}
};

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

