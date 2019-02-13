#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra `#-Wfatal-errors` -D_TEST_BOOST_MULTI_DETAIL_MEMORY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "layout.hpp"

#include<memory>

namespace boost{
namespace multi{

// https://en.cppreference.com/w/cpp/memory/destroy_at
template<class Alloc, class T, typename AT = std::allocator_traits<Alloc> > 
void destroy_at(Alloc& a, T* p){AT::destroy(a, p);}//p->~T();}

//https://en.cppreference.com/w/cpp/memory/destroy
template<class Alloc, class ForwardIt, typename AT = typename std::allocator_traits<Alloc> >
void destroy(Alloc& a, ForwardIt first, ForwardIt last){
	for(; first != last; ++first) 
		AT::destroy(a, std::addressof(*first));//destroy_at(std::addressof(*first));
}

// https://en.cppreference.com/w/cpp/memory/destroy_n
template<class Alloc, class ForwardIt, class Size, typename AT = typename std::allocator_traits<Alloc> >
ForwardIt destroy_n(Alloc& a, ForwardIt first, Size n){
	for(; n > 0; ++first, --n)
		AT::destroy(a, std::addressof(*first));//std::destroy_at(std::addressof(*first));
	return first;
}

template<class Alloc, class InputIt, class ForwardIt, typename AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt d){
	ForwardIt current = d;
	try{
		for(; first != last; ++first, ++current)
			AT::construct(a, std::addressof(*current), *first);
		return current;
	}catch(...){
		for(; d != current; ++d) AT::destroy(a, std::addressof(*d));
		throw;
	}
}

template<class Alloc, class InputIt, class Size, class ForwardIt, typename AT = std::allocator_traits<Alloc> >
ForwardIt uninitialized_copy_n(Alloc& a, InputIt first, Size count, ForwardIt d){
	ForwardIt current = d;
	try{
		for(; count > 0; ++first, (void) ++current, --count)
			AT::construct(a, std::addressof(*current), *first); // ::new (static_cast<void*>(std::addressof(*current))) Value(*first);
	}catch(...){
		for(; d != current; ++d) AT::destroy(a, std::addressof(*d));
		throw;
	}
	return current;
}

template<class Alloc, class ForwardIt, class Size, class T, typename AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_fill_n(Alloc& a, ForwardIt first, Size count, const T& value){
	ForwardIt current = first;
	try{
		for(; count > 0; ++current, (void) --count)
			AT::construct(a, std::addressof(*current), value); // ::new (static_cast<void*>(std::addressof(*current))) Value(value);
		return current;
	}catch(...){
		for(; first != current; ++first) AT::destroy(a, std::addressof(*first));
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class AT = typename std::allocator_traits<Alloc> >
ForwardIt uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n){
	using T = typename std::iterator_traits<ForwardIt>::value_type;
	ForwardIt current = first;
	try{
		for(; n > 0; (void) ++current, --n) 
			AT::construct(a, std::addressof(*current), T()); // ::new (static_cast<void*>(std::addressof(*current))) T();
		return current;
    }catch(...){destroy(a, first, current); throw;}
}

}}

namespace boost{
namespace multi{

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

