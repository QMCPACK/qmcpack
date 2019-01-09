#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra `#-Wfatal-errors` -D_TEST_BOOST_MULTI_DETAIL_MEMORY $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "layout.hpp"

#include<memory>

namespace boost{
namespace multi{

#if __cplusplus < 201703L
// https://en.cppreference.com/w/cpp/memory/destroy_at
template<class T> void destroy_at(T* p){p->~T();}
// https://en.cppreference.com/w/cpp/memory/destroy_n
template<class ForwardIt, class Size>
ForwardIt destroy_n(ForwardIt first, Size n){
	for(; n > 0; (void)++first, --n)
		destroy_at(std::addressof(*first));
	return first;
}
// https://en.cppreference.com/w/cpp/memory/destroy
template<class ForwardIt>
void destroy(ForwardIt first, ForwardIt last){
  for(; first != last; ++first) destroy_at(std::addressof(*first));
}
#else
using std::destroy_at;
using std::destroy_n;
using std::destroy;
#endif

#if __cplusplus < 201703L
template<class ForwardIt, class Size>
ForwardIt uninitialized_default_construct_n(ForwardIt first, Size n){
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    ForwardIt current = first;
    try{
        for(; n > 0 ; (void) ++current, --n) ::new (static_cast<void*>(std::addressof(*current))) T;
        return current;
    }catch(...){
        destroy(first, current); throw;
    }
}
template<class ForwardIt, class Size>
ForwardIt uninitialized_value_construct_n(ForwardIt first, Size n){
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    ForwardIt current = first;
    try{
        for(; n > 0 ; (void) ++current, --n) ::new (static_cast<void*>(std::addressof(*current))) T();
        return current;
    }catch(...){
        destroy(first, current); throw;
    }
}
#else
using std::uninitialized_default_construct_n;
using std::uninitialized_value_construct_n;
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
			uninitialized_copy<N-1>(first->begin(), first->end(), dest->begin());
			++first;
			++dest;
		}
		return dest;		
	}
};

template<>
struct uninitialized_copy_aux<1u>{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		using std::uninitialized_copy;
		return uninitialized_copy(first, last, dest);
	}
};

template<class Pointer>
struct pointer_traits : std::pointer_traits<Pointer>{
//	using default_allocator_type = 
//		decltype(default_allocator_of(std::declval<Pointer const&>()));
	using allocator_type = std::allocator<std::decay_t<typename pointer_traits::element_type>>;
	template<class P2>
	static allocator_type allocator_of(P2 const&){return {};}
};

template<class T>
auto default_allocator_of(T*){
	return std::allocator<std::decay_t<typename std::pointer_traits<T*>::element_type>>{};
}

template<class Ptr>
auto default_allocator_of(Ptr){
	return std::allocator<std::decay_t<typename std::pointer_traits<Ptr>::element_type>>{};
}

}}

#if _TEST_BOOST_MULTI_DETAIL_MEMORY

namespace multi = boost::multi;

int main(){

}
#endif
#endif

#if 0
template<dimensionality_type N> struct uninitialized_copy_from_initializer_list_aux;

template<dimensionality_type N, class InputIt, class ForwardIt>
ForwardIt uninitialized_copy_from_initializer_list(InputIt first, InputIt last, ForwardIt dest){
	return uninitialized_copy_from_initializer_list_aux<N>::call(first, last, dest);
}

template<dimensionality_type N>
struct uninitialized_copy_from_il;

template<>
struct uninitialized_copy_from_il<1u>{
	template<class InpIt, class FwdIt>
	static auto call(InpIt first, InpIt last, FwdIt dest){
		while(first != last){
		//	construct_from_il(std::addressof(*dest), *first);
			using T = typename std::iterator_traits<FwdIt>::value_type;
			::new(static_cast<void*>(std::addressof(*dest))) T(*first);
		    ++first;
		    ++dest;
		}
		return dest;
	}
};

template<dimensionality_type N>
struct uninitialized_copy_from_il{
	template<class InIt, class FwdIt>
	static auto call(InIt first, InIt last, FwdIt dest){
		while(first != last){
			uninitialized_copy_from_il<N-1>::call(
				first->begin(), first->end(), dest->begin()
			);
			++first;
			++dest;
		}
		return dest;
	}
};

template<dimensionality_type N>
struct uninitialized_copy_from_initializer_list_aux{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			uninitialized_copy_from_initializer_list<N-1>(
				first->begin(), first->end(), dest->begin()
			);
			++first;
			++dest;
		}
		return dest;		
	}
};

template<>
struct uninitialized_copy_from_initializer_list_aux<1u>{
	template<class InputIt, class ForwardIt>
	static auto call(InputIt first, InputIt last, ForwardIt dest){
		while(first != last){
			construct_from_initializer_list(std::addressof(*dest), *first);
		//	using T = typename std::iterator_traits<ForwardIt>::value_type;
		//	::new(static_cast<void*>(std::addressof(*dest))) T(*first);
		    ++first;
		    ++dest;
		}
		return dest;
	}
};
#endif

