// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2020-2021

#ifndef MULTI_MEMORY_POINTER_TRAITS_HPP
#define MULTI_MEMORY_POINTER_TRAITS_HPP

#include<memory>
#include<type_traits> // std::conditional_t

namespace boost::multi {

template<std::size_t I> struct priority_me : std::conditional_t<I==0, std::true_type, struct priority_me<I-1>>{};

template<class Pointer>  auto dat_aux(priority_me<0>, Pointer ) -> std::allocator<typename std::iterator_traits<Pointer>::value_type>;
template<class T>        auto dat_aux(priority_me<1>, T*      ) -> std::allocator<typename std::iterator_traits<T*>::value_type>;
template<class FancyPtr> auto dat_aux(priority_me<2>, FancyPtr) -> typename FancyPtr::default_allocator_type;

template<class Pointer>
struct pointer_traits/*, typename Pointer::default_allocator_type>*/ : std::pointer_traits<Pointer>{
	using default_allocator_type = decltype(dat_aux(priority_me<2>{}, std::declval<Pointer>()));
};

} // end namespace boost::multi
#endif
