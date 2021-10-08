// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2020-2021

#ifndef MULTI_POINTER_TRAITS_HPP
#define MULTI_POINTER_TRAITS_HPP

#include<memory>
#include<type_traits> // std::conditional_t

namespace boost{
namespace multi{

template<std::size_t I> struct Priority : std::conditional_t<I==0, std::true_type, struct Priority<I-1>>{}; 

template<class Pointer>  auto dat_aux(Priority<0>, Pointer ) -> std::allocator<typename std::iterator_traits<Pointer>::value_type>;
template<class T>        auto dat_aux(Priority<1>, T*      ) -> std::allocator<typename std::iterator_traits<T*>::value_type>;
template<class FancyPtr> auto dat_aux(Priority<2>, FancyPtr) -> typename FancyPtr::default_allocator_type;

template<class Pointer>
struct pointer_traits/*, typename Pointer::default_allocator_type>*/ : std::pointer_traits<Pointer>{
	using default_allocator_type = decltype(dat_aux(Priority<2>{}, std::declval<Pointer>()));
};

} // end namespace multi
} // end namespace boost
#endif

