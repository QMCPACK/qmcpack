// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_POINTER_TRAITS_HPP
#define BOOST_MULTI_DETAIL_POINTER_TRAITS_HPP
// #pragma once

#include <cstddef>      // for size_t
#include <iterator>     // for iterator_traits
#include <memory>       // for allocator, pointer_traits
#include <type_traits>  // for conditional_t, declval, true_type

namespace boost::multi {

namespace detail {
template<std::size_t N> struct priority_me : std::conditional_t<N == 0, std::true_type, priority_me<N-1>>{};

template<class Pointer>  auto dat_aux(priority_me<0>, Pointer ) -> std::allocator<typename std::iterator_traits<Pointer>::value_type>;
template<class T>        auto dat_aux(priority_me<1>, T*      ) -> std::allocator<typename std::iterator_traits<T*>::value_type>;
template<class FancyPtr> auto dat_aux(priority_me<2>, FancyPtr) -> typename FancyPtr::default_allocator_type;
}  // end namespace detail

template<class Pointer>
struct pointer_traits/*, typename Pointer::default_allocator_type>*/ : std::pointer_traits<Pointer>{
	using default_allocator_type = decltype(detail::dat_aux(detail::priority_me<2>(), std::declval<Pointer>()));
};

}  // end namespace boost::multi
#endif  // BOOST_MULTI_DETAIL_POINTER_TRAITS_HPP
