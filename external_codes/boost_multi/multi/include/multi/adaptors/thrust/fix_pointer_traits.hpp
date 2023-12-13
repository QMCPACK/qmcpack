// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2021-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_THRUST_FIX_POINTER_TRAITS_HPP_
#define MULTI_ADAPTORS_THRUST_FIX_POINTER_TRAITS_HPP_
#pragma once

#include <thrust/memory.h>

// #include <thrust/detail/type_traits/pointer_traits.h>

// begin of nvcc trhust 11.5 workaround : https://github.com/NVIDIA/thrust/issues/1629
namespace thrust {

template<typename Element, typename Tag, typename Reference, typename Derived> class pointer;
template<class T> struct pointer_traits;

}  // end namespace thrust

namespace std {

template<class... As> struct pointer_traits<thrust::pointer<As...>>
: thrust::detail::pointer_traits<thrust::pointer<As...>> {
	template<class T>
	using rebind = typename thrust::detail::pointer_traits<thrust::pointer<As...>>::template rebind<T>::other;
};

}  // end namespace std
// end of nvcc thrust 11.5 workaround

#endif
