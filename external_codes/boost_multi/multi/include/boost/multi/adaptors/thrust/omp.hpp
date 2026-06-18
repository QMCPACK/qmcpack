// Copyright 2018-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP

#include "boost/multi/array.hpp"

#include <thrust/system/omp/memory.h> // for ::thrust::omp::allocator
#include <type_traits>

namespace boost::multi::thrust::omp {
	template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::omp::allocator<T>>;
}  // end namespace boost::multi::thrust::omp

namespace thrust {

// template<class It> struct iterator_system;  // not needed in cuda 12.0, doesn't work on cuda 12.5

namespace detail {

//     template<class MultiArrayIterator>  // ::boost::multi::dimensionality_type D, class Pointer, bool IsConst, bool IsMove, typename Stride>
// struct iterator_system_impl<MultiArrayIterator, std::void_t<typename MultiArrayIterator::element_ptr> > {
//  using type = typename ::thrust::iterator_system<typename MultiArrayIterator::element_ptr>::type;
// };

}  // end namespace detail

template<class T, ::boost::multi::dimensionality_type D, class Pointer, bool IsConst, bool IsMove, typename Stride>
struct iterator_system<::boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride> > {
    using type = typename ::thrust::iterator_system<typename boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride>::element_ptr>::type;
};

template<typename Pointer, class LayoutType>
struct iterator_system<::boost::multi::elements_iterator_t<Pointer, LayoutType> > {  // TODO(correaa) might need changes for IsConst templating
	using type = typename ::thrust::iterator_system<typename ::boost::multi::elements_iterator_t<Pointer, LayoutType>::pointer>::type;
};

// namespace detail {
// template<class T1, class T2, class LO>
// struct pointer_traits<
//  boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, thrust::pointer<T2, thrust::cuda_cub::tag, thrust::tagged_reference<T2, thrust::cuda_cub::tag>, thrust::use_default>, LO>, 
//      LO
//  >
// >
// {
//  using Ptr = boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, thrust::pointer<T2, thrust::cuda_cub::tag, thrust::tagged_reference<T2, thrust::cuda_cub::tag>, thrust::use_default>, LO>, 
//      LO
//  >;
//  using pointer = Ptr;
//  using reference = thrust::tagged_reference<T2, thrust::cuda_cub::tag>;
//   typedef typename pointer_element<Ptr>::type    element_type;
//   typedef typename pointer_difference<Ptr>::type difference_type;

//   template<typename U>
//     struct rebind
//   {
//     typedef typename rebind_pointer<Ptr,U>::type other;
//   };

// //  __host__ __device__
// //   inline static pointer pointer_to(typename pointer_traits_detail::pointer_to_param<element_type>::type r)
// //   {
// //     // XXX this is supposed to be pointer::pointer_to(&r); (i.e., call a static member function of pointer called pointer_to)
// //     //     assume that pointer has a constructor from raw pointer instead

// //     return pointer(&r);
// //   }

//   // thrust additions follow
//   //typedef typename pointer_raw_pointer<Ptr>::type raw_pointer;
//  using raw_pointer = boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, T2*, LO>, 
//      LO
//  >;

//   __host__ __device__
//   inline static raw_pointer get(pointer ptr)
//   {
//  return reinterpret_cast<raw_pointer&>(ptr); //     return ptr.get();
//   }
// };
// }

}  // end namespace thrust

#endif  // BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP
