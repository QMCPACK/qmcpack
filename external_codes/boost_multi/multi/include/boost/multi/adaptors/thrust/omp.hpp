// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP

// CCCL's OMP backend (e.g. thrust/system/omp/detail/scan.h in CCCL 3.3+) uses
// `numeric_limits<T>::max()`; on Windows the <windows.h> `max(...)` function-macro mangles it
// ("illegal token '(' on right side of '::'"). Define NOMINMAX before any system/thrust include.
// Windows-only and non-overriding, so it is a no-op for the Linux CUDA 12/13 builds.
#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif

#include <boost/multi/array.hpp>

#include <thrust/system/omp/memory.h>  // for ::thrust::omp::allocator
#include <type_traits>

namespace boost::multi::thrust::omp {
	template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::omp::allocator<T>>;
}  // end namespace boost::multi::thrust::omp

namespace thrust {

template<class T, ::boost::multi::dimensionality_type D, class Pointer, bool IsConst, bool IsMove, typename Stride, class SubLayout>
struct iterator_system<::boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride, SubLayout> > {
    using type = typename ::thrust::iterator_system<typename boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride, SubLayout>::element_ptr>::type;
};

template<typename Pointer, class LayoutType>
struct iterator_system<::boost::multi::detail::elements_iterator_t<Pointer, LayoutType> > {  // TODO(correaa) might need changes for IsConst templating
	using type = typename ::thrust::iterator_system<typename ::boost::multi::detail::elements_iterator_t<Pointer, LayoutType>::pointer>::type;
};

}  // end namespace thrust

#endif  // BOOST_MULTI_ADAPTORS_THRUST_OMP_HPP
