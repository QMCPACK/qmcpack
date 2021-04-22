// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2021

#pragma once

#include "../../array.hpp"

#include "./thrust/cuda/managed.hpp"

#include <thrust/device_allocator.h>
#include <thrust/system/cuda/memory.h> // ::thrust::cuda::allocator



namespace boost{
namespace multi{
namespace thrust{

template<class T, multi::dimensionality_type D> using device_array = multi::array<T, D, ::thrust::device_allocator<T>>;
template<class T, multi::dimensionality_type D> using host_array   = multi::array<T, D                               >;

namespace device{

template<class T, multi::dimensionality_type D> using array = device_array<T, D>;

}

namespace host{

template<class T, multi::dimensionality_type D> using array = host_array<T, D>;

}

namespace cuda{

template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::cuda::allocator<T>>;

namespace managed{

template<class T, multi::dimensionality_type D> using array = multi::array<T, D, boost::multi::thrust::cuda::managed::allocator<T>>;

}

}

}}}

