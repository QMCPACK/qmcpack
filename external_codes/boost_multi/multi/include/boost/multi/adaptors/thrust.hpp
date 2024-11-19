// Copyright 2021-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_THRUST_HPP_
#define BOOST_MULTI_ADAPTORS_THRUST_HPP_
#pragma once

#include <boost/multi/array.hpp>

#include <thrust/device_allocator.h>

#include <thrust/universal_allocator.h>

#if !defined(MULTI_USE_HIP)
#include <thrust/system/cuda/memory.h> // for ::thrust::cuda::allocator
#else
#include <thrust/system/hip/memory.h>  // for ::thrust::hip::allocator
#endif
// #include <thrust/system/cuda/memory.h>  // ::thrust::cuda::allocator

// #include <thrust/detail/type_traits/pointer_traits.h>

// #include <utility>  // std::copy

#include <thrust/detail/pointer.h>                             // for pointer

#include <thrust/mr/allocator.h>                               // for allocator (ptr only), stateless_resource_allocator
#include <thrust/mr/memory_resource.h>                         // for memory_resource

#if !defined(MULTI_USE_HIP)
#include <thrust/system/cuda/detail/execution_policy.h>        // for tag
#include <thrust/system/cuda/memory_resource.h>                // for universal_memory_resource
#include <thrust/system/cuda/pointer.h>                        // for universal_pointer

#include <cuda_runtime_api.h>                                  // for cudaGetDevice, cudaMemPrefetchAsync, cudaPointerGetAttributes
#include <driver_types.h>                                      // for cudaErrorInvalidValue, cudaPointerAttributes, cudaSuccess, cudaErrorInvalidDevice, cudaMemoryTypeManaged
#else
// #include <thrust/system/hip/detail/execution_policy.h>        // for tag
// #include <thrust/system/hip/memory_resource.h>                // for universal_memory_resource
// #include <thrust/system/hip/pointer.h>                        // for universal_pointer

// #include <hip_runtime_api.h>                                  // for cudaGetDevice, cudaMemPrefetchAsync, cudaPointerGetAttributes
#endif

#include <boost/multi/adaptors/thrust/fix_pointer_traits.hpp>

#include<cassert>
#include <iterator>                                            // for iterator_traits
#include <memory>                                              // for allocator_traits, allocator, pointer_traits
// #include <thrust/iterator/detail/iterator_traits.inl>          // for iterator_system
#include <type_traits>                                         // for decay_t

// #include <boost/multi/adaptors/thrust/fix_copy.hpp>

// // begin of nvcc trhust 11.5 workaround : https://github.com/NVIDIA/thrust/issues/1629
// namespace thrust {

// template<typename Element, typename Tag, typename Reference, typename Derived> class pointer;
// template<class T> struct pointer_traits;

// }  // end namespace thrust

// namespace std {

// template<class... As> struct pointer_traits<thrust::pointer<As...>>
// : thrust::detail::pointer_traits<thrust::pointer<As...>> {
//  template<class T>
//  using rebind = typename thrust::detail::pointer_traits<thrust::pointer<As...>>::template rebind<T>::other;
// };

// }  // end namespace std
// // end of nvcc thrust 11.5 workaround

#if !defined(MULTI_USE_HIP)
#define HICUP cuda
#define HICUP_(NAME)  cuda ## NAME
#else
#define HICUP hip
#define HICUP_(NAME)  hip ## NAME
#endif

namespace boost::multi { template <class Alloc> struct allocator_traits; }

namespace boost::multi {

template<class T>
struct pointer_traits<::thrust::pointer<T, ::thrust::HICUP::tag, T&>> : std::pointer_traits<::thrust::pointer<T, ::thrust::HICUP::tag, T&>> {
	using default_allocator_type = ::thrust::universal_allocator<std::decay_t<T>>;
};

} // end namespace boost::multi

namespace boost::multi {

template<class TT>
struct allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>>
: std::allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>> {
 private:
	using Alloc = ::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>;
	using base = std::allocator_traits<Alloc>;

 public:
	using typename base::pointer;
	using typename base::size_type;
	using typename base::const_void_pointer;

	using base::allocate;
	[[nodiscard]] static constexpr auto allocate(Alloc& a, size_type n, const_void_pointer hint) -> pointer {
		auto ret = allocator_traits::allocate(a, n);
		if(not hint) {
			prefetch_to_device(ret, n*sizeof(TT), get_current_device());
			return ret;
		}
		prefetch_to_device(ret, n*sizeof(TT), get_device(hint));
		return ret;
	}

 private:
	using device_index = int;
	static auto get_current_device() -> device_index {
		int device;
		switch(HICUP_(GetDevice)(&device)) {
			case HICUP_(Success)          : break;
			case HICUP_(ErrorInvalidValue): assert(0);
			default: assert(0);
		}
		return device;
	}
	static void prefetch_to_device(const_void_pointer p, size_type byte_count, device_index d) {
		switch(HICUP_(MemPrefetchAsync)(raw_pointer_cast(p), byte_count, d)) {
			case HICUP_(Success)           : break;
			case HICUP_(ErrorInvalidValue) : assert(0); break;
			case HICUP_(ErrorInvalidDevice): assert(0); break;
			default: assert(0);
		}
	}

	static auto get_device(const_void_pointer p) -> device_index {
		#if defined(__HIPCC__)
		hipPointerAttribute_t attr{};
		#else  // #if defined(__NVCC__)
		cudaPointerAttributes attr{};
		#endif
		switch(HICUP_(PointerGetAttributes)(&attr, raw_pointer_cast(p))) {
			case HICUP_(Success): break;
			case HICUP_(ErrorInvalidDevice): assert(0); break;
			case HICUP_(ErrorInvalidValue): assert(0); break;
			default: assert(0);  // 71 enumeration values not handled in switch: 'hipErrorOutOfMemory', 'hipErrorNotInitialized', 'hipErrorDeinitialized'...
		}
		assert(attr.type == HICUP_(MemoryTypeManaged));
		return attr.device;
	}
};

}  // end namespace ::boost::multi

// this is important for algorithms to dispatch to the right thrust executor
namespace thrust {

// template<class It> struct iterator_system;  // not needed in cuda 12.0, doesn't work on cuda 12.5

template<class T, ::boost::multi::dimensionality_type D, class Pointer, bool IsConst>
struct iterator_system<::boost::multi::array_iterator<T, D, Pointer, IsConst>>{
	using type = typename ::thrust::iterator_system<typename boost::multi::array_iterator<T, D, Pointer, IsConst>::element_ptr>::type;
};

template<typename Pointer, class LayoutType>
struct iterator_system<::boost::multi::elements_iterator_t<Pointer, LayoutType>> {  // TODO(correaa) might need changes for IsConst templating
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

}  // end namespace ::thrust

namespace boost::multi {
namespace thrust {

// defines multi::thrust::device_array
// defines multi::thrust::host_array

template<typename T, multi::dimensionality_type D, class Alloc = ::thrust::device_allocator   <T>> using device_array    = multi::array<T, D, Alloc>;
template<typename T, multi::dimensionality_type D, class Alloc = ::thrust::universal_allocator<T>> using universal_array = multi::array<T, D, Alloc>;
template<typename T, multi::dimensionality_type D, class Alloc = std::allocator               <T>> using host_array      = multi::array<T, D, Alloc>;

// defines multi::thrust::device::array
// defines multi::thrust::host  ::array
namespace device    {template<class T, multi::dimensionality_type D> using array = device_array   <T, D>;}  // end namespace device
namespace universal {template<class T, multi::dimensionality_type D> using array = universal_array<T, D>;}  // end namespace universal
namespace host      {template<class T, multi::dimensionality_type D> using array = host_array     <T, D>;}  // end namespace host

// defines multi::thrust::cuda::array
// defines multi::thrust::cuda::managed::array
namespace cuda {
	template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::HICUP::allocator<T>>;

	// namespace managed {
	//  template<class T, multi::dimensionality_type D> using array = multi::array<T, D, boost::multi::thrust::hip::managed::allocator<T>>;
	// }  // end namespace managed
}  // end namespace cuda

namespace  mr {template<class T, multi::dimensionality_type D, class MR> using array = array<T, D, ::thrust::mr::allocator<T, MR>>;}
namespace pmr {
	template<class T, multi::dimensionality_type D, class Pointer> using array = mr::array<T, D, ::thrust::mr::memory_resource<Pointer>>;
	template<class T, multi::dimensionality_type D> using universal_array = pmr::array<T, D, ::thrust::universal_ptr<void>>;
}  // end namespace pmr

namespace cuda {

template<class T, multi::dimensionality_type D> using universal_array = multi::array<T, D, ::thrust::HICUP::universal_allocator<T> >;

namespace universal {
	template<class T, multi::dimensionality_type D> using array = multi::thrust::cuda::universal_array<T, D>;
}

namespace pmr {
	template<class T, multi::dimensionality_type D> using universal_array = ::boost::multi::thrust::pmr::array<T, D, ::thrust::HICUP::universal_pointer<void>>;
}  // end namespace pmr
}  // end namespace cuda


}  // end namespace thrust
}  // end namespace boost::multi

namespace boost::multi {

template<class Q, class R>
constexpr auto default_allocator_of(::thrust::pointer<Q, ::thrust::HICUP::tag, Q&> /*unused*/) {
	return ::thrust::HICUP::universal_allocator<typename std::iterator_traits<::thrust::pointer<Q, ::thrust::HICUP::tag, Q&>>::value_type>{};
}

}

#undef HICUP
#undef HICUP_

#endif
