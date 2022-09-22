// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2021-2022 Alfredo A. Correa

#pragma once

#include "../array.hpp"

#include "./thrust/cuda/managed.hpp"

#include <thrust/universal_allocator.h>

#include <thrust/device_allocator.h>
#include <thrust/system/cuda/memory.h>  // ::thrust::cuda::allocator

//#include <thrust/system/cuda/experimental/pinned_allocator.h>
//#include <thrust/host_vector.h>

#include <thrust/detail/type_traits/pointer_traits.h>

#include <utility>  // std::copy

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
// end of nvcc trhust 11.5 workaround

namespace boost::multi {

template<class T>
struct pointer_traits<::thrust::pointer<T, ::thrust::cuda_cub::tag, T&>> : std::pointer_traits<::thrust::pointer<T, ::thrust::cuda_cub::tag, T&>> {
	using default_allocator_type = ::thrust::universal_allocator<std::decay_t<T>>;
};

} // end namespace boost::multi

namespace boost::multi {

template<class TT>
struct allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::cuda::universal_memory_resource>>
: std::allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::cuda::universal_memory_resource>> {
 private:
	using Alloc = ::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::cuda::universal_memory_resource>;
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
		switch(cudaGetDevice(&device)) {
			case cudaSuccess          : break;
			case cudaErrorInvalidValue: assert(0);
		}
		return device;
	}
	static void prefetch_to_device(const_void_pointer p, size_type byte_count, device_index d) {
		switch(cudaMemPrefetchAsync(raw_pointer_cast(p), byte_count, d)) {
			case cudaSuccess           : break;
			case cudaErrorInvalidValue : assert(0); break;
			case cudaErrorInvalidDevice: assert(0); break;
		}
	}

	static auto get_device(const_void_pointer p) -> device_index {
		cudaPointerAttributes attr{};
		switch(cudaPointerGetAttributes(&attr, raw_pointer_cast(p))) {
			case cudaSuccess: break;
			case cudaErrorInvalidDevice: assert(0); break;
			case cudaErrorInvalidValue: assert(0); break;
		}
		assert(attr.type == cudaMemoryTypeManaged);
		return attr.device;
	}
};

}

// this is important for algorithms to dispatch to the right thrust executor
namespace thrust {

template<class It> struct iterator_system;

template<class T, boost::multi::dimensionality_type D, class Pointer>
struct iterator_system<boost::multi::array_iterator<T, D, Pointer>>{
	using type = typename thrust::iterator_system<typename boost::multi::array_iterator<T, D, Pointer>::element_ptr>::type;
};

template<typename Pointer, class LayoutType>
struct iterator_system<boost::multi::elements_iterator_t<Pointer, LayoutType>> {
	using type = typename thrust::iterator_system<typename boost::multi::elements_iterator_t<Pointer, LayoutType>::pointer>::type;
};

}

namespace boost::multi {
namespace thrust {

// defines multi::thrust::device_array
// defines multi::thrust::host_array

template<class T, multi::dimensionality_type D> using device_array    = multi::array<T, D, ::thrust::device_allocator   <T>>;
template<class T, multi::dimensionality_type D> using universal_array = multi::array<T, D, ::thrust::universal_allocator<T>>;
template<class T, multi::dimensionality_type D> using host_array      = multi::array<T, D                                  >;

// defines multi::thrust::device::array
// defines multi::thrust::host  ::array
namespace device    {template<class T, multi::dimensionality_type D> using array = device_array   <T, D>;}  // end namespace device
namespace universal {template<class T, multi::dimensionality_type D> using array = universal_array<T, D>;}  // end namespace universal
namespace host      {template<class T, multi::dimensionality_type D> using array = host_array     <T, D>;}  // end namespace host

// defines multi::thrust::cuda::array
// defines multi::thrust::cuda::managed::array
namespace cuda {
	template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::cuda::allocator<T>>;

	namespace managed {
		template<class T, multi::dimensionality_type D> using array = multi::array<T, D, boost::multi::thrust::cuda::managed::allocator<T>>;
	}  // end namespace managed
}  // end namespace cuda

}  // end namespace thrust
}  // end namespace boost::multi

namespace boost::multi {

template<class Q, class R>
constexpr auto default_allocator_of(::thrust::pointer<Q, ::thrust::cuda_cub::tag, Q&> /*unused*/) {
	return ::thrust::cuda::universal_allocator<typename std::iterator_traits<::thrust::pointer<Q, ::thrust::cuda_cub::tag, Q&>>::value_type>{};
}

// copy_n
#if 1
template<class Q1, class L1, class Size, class Q2, class R2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<                  Q1*                                                , L1>   first, Size count,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if constexpr(L1::dimensionality == 1 and L2::dimensionality == 1) {
			if(first.layout().stride() == 1 and d_first.layout().stride() == 1) {
				auto s = cudaMemcpy  (raw_pointer_cast(d_first.current()),                                                                 first.current(),                                                               sizeof(Q2)* static_cast<std::size_t>(count), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
			} else {
				auto s = cudaMemcpy2D(raw_pointer_cast(d_first.current()), static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), first.current(), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), sizeof(Q2), static_cast<std::size_t>(count), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
			}
			return d_first + count;
		} else if constexpr(L1::dimensionality == 2 and L1::dimensionality == 2) {
			if(std::get<1>(first.layout().strides()) == 1 and std::get<1>(d_first.layout().strides()) == 1 and count%std::get<1>(first.layout().sizes()) == 0) {
				auto s = cudaMemcpy2D(raw_pointer_cast(d_first.current()), static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), first.current(), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), static_cast<std::size_t>(std::get<1>(first.layout().sizes()))*sizeof(Q2), static_cast<std::size_t>(count/std::get<1>(first.layout().sizes())), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
				return d_first + count;
			}  // else fallthrough
		}
		cudaHostRegister(
			const_cast<void*>(static_cast<void const*>(first.base())),
			static_cast<std::size_t>                  (first.layout().hull_size()*sizeof(Q1)),
			cudaHostRegisterPortable
		);
		auto ret = ::thrust::copy_n(
			::thrust::cuda::par,
			first, count, d_first
		);
		cudaHostUnregister(
			const_cast<void*>(static_cast<void const*>(first.base()))
		);
		return ret;
	} else {
		return ::thrust::copy_n(first, count, d_first);
	}
	return d_first + count;
}

template<class Q1, class R1, class L1, class Size, class Q2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1, ::thrust::use_default>, L1>   first, Size count,
	boost::multi::elements_iterator_t<                  Q2*                                                    , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                                                    , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if constexpr(L1::dimensionality == 1 and L2::dimensionality == 1) {
			if(first.layout().stride() == 1 and d_first.layout().stride() == 1) {
				auto s = cudaMemcpy  (                 d_first.current() ,                                                                 raw_pointer_cast(first.current()),                                                               sizeof(Q2)* static_cast<std::size_t>(count), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
			} else {
				auto s = cudaMemcpy2D(                 d_first.current() , static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), raw_pointer_cast(first.current()), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), sizeof(Q2), static_cast<std::size_t>(count), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
			}
			return d_first + count;
		} else if constexpr(L1::dimensionality == 2 and L1::dimensionality == 2) {
			if(std::get<1>(first.layout().strides()) == 1 and std::get<1>(d_first.layout().strides()) == 1 and count%std::get<1>(first.layout().sizes()) == 0) {
				auto s = cudaMemcpy2D(                 d_first.current() , static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), raw_pointer_cast(first.current()), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), static_cast<std::size_t>(std::get<1>(first.layout().sizes()))*sizeof(Q2), static_cast<std::size_t>(count/std::get<1>(first.layout().sizes())), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
				return d_first + count;
			}
		}
		cudaHostRegister(
			const_cast<void*>(static_cast<void const*>(d_first.base())),
			static_cast<std::size_t>                  (d_first.layout().hull_size()*sizeof(Q1)),
			cudaHostRegisterPortable
		);
		auto ret = ::thrust::copy_n(
			::thrust::cuda::par,
			first, count, d_first
		);
		cudaHostUnregister(
			const_cast<void*>(static_cast<void const*>(d_first.base()))
		);
		return ret;
	} else {
		return ::thrust::copy_n(first, count, d_first);
	}
	return d_first + count;
}

template<class Q1, class L1, class Size, class Q2, class R2, class L2>
auto uninitialized_copy_n(
	boost::multi::elements_iterator_t<                  Q1*                                                    , L1>   first, Size count,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2, ::thrust::use_default>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2, ::thrust::use_default>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy_n(first, count, d_first);
	} else {
		return ::thrust::uninitialized_copy_n(first, count, d_first);
	}
}

template<class Q1, class R1, class L1, class Size, class Q2, class L2>
auto uninitialized_copy_n(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1, ::thrust::use_default>, L1>   first, Size count,
	boost::multi::elements_iterator_t<                  Q2*                                                    , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                                                    , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy_n(first, count, d_first);
	} else {
		return ::thrust::uninitialized_copy_n(first, count, d_first);
	}
}

template<class Q1, class L1, class Q2, class R2, class L2>
auto copy(
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   first,
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   last ,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> {
	return boost::multi::copy_n(first, last - first, d_first);
}

template<class Q1, class R1, class L1, class Q2, class L2>
auto copy(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1>, L1>   first,
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1>, L1>   last ,
	boost::multi::elements_iterator_t<                  Q2*                             , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                             , L2> {
	return boost::multi::copy_n(first, last - first, d_first);
}

template<class Q1, class L1, class Q2, class R2, class L2>
auto uninitialized_copy(
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   first,
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   last ,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::cuda_cub::tag, R2>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy(first, last, d_first);
	} else {
		return ::thrust::uninitialized_copy(first, last, d_first);
	}
}

template<class Q1, class R1, class L1, class Q2, class L2>
auto uninitialized_copy(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1>, L1>   first,
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::cuda_cub::tag, R1>, L1>   last ,
	boost::multi::elements_iterator_t<                  Q2*                             , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                             , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy(first, last, d_first);
	} else {
		return ::thrust::uninitialized_copy(first, last, d_first);
	}
}

#endif

}  // end namespace boost::multi
