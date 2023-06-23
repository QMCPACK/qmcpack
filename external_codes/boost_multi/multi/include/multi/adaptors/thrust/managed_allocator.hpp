#ifndef MULTI_MEMORY_ADAPTORS_THRUST_MANAGED_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_THRUST_MANAGED_ALLOCATOR_HPP
// Copyright 2018-2022 Alfredo A. Correa

//#include "../../../adaptors/cuda/allocator.hpp"
//#include "../../../adaptors/cuda/managed/ptr.hpp"

//#include "../../../adaptors/cuda/managed/clib.hpp" // cuda::malloc
//#include "../../../adaptors/cuda/managed/malloc.hpp"

#include<thrust/system/cuda/pointer.h>

#include<driver_types.h>     // cudaError_t
#include<cuda_runtime_api.h> // cudaGetErrorString

//#include<cassert>
//#include<cstddef>
//#include<iostream> // debug
//#include<limits>
//#include<new> // bad_alloc

namespace boost {
namespace multi {

namespace thrust {
namespace cuda {

struct bad_alloc : std::bad_alloc{};

template<class T>
class managed_allocator {
	static_assert( std::is_same<T, std::decay_t<T>>{}, "!" );

 public:
//	using allocator_type = managed_allocator<T>;
	using value_type = T;
	using pointer = ::thrust::cuda::pointer<T>;
	using size_type = ::size_t; // as specified by CudaMalloc
	using const_void_pointer = ::thrust::cuda::pointer<void const>;

	template<class TT> using rebind = managed_allocator<TT>;  //, PrefetchDevice>;

 private:
	static inline cudaError_t Malloc(void** p, size_type bytes) {return cudaMalloc(p, bytes);}
	static inline void* malloc(size_type bytes) {
		void* ret;
		switch(auto e = Malloc(&ret, bytes)){
			case cudaSuccess               : break;
			case cudaErrorMemoryAllocation : ret = nullptr; break;
			default                        : throw std::runtime_error("cannot allocate "+std::to_string(bytes)+" bytes in '"+__PRETTY_FUNCTION__+"' because "+cudaGetErrorString(e));
		}
		return ret;
	}

 public:
	pointer allocate(size_type n) {
		if(n == 0) return pointer{nullptr};
		auto ret = static_cast<T*>(malloc(n*sizeof(T)));
		if(!ret) throw bad_alloc{};
		return pointer{ret};
	}
	pointer allocate(size_type n, const_void_pointer hint) {
		MULTI_MARK_SCOPE("thrust::managed_allocate");

		auto const ret = allocate(n);
		// if(not hint){
		// 	if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), /*device*/ 0) != cudaSuccess) throw std::runtime_error{"cannot prefetch"};
		// 	return ret;
		// }
		// cudaPointerAttributes attr; if(cudaPointerGetAttributes(&attr, raw_pointer_cast(hint))!=cudaSuccess) throw std::runtime_error{"cannot use attributes for hint"};
		// switch(attr.type){
		// 	case cudaMemoryTypeUnregistered:{//std::cout<< n <<" cudaMemoryTypeUnregistered"<< attr.device <<" "<< attr.device <<" cpuid:"<< cudaCpuDeviceId <<std::endl;
		// 		if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
		// 		return ret;
		// 	}
		// 	case cudaMemoryTypeHost        :{//std::cout<< n <<" cudaMemoryTypeHost "<< attr.device <<" "<< cudaCpuDeviceId <<std::endl;
		// 		if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
		// 		return ret;
		// 	}
		// 	case  cudaMemoryTypeDevice     :{//std::cout<< n <<" cudaMemoryTypeDevice "<< attributes.device <<" "<< attributes.device<<std::endl;
		// 		if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
		// 		return ret;
		// 	}
		// 	case  cudaMemoryTypeManaged    :{//std::cout<< n <<" cudaMemoryTypeManaged "<< attr.device <<" "<< attr.device <<std::endl;
		// 		if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device /*0?*/) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
		// 		return ret;
		// 	}
		// }
		return ret;
	}
	void deallocate(pointer p, size_type) {
		MULTI_MARK_SCOPE("thrust::managed_deallocate");
	//	cuda::managed::free(static_cast<managed::ptr<void>>(p));
		cudaFree(raw_pointer_cast(p));
	}
	template<class P, class... Args>
	void construct(P p, Args&&... args){
		::new(raw_pointer_cast(p)) T(std::forward<Args>(args)...);
	}
	template<class P, class... Args>
	void construct(P* p, Args&&... args){
		::new(raw_pointer_cast(p)) T(std::forward<Args>(args)...);
	}
	template<class P> void destroy(P p){raw_pointer_cast(p)->~T();}
	template<class P> void destroy(P* p){p->~T();}

	constexpr bool operator==(managed_allocator const&) const {return true ;}
	constexpr bool operator!=(managed_allocator const&) const {return false;}

	template<class InputIt, class ForwardIt>
	constexpr ForwardIt alloc_uninitialized_copy(InputIt first, InputIt last, ForwardIt d_first) const{
		return ForwardIt{adl_uninitialized_copy(first, last, d_first)};
	}
	template<class InputIt, class Size, class ForwardIt>
	constexpr ForwardIt alloc_uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first) const{
		return ForwardIt{adl_uninitialized_copy_n(first, count, d_first)};
	}
	template<class ForwardIt, class Size>
	constexpr ForwardIt alloc_uninitialized_default_construct_n(ForwardIt first, Size n) const{
		return ForwardIt{adl_uninitialized_default_construct_n(first, n)};
	}
	template<class ForwardIt, class Size>
	constexpr ForwardIt alloc_destroy_n(ForwardIt first, Size n) const{return ForwardIt{destroy_n(first, n)};}
};

}  // end namespace cuda
}  // end namespace thrust
}  // end namespace multi
}  // end namespace boost

#endif
