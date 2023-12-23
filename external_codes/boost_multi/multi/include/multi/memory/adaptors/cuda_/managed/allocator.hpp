#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lcudart&&$0x&&rm $0x;exit
#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP

#include "../../../adaptors/cuda/allocator.hpp"
#include "../../../adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/cuda/managed/clib.hpp" // cuda::malloc
#include "../../../adaptors/cuda/managed/malloc.hpp"

#include<cassert>
#include<cstddef>
#include<iostream> // debug
#include<limits>
#include<new> // bad_alloc

namespace boost{namespace multi{
namespace memory{namespace cuda{

namespace managed{
	struct bad_alloc : std::bad_alloc{};

	template<class T, class PrefetchDevice>
	class allocator : cuda::allocator<T>{
		static_assert( std::is_same<T, std::decay_t<T>>{}, "!" );
	public:
		using value_type = T;
		using pointer = managed::ptr<T>;
		using size_type = ::size_t; // as specified by CudaMalloc
		using const_void_pointer = managed::ptr<void const>;
		template<class TT> using rebind = managed::allocator<TT, PrefetchDevice>;
		pointer allocate(typename allocator::size_type n){
			if(n == 0) return pointer{nullptr};
			auto ret = static_cast<pointer>(cuda::managed::malloc(n*sizeof(T)));
			if(!ret) throw bad_alloc{};
			if(PrefetchDevice::value != -99)
				if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), PrefetchDevice::value) != cudaSuccess) throw std::runtime_error{"cannot prefetch for some reason"};
	//		++allocator::n_allocations; allocator::bytes_allocated+=sizeof(T)*n;
			return ret;
		}
		pointer allocate(typename allocator::size_type n, const_void_pointer hint){
			MULTI_MARK_SCOPE("cuda::managed::allocate");
			
			auto const ret = allocate(n);
			if(not hint){
				if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), /*device*/ 0) != cudaSuccess) throw std::runtime_error{"cannot prefetch"};
				return ret;
			}
			cudaPointerAttributes attr; if(cudaPointerGetAttributes(&attr, raw_pointer_cast(hint))!=cudaSuccess) throw std::runtime_error{"cannot use attributes for hint"};
			switch(attr.type){
				case cudaMemoryTypeUnregistered:{//std::cout<< n <<" cudaMemoryTypeUnregistered"<< attr.device <<" "<< attr.device <<" cpuid:"<< cudaCpuDeviceId <<std::endl;
					if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
					return ret;
				}
				case cudaMemoryTypeHost        :{//std::cout<< n <<" cudaMemoryTypeHost "<< attr.device <<" "<< cudaCpuDeviceId <<std::endl;
					if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
					return ret;
				}
				case  cudaMemoryTypeDevice     :{//std::cout<< n <<" cudaMemoryTypeDevice "<< attributes.device <<" "<< attributes.device<<std::endl;
					if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
					return ret;
				}
				case  cudaMemoryTypeManaged    :{//std::cout<< n <<" cudaMemoryTypeManaged "<< attr.device <<" "<< attr.device <<std::endl;
					if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device /*0?*/) != cudaSuccess) throw std::runtime_error{"could not prefetch in managed memory"};
					return ret;
				}
			}
			return ret;
		}
		void deallocate(pointer p, size_type){
			MULTI_MARK_SCOPE("cuda::managed::deallocate");
			cuda::managed::free(static_cast<managed::ptr<void>>(p));
		}
		template<class P, class... Args>
		void construct(P p, Args&&... args){
			::new(p.rp_) T(std::forward<Args>(args)...);
		}
		template<class P, class... Args>
		void construct(P* p, Args&&... args){
			::new(p) T(std::forward<Args>(args)...);
		}
		template<class P> void destroy(P p){p.rp_->~T();}
		template<class P> void destroy(P* p){p->~T();}
		constexpr bool operator==(allocator<T> const&) const{return true;}
		constexpr bool operator!=(allocator<T> const&) const{return false;}
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
}

}}}}

#if not __INCLUDE_LEVEL__

#include<memory>
#include<iostream>
#include "../../../../array.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

int main(){

	multi::array<double, 1, multi::memory::cuda::managed::allocator<double> > A(32);
	A[17] = 3.;
	assert( A[17] == 3. );

}
#endif
#endif

