#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lcudart&&$0x&&rm $0x;exit
#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP

#include "../../../adaptors/cuda/allocator.hpp"
#include "../../../adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/cuda/managed/clib.hpp" // cuda::malloc
#include "../../../adaptors/cuda/managed/malloc.hpp"

#include<new> // bad_alloc
#include<cassert>
#include<iostream> // debug

#include<cstddef>

namespace boost{namespace multi{
namespace memory{namespace cuda{

namespace managed{
	struct bad_alloc : std::bad_alloc{};

	template<class T=void>
	class allocator : cuda::allocator<T>{
		static_assert( std::is_same<T, std::decay_t<T>>{}, "!" );
	public:
		using value_type = T;
		using pointer = managed::ptr<T>;
		using size_type = ::size_t; // as specified by CudaMalloc
		pointer allocate(typename allocator::size_type n, const void* = 0){
			if(n == 0) return pointer{nullptr};
			auto ret = static_cast<pointer>(cuda::managed::malloc(n*sizeof(T)));
			if(!ret) throw bad_alloc{};
			++allocator::n_allocations; allocator::bytes_allocated+=sizeof(T)*n;
			return ret;
		}
		void deallocate(pointer p, size_type){cuda::free(static_cast<managed::ptr<void>>(p));}
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

