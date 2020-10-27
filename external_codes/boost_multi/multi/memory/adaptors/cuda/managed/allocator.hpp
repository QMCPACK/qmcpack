#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
nvcc -D_TEST_MULTI_MEMORY_CUDA_MANAGED_ALLOCATOR -x cu $0 -o $0x&&$0x&&rm $0x;exit
#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_ALLOCATOR_HPP

#include<cuda_runtime.h> // cudaMalloc

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
		template<class InputIt, class Size, class ForwardIt>
		ForwardIt alloc_uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first) const{
			return adl_uninitialized_copy_n(first, count, d_first);
		}
		template<class ForwardIt, class Size>
		ForwardIt alloc_uninitialized_default_construct_n(ForwardIt first, Size n) const{
			return adl_uninitialized_default_construct_n(first, n);
		}
		template<class ForwardIt, class Size>
		ForwardIt alloc_destroy_n(ForwardIt first, Size n) const{return destroy_n(first, n);}
	};
}

}}}}

#ifdef _TEST_MULTI_MEMORY_CUDA_MANAGED_ALLOCATOR

#include<memory>
#include<iostream>
#include "../../../../array.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

int main(){

	multi::array<double, 1, multi::memory::cuda::managed::allocator<double> > A(32);
	A[17] = 3.;
	assert( A[17] == 3. );
#if 0
	{
		multi::static_array<double, 1> A(32, double{}); A[17] = 3.;
		multi::static_array<double, 1, cuda::allocator<double>> A_gpu = A;
		assert( A_gpu[17] == 3 );

		multi::static_array<double, 1, cuda::managed::allocator<double>> A_mgpu = A;
		assert( A_mgpu[17] == 3 );
		
		multi::static_array<double, 1, cuda::managed::allocator<double>> AA_mgpu = A_gpu;
		assert( A_mgpu[17] == 3 );
	}
	{
		multi::array<double, 1> A(32, double{}); A[17] = 3.;
		multi::array<double, 1, cuda::allocator<double>> A_gpu = A;
		assert( A_gpu[17] == 3 );

		multi::array<double, 1, cuda::managed::allocator<double>> A_mgpu = A;
		assert( A_mgpu[17] == 3 );
		
		multi::array<double, 1, cuda::managed::allocator<double>> AA_mgpu = A_gpu;
		assert( A_mgpu[17] == 3 );
	}
	{
		multi::array<double, 2> A_cpu({32, 64}, double{}); A_cpu[17][22] = 3.;
		multi::array<double, 2, cuda::allocator<double>> A_gpu = A_cpu;
		assert( A_gpu[17][22] == 3 );

		multi::array<double, 2, cuda::managed::allocator<double>> A_mgpu = A_cpu;
		assert( A_mgpu[17][22] == 3 );
		
		multi::array<double, 2, cuda::managed::allocator<double>> AA_mgpu = A_gpu;
		assert( AA_mgpu[17][22] == 3 );
	}
	{
		multi::static_array<double, 1> A1(32, double{}); A1[17] = 3.;
		multi::static_array<double, 1, cuda::managed::allocator<double>> A1_gpu = A1;
		assert( A1_gpu[17] == 3 );
	}

	{
		multi::array<double, 1> A1(32, double{}); A1[17] = 3.;
		multi::array<double, 1, cuda::allocator<double>> A1_gpu = A1;
		assert( A1_gpu[17] == 3 );
	}
	{
		multi::static_array<double, 2> A2({32, 64}, double{}); A2[2][4] = 8.;
		multi::static_array<double, 2, cuda::allocator<double>> A2_gpu = A2;
		assert( A2_gpu[2][4] == 8. );
	}
	{
		multi::array<double, 2> A2({32, 64}, double{}); A2[2][4] = 8.;
		multi::static_array<double, 2, cuda::allocator<double>> A2_gpu = A2;
		assert( A2_gpu[2][4] == 8. );
	}
	{
		multi::array<double, 2> A2({32, 64}, double{}); A2[2][4] = 8.;
		multi::array<double, 2, cuda::allocator<double>> A2_gpu = A2;
		assert( A2_gpu[2][4] == 8. );
	}
	{
		multi::array<double, 2> A2({32, 8000000}, double{}); A2[2][4] = 8.;
		multi::array<double, 2, cuda::allocator<double>> A2_gpu = A2;
		int s; std::cin >> s;
		assert( A2_gpu[2][4] == 8. );
	}
	{
		static_assert(std::is_same<std::allocator_traits<cuda::allocator<double>>::difference_type, std::ptrdiff_t>{}, "!");
		static_assert(std::is_same<std::allocator_traits<cuda::allocator<double>>::pointer, cuda::ptr<double>>{}, "!");
		static_assert(
			std::is_same<
				std::allocator_traits<cuda::allocator<int>>::rebind_alloc<double>,
				cuda::allocator<double>
			>{}, "!"
		);
		cuda::allocator<double> calloc;
		assert(calloc == calloc);
		cuda::ptr<double> p = calloc.allocate(100);
		p[33] = 123.;
		p[99] = 321.;
		p[33]+=1;
		double p33 = p[33];
		assert( p33 == 124. );
		assert( p[33] == 124. );
		assert( p[33] == p[33] );
		swap(p[33], p[99]);
		assert( p[99] == 124. );
		assert( p[33] == 321. );
		std::cout << p[33] << std::endl;
		calloc.deallocate(p, 100);
		p = nullptr;
		cout<<"n_alloc/dealloc "<< cuda::allocation_counter::n_allocations <<"/"<< cuda::allocation_counter::n_deallocations <<"\n"
			<<"bytes_alloc/dealloc "<< cuda::allocation_counter::bytes_allocated <<"/"<< cuda::allocation_counter::bytes_deallocated <<"\n";
	}
#endif
}
#endif
#endif

