//#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
//$CXXX $CXXFLAGS $0 -o $0x -lcudart&&$0x&&rm $0x;exit
//#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_CACHED_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_CACHED_ALLOCATOR_HPP

#include "../../../adaptors/cuda/allocator.hpp"
#include "../../../adaptors/cuda/cached/ptr.hpp"

#include "../../../adaptors/cuda/cached/clib.hpp"  // cuda::malloc
#include "../../../adaptors/cuda/cached/malloc.hpp"

#include<cassert>
#include<cstddef>
#include<iostream>  // debug
#include<limits>
#include<new>  // bad_alloc
#include<list>
#include<unordered_map>

namespace boost::multi::memory::cuda {

namespace cached {

template <typename PointerType>
class allocator_cache {

	struct block {
		size_t size;
		PointerType loc;
	};

	static const size_t max_size = 4ul*1024ul*1024ul*1024ul;
	static const unsigned max_entries = 200;
	std::list<block> blocks_;
	std::unordered_multimap<size_t, typename decltype(blocks_)::iterator> map_;
	size_t mem_used;

 public:

	allocator_cache()
	: mem_used(0) {}

	auto put(size_t size, PointerType loc) {

		if(size >= max_size) return false;

		while(size + mem_used > max_size or map_.size() >= max_entries) {
			assert(map_.size() > 0);

			cuda::cached::free(blocks_.back().loc);
			mem_used -= blocks_.back().size;
			auto range = map_.equal_range(blocks_.back().size);
			for(auto it = range.first; it != range.second; ++it) {
				if(it->second == --blocks_.end()) {
					map_.erase(it);
					break;
				}
			};
			blocks_.pop_back();
		}

		blocks_.emplace_front(block{size, loc});
		map_.emplace(size, blocks_.begin());
		assert(map_.size() == blocks_.size());
		mem_used += size;

		return true;
	}

	PointerType get(size_t size) {
		PointerType loc;
		auto pos = map_.find(size);
		if(pos != map_.end()) {
			auto block_pos = pos->second;
			loc = block_pos->loc;
			blocks_.erase(block_pos);
			map_.erase(pos);
			mem_used -= size;
		} else {
			loc = nullptr;
		}
		assert(map_.size() == blocks_.size());
		return loc;
	}
};

auto & cache() {
	static allocator_cache<cached::ptr<void>> alloc_cache;
	return alloc_cache;
}

struct bad_alloc : std::bad_alloc {};

template<class T, class PrefetchDevice>
class allocator : cuda::allocator<T> {
	static_assert( std::is_same<T, std::decay_t<T>>{}, "!" );

 public:
	using value_type = T;
	using pointer = cached::ptr<T>;
	using size_type = ::size_t; // as specified by CudaMalloc
	using const_void_pointer = cached::ptr<void const>;
	template<class TT> using rebind = cached::allocator<TT, PrefetchDevice>;

	pointer allocate(typename allocator::size_type n) {
		MULTI_MARK_SCOPE("cuda::cached::allocate");

		if(n == 0) return pointer{nullptr};

		auto ret = static_cast<pointer>(cache().get(n*sizeof(T)));
		if(ret == pointer{nullptr}) {
			ret = static_cast<pointer>(cuda::cached::malloc(n*sizeof(T)));
			if(!ret) throw bad_alloc{};
		//  ++allocator::n_allocations; allocator::bytes_allocated+=sizeof(T)*n;
		}
		if(PrefetchDevice::value != -99) {
			auto const code = cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), PrefetchDevice::value);
			if(code != cudaSuccess) {
				throw std::runtime_error{"cannot prefetch for reason "+std::to_string(code)+" device is "+std::to_string(PrefetchDevice::value)};
			}
		}
		return ret;
	}

	pointer allocate(typename allocator::size_type n, const_void_pointer hint){
		auto const ret = allocate(n);
		if(not hint) {
			if(cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), /*device*/ 0) != cudaSuccess) {throw std::runtime_error{"cannot prefetch"};}
			return ret;
		}

		cudaPointerAttributes attr; if(cudaPointerGetAttributes(&attr, raw_pointer_cast(hint))!=cudaSuccess) {throw std::runtime_error{"cannot use attributes for hint"};}
		switch(attr.type) {
			case cudaMemoryTypeUnregistered: {//std::cout<< n <<" cudaMemoryTypeUnregistered"<< attr.device <<" "<< attr.device <<" cpuid:"<< cudaCpuDeviceId <<std::endl;
				auto code = cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId);
				if(code != cudaSuccess) {
					throw std::runtime_error{"1. cannot prefetch for reason "+std::to_string(code)+" device is "+std::to_string(cudaCpuDeviceId)};
				}
				return ret;
			}
			case cudaMemoryTypeHost        : {//std::cout<< n <<" cudaMemoryTypeHost "<< attr.device <<" "<< cudaCpuDeviceId <<std::endl;
				auto code = cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), cudaCpuDeviceId);
				if(code != cudaSuccess) {
					throw std::runtime_error{"2. cannot prefetch for reason "+std::to_string(code)+" device is "+std::to_string(cudaCpuDeviceId)};
				}
				return ret;
			}
			case cudaMemoryTypeDevice      : {//std::cout<< n <<" cudaMemoryTypeDevice "<< attributes.device <<" "<< attributes.device<<std::endl;
				auto code = cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device);
				if(code != cudaSuccess) {
					throw std::runtime_error{"3. cannot prefetch for reason "+std::to_string(code)+" device is "+std::to_string(attr.device)};
				}
				return ret;
			}
			case  cudaMemoryTypeManaged    : {//std::cout<< n <<" cudaMemoryTypeCached "<< attr.device <<" "<< attr.device <<std::endl;
				auto code = cudaMemPrefetchAsync(raw_pointer_cast(ret), n*sizeof(T), attr.device);
				if(code != cudaSuccess) {
					throw std::runtime_error{"4. cannot prefetch for reason "+std::to_string(code)+" device is "+std::to_string(attr.device)};
				}
				return ret;
			}
		}
		return ret;
	}

	void deallocate(pointer p, size_type n) {
		MULTI_MARK_SCOPE("cuda::cached::deallocate");

		if(not cache().put(n*sizeof(T), static_cast<cached::ptr<void>>(p))) {
			cuda::cached::free(static_cast<cached::ptr<void>>(p));
		}
	}

	template<class P, class... Args>
	void construct(P p, Args&&... args) {
		::new(p.rp_) T(std::forward<Args>(args)...);
	}
	
	template<class P, class... Args>
	void construct(P* p, Args&&... args) {
		::new(p) T(std::forward<Args>(args)...);
	}

	template<class P> void destroy(P  p) {p.rp_->~T();}
	template<class P> void destroy(P* p) {p->~T();}

	constexpr bool operator==(allocator<T> const&) const {return true;}
	constexpr bool operator!=(allocator<T> const&) const {return false;}

	template<class InputIt, class ForwardIt>
	constexpr ForwardIt alloc_uninitialized_copy(InputIt first, InputIt last, ForwardIt d_first) const {
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

}  // end namespace boost::multi::memory::cuda

//#if not __INCLUDE_LEVEL__

//#include<memory>
//#include<iostream>
//#include "../../../../array.hpp"

//namespace multi = boost::multi;
//namespace cuda = multi::memory::cuda;

//int main(){

//	multi::array<double, 1, multi::memory::cuda::cached::allocator<double> > A(32);
//	A[17] = 3.;
//	assert( A[17] == 3. );

//}
//#endif
#endif
