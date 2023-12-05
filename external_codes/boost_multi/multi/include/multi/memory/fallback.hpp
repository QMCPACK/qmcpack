#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
#ifndef MULTI_MEMORY_FALLBACK_HPP
#define MULTI_MEMORY_FALLBACK_HPP

#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource>=201603L)
#include<memory_resource>
#endif

#include "../memory/block.hpp"
#include "../memory/allocator.hpp"

#include<stdlib.h> // aligned_alloc, in c++17 this will be <cstdlib>
#include <cstddef> // std::max_align_t

namespace boost {
namespace multi {
namespace memory {

template<class Ptr = memory::byte*> struct resource;

template<class T>
struct resource<T*>{
	auto allocate(std::size_t size, std::size_t alignment = alignof(std::max_align_t)){
		return ::aligned_alloc(alignment, size);
	}
	void deallocate(void* ptr, std::size_t = alignof(std::max_align_t)){::free(ptr);}
};

template<class Ptr = memory::byte*>
inline auto get_default_resource(){
	static resource<Ptr> instance_;
	return &instance_;
}

template<class MemoryResource1, class MemoryResource2
#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource>=201603L)
	= std::pmr::memory_resource
#else
	= memory::resource<>
#endif
>
class fallback : public MemoryResource1{
	MemoryResource2* back_ = nullptr;
	long fallbacks_ = 0;
public:
	long fallbacks() const{return fallbacks_;}
	fallback() = default;

	fallback(MemoryResource1 const& mr, MemoryResource2* back) : MemoryResource1{mr}, back_{back} {}

	// cppcheck-suppress noExplicitConstructor ; allocators are pointers to memory resources
	fallback(MemoryResource1 const& mr)
	: fallback{
		mr,
#if(__cpp_lib_memory_resource>=201603L)
		std::pmr::get_default_resource()
#else
		nullptr//memory::get_default_resource<>()
#endif
	} {}


	typename fallback::void_pointer 
	allocate(size_type required_bytes, typename fallback::size_type align = alignof(std::max_align_t)) try{
		return MemoryResource1::allocate(required_bytes, align);
	}catch(...){
		++fallbacks_;
		return back_->allocate(required_bytes, align);
	}
	void deallocate(typename fallback::void_pointer p, typename fallback::size_type discarded_bytes) try{
		MemoryResource1::deallocate(p, discarded_bytes);
	}catch(...){back_->deallocate(p, discarded_bytes);}
};

template<class T, class MR1, class MR2
#if(__cpp_lib_memory_resource>=201603L)
	= std::pmr::memory_resource //= std::allocator<std::byte>
#else
	= memory::resource<>
#endif
> 
using fallback_allocator = memory::allocator<T, fallback<MR1, MR2>>;//, alignof(T)>>;

}  // end namespace memory
}  // end namespace multi
}  // end namespace boost

#if not __INCLUDE_LEVEL__ //  _TEST_MULTI_MEMORY_FALLBACK

#include "../../multi/array.hpp"
#include "../memory/stack.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
namespace memory = multi::memory;
using std::cout;

int main() {
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	memory::monotonic<char*> m(buffer.data(), buffer.size());
	auto p1 = m.allocate(100*sizeof(double), alignof(double));
	try {
		auto p2 = m.allocate(200*sizeof(double), alignof(double));
		m.deallocate(p2, 200*sizeof(double));
	} catch(std::bad_alloc& e) {std::cout <<"throw "<< e.what() << std::endl;}
	m.deallocate(p1, 100*sizeof(double));
}
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	memory::fallback<memory::monotonic<char*>> m(buffer.data(), buffer.size());//, boost::multi::memory::get_default_resource());
	auto p1 = m.allocate(100*sizeof(double), alignof(double));
	auto p2 = m.allocate(200*sizeof(double), alignof(double));
	m.deallocate(p2, 200*sizeof(double));
	m.deallocate(p1, 100*sizeof(double));
	assert( m.fallbacks() == 1 );
}
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	memory::stack<char*> s(buffer.data(), buffer.size());
	auto p1 = s.allocate(1*sizeof(double), alignof(double));
	auto p2 = s.allocate(100*sizeof(double), alignof(double));
	s.deallocate(p2, 100*sizeof(double));
	s.deallocate(p1, 1*sizeof(double));
	assert( s.max_needed() == 101*sizeof(double) );
}
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	memory::fallback<memory::stack<char*>> s(buffer.data(), buffer.size());
	auto p1 = s.allocate(10000*sizeof(double), alignof(double));
	s.deallocate(p1, 10000*sizeof(double));
	assert( s.fallbacks() == 1 );
}
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	memory::fallback<memory::stack<char*>> s(buffer.data(), buffer.size());
	auto p2 = s.allocate(255*sizeof(double), alignof(double));
	auto p3 = s.allocate(255*sizeof(double), alignof(double));
	s.deallocate(p3, 255*sizeof(double));
	s.deallocate(p2, 255*sizeof(double));
	assert( s.fallbacks() == 1 );
	assert( s.max_needed() == (255+255)*sizeof(double) );
}
{
	std::size_t guess = 1000;
	for(int i = 0; i != 3; ++i) {
		std::vector<char> buffer(guess);
		memory::fallback<memory::stack<char*>> f({buffer.data(), (std::ptrdiff_t)buffer.size()});
		std::vector<double, memory::fallback_allocator<double, memory::stack<char*>>> v(10, &f);
		std::vector<int, memory::fallback_allocator<int, memory::stack<char*>>> w(1000, &f);
		cout<<"iteration: "<< i <<"-> fallbacks "<< f.fallbacks() <<" max_needed " << f.max_needed() <<std::endl;
		guess = std::max(guess, (std::size_t)f.max_needed());
	}
}
cout <<"----------"<< std::endl;
{
	memory::fallback<memory::stack<char*>> f;//({buffer.data(), (std::ptrdiff_t)buffer.size()}); // TODO error in nvcc
	for(int i = 0; i != 3; ++i) {
		std::vector<char> buffer(f.max_needed());
		f = {{buffer.data(), (std::ptrdiff_t)buffer.size()}};
		std::vector<double, memory::fallback_allocator<double, memory::stack<char*>>> v(10, &f);
		std::vector<int, memory::fallback_allocator<int, memory::stack<char*>>> w(1000, &f);
		cout<<"iteration: "<< i <<"-> fallbacks "<< f.fallbacks() <<" max_needed " << f.max_needed() <<std::endl;
	}
}
cout <<"----------"<< std::endl;
{
	memory::fallback<memory::stack<memory::byte*>> f;//({buffer.data(), (std::ptrdiff_t)buffer.size()});
	for(int i = 0; i != 3; ++i) {
		std::vector<memory::byte> buffer(f.max_needed());
		f = {{buffer.data(), (std::ptrdiff_t)buffer.size()}};
		std::vector<double, memory::fallback_allocator<double, memory::stack<memory::byte*>>> v(10, &f);
		std::vector<int, memory::fallback_allocator<int, memory::stack<memory::byte*>>> w(1000, &f);
		cout<<"iteration: "<< i <<"-> fallbacks "<< f.fallbacks() <<" max_needed " << f.max_needed() <<std::endl;
	}
}
cout <<"----------"<< std::endl;
{
	memory::fallback<memory::stack<>> f;
	for(int i = 0; i != 3; ++i) {
		std::vector<memory::byte> buffer(f.max_needed());
		f = {{buffer.data(), f.max_needed()}};
		std::vector<double, memory::fallback_allocator<double, memory::stack<>>> v(10, &f);
		std::vector<int, memory::fallback_allocator<int, memory::stack<>>> w(1000, &f);
		cout<<"iteration: "<< i <<"-> fallbacks "<< f.fallbacks() <<" max_needed " << f.max_needed() <<std::endl;
	}
}
}
#endif
#endif
