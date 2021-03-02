#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef BOOST_MULTI_MEMORY_ALLOCATOR_HPP
#define BOOST_MULTI_MEMORY_ALLOCATOR_HPP

#include "../detail/memory.hpp"
#include "../config/NODISCARD.hpp"

#include<cassert>

#if defined(__cpp_lib_memory_resource) and __cpp_lib_memory_resource>=201603L
#include<memory_resource>
#endif

namespace boost{
namespace multi{
namespace memory{

template<class T, class Memory
#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource>=201603L)
	= std::pmr::memory_resource //= std::allocator<std::byte>
#endif
	, typename Constructor = std::allocator<T>
>
class allocator{
	using memory_type = Memory;
	memory_type* mp_ 
#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource>=201603L)
		= std::pmr::get_default_resource()
#endif
	;
	using constructor_type = typename Constructor::template rebind<T>::other;
	constructor_type ctor_;
public:
	using value_type = T;
	using pointer = typename std::pointer_traits<decltype(std::declval<memory_type*>()->allocate(0, 0))>::template rebind<value_type>;
	using difference_type = typename std::pointer_traits<pointer>::difference_type;
	using size_type =  std::make_unsigned_t<difference_type>;
//	static_assert( std::is_same<pointer, typename constructor_type::pointer>{}, 
//		"incompatible pointer types");
	allocator() = default;
	allocator(memory_type* mp, constructor_type const& ctor = {}) : mp_{mp}, ctor_{ctor}{}
	allocator(constructor_type const& ctor) : mp_{}, ctor_{ctor}{}
	allocator(allocator const& o) = default;//: mp_{o.mp_}, ctor_{o.ctor_}{}
	allocator& operator=(allocator const&) = default;
	bool operator==(allocator const& o) const{return mp_ == o.mp_;}
	bool operator!=(allocator const& o) const{return mp_ != o.mp_;}
	using void_pointer = typename std::pointer_traits<decltype(std::declval<memory_type*>()->allocate(0, 0))>::template rebind<void>;
	NODISCARD("because otherwise it will generate a memory leak")
	pointer allocate(size_type n){
		static_assert( alignof(std::max_align_t) >= 16 , "for cuda");
		return static_cast<pointer>(static_cast<void_pointer>(mp_->allocate(
			n*sizeof(value_type),
			alignof(std::max_align_t) // this takes into account 
		)));
	}
	void deallocate(pointer p, size_type n){
		mp_->deallocate(p, n*sizeof(value_type));
	}
	template<class... Args>
	void construct(pointer p, Args&&... args){
		allocator_traits<Constructor>::construct(ctor_, p, std::forward<Args>(args)...);
	}
	decltype(auto) destroy(pointer p){
		allocator_traits<Constructor>::destroy(ctor_, p);
	}
};

}}}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MULTI_MEMORY_ALLOCATOR
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi memory allocator"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../memory/monotonic.hpp"
#include<boost/align/is_aligned.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_memory_allocator){

	alignas(double) char buffer[280*sizeof(double)];
	multi::memory::monotonic<char*> m(buffer);
	multi::memory::allocator<double, multi::memory::monotonic<char*> > A(&m);
	double* p = A.allocate(1);
	A.construct(p, 8.);
	BOOST_REQUIRE( *p == 8. );
	BOOST_REQUIRE( boost::alignment::is_aligned(p, alignof(double)) );

	double* arr = A.allocate(255);
	A.construct(arr, 81.);
	BOOST_REQUIRE( *arr == 81. );
}
#endif
#endif

