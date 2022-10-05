// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2022

#ifndef MULTI_MEMORY_ALLOCATOR_HPP
#define MULTI_MEMORY_ALLOCATOR_HPP

#include "../config/NODISCARD.hpp"
#include "../detail/memory.hpp"

#include<cassert>

#include<memory_resource>

namespace boost::multi::memory {

template<class T, class Memory
	= std::pmr::memory_resource
	, typename Constructor = std::allocator<T>
>
class allocator {
	using memory_type = Memory;
	memory_type* mp_ = std::pmr::get_default_resource();
	using constructor_type = typename std::allocator_traits<Constructor>::template rebind_alloc<T>;
	constructor_type ctor_;

 public:
	using value_type = T;
	using pointer = typename std::pointer_traits<decltype(std::declval<memory_type*>()->allocate(0, 0))>::template rebind<value_type>;
	using difference_type = typename std::pointer_traits<pointer>::difference_type;
	using size_type =  std::make_unsigned_t<difference_type>;

	allocator(allocator const& o) = default;
	allocator(allocator&&) noexcept = default;

	auto operator=(allocator const&) -> allocator& = default;
	auto operator=(allocator&&) noexcept -> allocator& = default;

	allocator() = default;
	~allocator() = default;

	allocator(memory_type* mp, constructor_type const& ctor) : mp_{mp}, ctor_{ctor} {}
	// cppcheck-suppress noExplicitConstructor ; allocator *is* a pointer to a heap
	allocator(memory_type* mp) : allocator(mp, constructor_type{}) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : to allow pointer syntax

	explicit allocator(constructor_type const& ctor) : mp_{}, ctor_{ctor} {}

	auto operator==(allocator const& o) const -> bool {return mp_ == o.mp_;}
	auto operator!=(allocator const& o) const -> bool {return mp_ != o.mp_;}

	using void_pointer = typename std::pointer_traits<decltype(std::declval<memory_type*>()->allocate(0, 0))>::template rebind<void>;

	NODISCARD("because otherwise it will generate a memory leak")
	auto allocate(size_type n) -> pointer {
		return static_cast<pointer>(static_cast<void_pointer>(mp_->allocate(
			n*sizeof(value_type),
			std::max(
				alignof(std::max_align_t),
				std::size_t{16}  // for cuda
			)
		)));
	}
	void deallocate(pointer p, size_type n) {
		mp_->deallocate(p, n*sizeof(value_type));
	}
	template<class... Args>
	void construct(pointer p, Args&&... args) {
		allocator_traits<Constructor>::construct(ctor_, p, std::forward<Args>(args)...);
	}
	void destroy(pointer p) {allocator_traits<Constructor>::destroy(ctor_, p);}
};

}  // end namespace boost::multi::memory

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi memory allocator"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../memory/monotonic.hpp"
//#include<boost/align/is_aligned.hpp>

//namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(multi_memory_allocator) {
//	alignas(double) char buffer[280*sizeof(double)];  // flawfinder: ignore, test for legacy type
//	multi::memory::monotonic<char*> m(buffer);
//	multi::memory::allocator<double, multi::memory::monotonic<char*> > A(&m);
//	double* p = A.allocate(1);
//	A.construct(p, 8.);
//	BOOST_REQUIRE( *p == 8. );
//	BOOST_REQUIRE( boost::alignment::is_aligned(p, alignof(double)) );

//	double* arr = A.allocate(255);
//	A.construct(arr, 81.);
//	BOOST_REQUIRE( *arr == 81. );
//}

//BOOST_AUTO_TEST_CASE(multi_memory_allocator) {
//	alignas(double) std::array<char, 280*sizeof(double)> buffer;  // char buffer[280*sizeof(double)];
//	multi::memory::monotonic<char*> m(buffer.data(), buffer.size());
//	multi::memory::allocator<double, multi::memory::monotonic<char*> > A(&m);
//	double* p = A.allocate(1);
//	A.construct(p, 8.);
//	BOOST_REQUIRE( *p == 8. );
//	BOOST_REQUIRE( boost::alignment::is_aligned(p, alignof(double)) );

//	double* arr = A.allocate(255);
//	A.construct(arr, 81.);
//	BOOST_REQUIRE( *arr == 81. );
//}
//#endif
#endif
