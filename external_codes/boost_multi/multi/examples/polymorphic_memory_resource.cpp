#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
${CXX:-c++} -std=c++17 $CXXFLAGS $0 -o $0x&&$0x&&rm $0x;exit
#endif
//  © Alfredo A. Correa 2020

#include "../../multi/array.hpp"

#include <iostream>
#include <memory_resource>  // polymorphic memory resource, monotonic buffer

namespace multi = boost::multi;

int main() {
	static_assert( sizeof(multi::array<char, 2>) < sizeof(multi::pmr::array<char, 2>) , "!");

	char buffer[13] = "____________";  // flawfinder: ignore , a small buffer on the stack or an allocation
	std::pmr::monotonic_buffer_resource pool{
		std::data(buffer), std::size(buffer), 
		std::pmr::null_memory_resource()
	};

	multi::pmr::array<char, 2> A({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> B({3, 2}, 'b', &pool);

	assert( A.get_allocator() == B.get_allocator() );
	assert( buffer == std::string{"aaaabbbbbb__"} );

	try {
		multi::pmr::array<char, 2> C({9, 9}, 'c', &pool); // there is no upstream resource so it throws
	} catch(std::bad_alloc&) {
		assert( buffer == std::string{"aaaabbbbbb__"} );
	}

	std::array<char, 99> buffer2;
	std::pmr::monotonic_buffer_resource pool2{
		buffer.data(), buffer.size(),
		std::pmr::null_memory_resource()
	};
	{
		multi::pmr::array<char, 2> D = A;
		assert(D.get_allocator() != A.get_allocator() );
		assert(D.get_allocator().resource() == std::pmr::get_default_resource() );
	}
	{
		multi::pmr::array<char, 2> D({2, 2}, &pool2);
		D = A;
		assert(D.get_allocator() != A.get_allocator() );
	}
	{
		multi::pmr::array<char, 2> D(&pool2);
		D = std::move(A);
		assert(D.get_allocator() != A.get_allocator() );
		assert(D[1][1] == 'a');
		assert(A.is_empty());
	}
	{
		multi::pmr::array<char, 2> D(&pool2);
		std::swap(D, B);
		assert(( D.get_allocator() == multi::pmr::array<char, 2>::allocator_type{&pool2} ));
		assert( B.is_empty() );
	}
}
