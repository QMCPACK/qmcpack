#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0.$X&&$0.$X&&rm $0.$X;exit
#endif
//  © Alfredo A. Correa 2020

#include "../../multi/array.hpp"

#include <iostream>
#include <memory_resource>  // polymorphic memory resource, monotonic buffer

namespace multi = boost::multi;

int main() {
	char buffer[13] = "____________"; // a small buffer on the stack or an allocation
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> A({2, 2}, 'a', &pool);
	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> B({3, 2}, 'b', &pool);
	std::cout << buffer << std::endl;
}

