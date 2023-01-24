// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>

#include<multi/array.hpp>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {

	{
		auto v = std::vector<double>(10, world.root()? 3.14 : 99.9);

		auto it = world.broadcast_n(v.begin(), 10);

		assert( it == v.end() );
		assert( std::find_if(v.cbegin(), v.cend(), [](auto& e) {return e != 3.14;}) == v.cend() );
	}
	{
		namespace multi = boost::multi;
		auto arr = multi::array<double, 1>(10, world.root()? 3.14 : 99.9);

		auto it = world.broadcast_n(arr.begin(), 10);

		assert( it == arr.end() );
		assert( std::find_if(arr.cbegin(), arr.cend(), [](auto& e) {return e != 3.14;}) == arr.cend() );
	}

	return 0;
} catch(...) {return 1;}
