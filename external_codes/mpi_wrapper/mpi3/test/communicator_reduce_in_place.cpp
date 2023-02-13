// Copyright 2018-2021 Alfredo A. Correa

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {

	{
		int n = 1;
		world.reduce_in_place_n(&n, 1, std::plus<>{});
		if(world.rank() == 0) {assert(n == world.size());}
	}
	{
		int n = 1;
		world.all_reduce_in_place_n(&n, 1, std::plus<>{});
		assert(n == world.size());
	}
	{
		int n = 1;
		world.all_reduce_n(&n, 1, std::plus<>{});
		assert(n == world.size());
	}
	{
		int n = 1;
		world.all_reduce_n(&n, 1);
		assert(n == world.size());
	}
	{
		int n = 1;
		auto const m = world.all_reduce_value(n);
		assert( n == 1 );
		assert(m == world.size());
	}
	{
		int n = 1;
		auto const m = (world += n);
		assert( n == 1 );
		assert(m == world.size());
	}
//	{
//		int n = 1;
//		auto const m = (world + n);
//		assert( n == 1 );
//		assert(m == world.size());
//	}
	return 0;
} catch(...) {return 1;}
