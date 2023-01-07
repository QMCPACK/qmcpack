// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"

namespace bmpi3 = boost::mpi3;

auto bmpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try{
	{
		using T = int;
		std::vector<T> v(10);
		if(world.is_root()) {iota(begin(v), end(v), 0);}

		world.broadcast(begin(v), end(v));
		assert( v[9] == T(9) );
	}
	{
		using T = double;
		std::vector<T> v(10);
		if(world.is_root()) {
			iota(begin(v), end(v), 0);
		}
		world.broadcast(begin(v), end(v));

		assert( v[9] == T(9) );
	}

	{
		using T = double;
		std::vector<T> v;
		if(world.is_root()) {v.resize(10); iota(begin(v), end(v), 0);}

		auto size = v.size();
		world.broadcast_n(&size, 1);

		v.resize(size);
		assert( v.size() == 10UL );

		world.broadcast_n(v.data(), v.size());

		assert( v[9] == T(9) );
	}
	{
		using T = double;
		std::vector<T> v;
		if(world.is_root()) {v.resize(10); iota(begin(v), end(v), 0);}

		world[0] & v;

		assert( v.size() == 10UL );
		assert( v[9] == T(9) );
	}

	return 0;
} catch(...) {return 1;}


