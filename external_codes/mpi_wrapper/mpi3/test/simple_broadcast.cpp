// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iterator>
#include <numeric>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	{
		using T = int;
		std::vector<T> v(10);
		if(world.is_root()) {iota(begin(v), end(v), 0);}  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

		world.broadcast(begin(v), end(v));
		BOOST_TEST( v[9] == T(9) );
	}
	{
		using T = double;
		std::vector<T> v(10);
		if(world.is_root()) {
			iota(begin(v), end(v), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		}
		world.broadcast(begin(v), end(v));

		BOOST_TEST( v[9] == T(9) );
	}

	{
		using T = double;
		std::vector<T> v;
		if(world.is_root()) { v.resize(10); iota(begin(v), end(v), 0); }  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

		auto size = v.size();
		world.broadcast_n(&size, 1);

		v.resize(size);
		BOOST_TEST( v.size() == 10UL );

		world.broadcast_n(v.data(), v.size());

		BOOST_TEST( v[9] == T(9) );
	}
	{
		using T = double;
		std::vector<T> v;
		if(world.is_root()) {v.resize(10); iota(begin(v), end(v), 0);}  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

		world[0] & v;

		BOOST_TEST( v.size() == 10UL );
		BOOST_TEST( v[9] == T(9) );
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
