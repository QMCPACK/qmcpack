// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <numeric>  // for std::iota
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	if(world.size()%2 == 1) {
		if(world.is_root()) {std::cerr<<"Must be called with an even number of processes\n";}
		return 1;
	}

	{
		std::vector<double> xsend(10); iota(begin(xsend), end(xsend), 0.0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		std::vector<double> xrecv(xsend.size(), -1.);

		auto last = world.send_receive(
			cbegin(xsend), cend(xsend), ((world.rank()/2)*2) + ((world.rank()+1)%2),
			begin(xrecv)
		);

		BOOST_TEST( last == end(xrecv) );
		BOOST_TEST( xrecv[5] == 5.0 );
	}
	{
		std::vector<double> xsend(20); iota(begin(xsend), end(xsend), 100.0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		std::vector<double> xrecv(xsend.size(), -1.);

		world.send(cbegin(xsend), cend(xsend), ((world.rank()/2)*2) + ((world.rank()+1)%2));
		world.receive(begin(xrecv), end(xrecv));

		BOOST_TEST( xrecv[5] == 105.0 );
	}
	{
		std::vector<double> xsend(20); iota(begin(xsend), end(xsend), 100.0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		std::vector<double> xrecv(xsend.size(), -1.0);

		world.send(cbegin(xsend), cend(xsend), ((world.rank()/2)*2) + ((world.rank()+1)%2));
		world.receive(begin(xrecv), end(xrecv));

		BOOST_TEST( xrecv[5] == 105.0 );
	}

	if(world.is_root()) { std::cerr<<"successfully completed\n"; }

	return boost::report_errors();
} catch(...) {
	return 1;
}
