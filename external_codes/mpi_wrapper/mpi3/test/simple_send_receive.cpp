// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/main.hpp"

#include<iostream>
#include<numeric>  // for std::iota
#include<vector>

namespace bmpi3 = boost::mpi3;

int bmpi3::main(int /*argc*/, char ** /*argv*/, bmpi3::communicator world) try {
	if(world.size()%2 == 1) {
		if(world.is_root()) {std::cerr<<"Must be called with an even number of processes"<<std::endl;}
		return 1;
	}

	{
		std::vector<double> xsend(10); iota(begin(xsend), end(xsend), 0.);
		std::vector<double> xrecv(xsend.size(), -1.);

		auto last = world.send_receive(
			cbegin(xsend), cend(xsend), (world.rank()/2)*2 + (world.rank()+1)%2,
			begin(xrecv)
		);

		assert( last == end(xrecv) );
		assert( xrecv[5] == 5. );
	}
	{
		std::vector<double> xsend(20); iota(begin(xsend), end(xsend), 100.);
		std::vector<double> xrecv(xsend.size(), -1.);

		world.send(cbegin(xsend), cend(xsend), (world.rank()/2)*2 + (world.rank()+1)%2);
		world.receive(begin(xrecv), end(xrecv));

		assert( xrecv[5] == 105.0 );
	}
	{
		std::vector<double> xsend(20); iota(begin(xsend), end(xsend), 100.);
		std::vector<double> xrecv(xsend.size(), -1.0);

		world.send(cbegin(xsend), cend(xsend), (world.rank()/2)*2 + (world.rank()+1)%2);
		world.receive(begin(xrecv), end(xrecv));

		assert( xrecv[5] == 105.0 );
	}

	if(world.is_root()) {std::cerr<<"successfully completed"<<std::endl;}
	return 0;
} catch(...) {return 1;}
