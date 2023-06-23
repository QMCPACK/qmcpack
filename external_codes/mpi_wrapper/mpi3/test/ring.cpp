// © Alfredo A. Correa 2021
#include "../../mpi3/main.hpp"
#include "../../mpi3/process.hpp"

#include<iostream>
#include<numeric>  // for iota
#include<vector>

namespace bmpi3 = boost::mpi3;

int bmpi3::main(int /*argc*/, char ** /*argv*/, bmpi3::communicator world) try {
	auto const size  = world.size(); assert(size != 0);

	mpi3::process&& next  = world[(world.rank() + size + 1) % size];
	mpi3::process&& prior = world[(world.rank() + size - 1) % size];

	int token;  // NOLINT(cppcoreguidelines-init-variables)
	if(not world.is_root()) {prior >> token;}
	else                    {token = -1;    }

	next << token;

	if(    world.is_root()) {prior >> token;}

	assert(token == -1);
	return 0;
} catch(...) {return 1;}
