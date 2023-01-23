// Copyright 2018-2021 Alfredo A. Correa
#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

auto mpi3::main(int/*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	cout <<"Before barrier, I am "<< world.rank() <<" of "<< world.size() << std::endl;
	world.barrier();
	cout <<"After barrier, I am "<< world.rank() <<" of "<< world.size() << std::endl;
	return 0;
} catch(...) {return 1;}
