// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>
#include <mpi3/process.hpp>

#include <future>

namespace mpi3 = boost::mpi3;

auto async_send(mpi3::communicator& comm, int val, int target) {  // NOLINT(bugprone-easily-swappable-parameters)
	return std::async([=, &comm] () { // was: [=] () mutable {  in original Joseph's S. posting
    	auto value = val + 1;
		comm[target] << value;  // same as comm.send(&value, &value + 1, target);
	});
}

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {
	if(world.rank() == 0) {
		auto fut = async_send(world, 41, 1);
		fut.wait();
	}

	if(world.rank() == 1) {
		int value{};
	 	world[0] >> value;  // same as world.receive(&value, &value + 1, 0);
		std::cout<< "value = "<< value << std::endl;
		assert( value == 42 );
	}
	return 0;
} catch(...) {return 1;}

