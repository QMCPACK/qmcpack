// Copyright 2022 Alfredo A. Correa
#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

#include<chrono>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	auto const my_rank = world.rank();
	mpi3::request r = world.ibarrier();

	using namespace std::literals::chrono_literals;
	std::this_thread::sleep_for(2s);

	std::cout<<"mpi process "<< my_rank <<" call ibarrier."<< std::endl;
	r.wait();
	assert( r.completed() );
	std::cout<<"mpi process "<< my_rank <<" the barrier is complete."<< std::endl;
	return 0;
} catch(...) {return 1;}
