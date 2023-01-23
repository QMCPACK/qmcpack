// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;

void syntax_test(mpi3::communicator& world) {
	{
		bool b = world.root();
		world.broadcast_value(b);
		assert( b == true );
	}
	{
		bool b = world.root();
		world.broadcast_n(&b, 1);
		assert( b == true );
	}
	{
		bool b = world.root();
		world[0] || b;
		assert( b == true );
	}
}

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {
	std::vector<std::size_t> sizes = {100L, 64L*1024L};//, 128L*1024L}; // TODO check larger number (fails with openmpi 4.0.5)
	std::size_t NUM_REPS = 5;

	using value_type = int;
	std::vector<value_type> buf(128L*1024L);

	for(std::size_t n=0; n != sizes.size(); ++n) {
	//	if(world.root()) cout<<"bcasting "<< sizes[n] <<" ints "<< NUM_REPS <<" times.\n";

		for(std::size_t reps = 0; reps != NUM_REPS; ++reps) {
			if(world.root()) {
				for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
					buf[i] = static_cast<value_type>(1000000.0 * static_cast<double>(n * NUM_REPS + reps) + static_cast<double>(i));
				}
			} else {
				for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
					buf[i] = static_cast<value_type>(-(n * NUM_REPS + reps) - 1);
				}
			}

			world.broadcast_n(buf.begin(), sizes[n]);
		//	world.broadcast(buf.begin(), buf.begin() + sizes[n], 0);

			for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
				assert( fabs(buf[i] - (1000000.0*static_cast<double>(n * NUM_REPS + reps) + static_cast<double>(i))) < 1e-4 );
			}
		}
	}
	syntax_test(world);
	return 0;
} catch(...) {return 1;}
