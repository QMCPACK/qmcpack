// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>
#include <cmath>
#include <cstddef>
#include <vector>

namespace mpi3 = boost::mpi3;

namespace {

void syntax_test(mpi3::communicator& world) {
	{
		bool b = world.root();
		world.broadcast_value(b);
		BOOST_TEST( b == true );
	}
	{
		bool b = world.root();
		world.broadcast_n(&b, 1);
		BOOST_TEST( b == true );
	}
	{
		bool b = world.root();
		world[0] || b;
		BOOST_TEST( b == true );
	}
}

}  // end namespace

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	std::vector<std::size_t> sizes    = {100L, 64L * 1024L};  //, 128L*1024L}; // TODO check larger number (fails with openmpi 4.0.5)
	std::size_t const NUM_REPS = 5;

	using value_type = int;
	std::vector<value_type> buf(128L * 1024L);

	for(std::size_t n = 0; n != sizes.size(); ++n) {
		//  if(world.root()) cout<<"bcasting "<< sizes[n] <<" ints "<< NUM_REPS <<" times.\n";

		for(std::size_t reps = 0; reps != NUM_REPS; ++reps) {
			if(world.root()) {
				for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
					buf[i] = static_cast<value_type>((1000000.0 * static_cast<double>((n * NUM_REPS) + reps)) + static_cast<double>(i));
				}
			} else {
				for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
					buf[i] = -(static_cast<value_type>((n * NUM_REPS) + reps)) - 1;
				}
			}

			world.broadcast_n(buf.begin(), sizes[n]);
			//  world.broadcast(buf.begin(), buf.begin() + sizes[n], 0);

			for(std::size_t i = 0; i != sizes[n]; ++i) {  // NOLINT(altera-unroll-loops)
				BOOST_TEST( std::abs(buf[i] - (1000000.0*static_cast<double>((n * NUM_REPS) + reps) + static_cast<double>(i))) < 1e-4 );
			}
		}
	}
	syntax_test(world);

	return boost::report_errors();
} catch(...) {
	return 1;
}
