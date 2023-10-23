//  -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include <mpi3/environment.hpp>

using std::cout;
namespace mpi3 = boost::mpi3;

int main() {
	assert(not mpi3::initialzed());
	assert(not mpi3::finalized());
	{
		mpi3::environment env;
		assert(mpi3::initialzed());
		assert(not mpi3::finalized());

		cout << "us " << env.get_world_instance().get_attribute_as<int>(mpi3::universe_size) << std::endl;

		return 0;
		auto self = env.self();
		assert(self.size() == 1);
		assert(self.rank() == 0);
		cout << "I am process " << self.rank() << " in communicator " << self.name() << std::endl;

		auto world = env.world();
		world.barrier();
		assert(world.size() == 4);
		assert(world.rank() < 4);
		cout << "I am process " << world.rank() << " in communicator " << world.name() << std::endl;
	}
	assert(mpi3::finalized());

	/* output:
	I am process 0 in communicator MPI_COMM_SELF
	I am process 0 in communicator MPI_COMM_SELF
	I am process 0 in communicator MPI_COMM_SELF
	I am process 0 in communicator MPI_COMM_SELF
	I am process 3 in communicator MPI_COMM_WORLD
	I am process 0 in communicator MPI_COMM_WORLD
	I am process 1 in communicator MPI_COMM_WORLD
	I am process 2 in communicator MPI_COMM_WORLD
	*/
}
