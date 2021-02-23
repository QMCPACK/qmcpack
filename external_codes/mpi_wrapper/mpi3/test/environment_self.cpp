#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 4 -H localhost,localhost,localhost,localhost,localhost $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#include "../../mpi3/environment.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;

int main(){
	mpi3::environment env;
	cout << "us " << env.get_world_instance().get_attribute_as<int>(mpi3::universe_size) << std::endl;

	return 0;
	auto self = env.self();
	assert( self.size() == 1 );
	assert( self.rank() == 0 );
	cout << "I am process " << self.rank() << " in communicator " << self.name() << std::endl;

	auto world = env.world();
	world.barrier();
	assert( world.size() == 4 );
	assert( world.rank() < 4 );
	cout << "I am process " << world.rank() << " in communicator " << world.name() << std::endl;

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

