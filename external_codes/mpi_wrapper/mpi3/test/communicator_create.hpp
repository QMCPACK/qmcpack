#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x `#-lboost_serialization` && time mpirun -np 15 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::communicator world2 = world;
	mpi3::group world2g = world2.group();
	mpi3::communicator world3 = world2.create(world2g);

	assert(world3.rank() == world.rank());

	return 0;
}

