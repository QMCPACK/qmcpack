#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && mpirun -np 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/shm/memory.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::shm::shared_memory_object shm(world);
	shm.truncate(1000); // collective
	mpi3::shm::mapped_region region(shm);
	if(world.rank() == 0){
	//	mpi3::shm::mapped_region region(shm);
		std::memset(region.get_address(), 1, region.get_size());
	}
	world.barrier();
	if(world.rank() == 1){
	//	mpi3::shm::mapped_region region(shm);
		char* mem = static_cast<char*>(region.get_address());
		for(std::size_t i = 0; i != region.get_size(); ++i){
			if(*mem++ != 1) assert(0);
		}
	}

	return 0;
}

