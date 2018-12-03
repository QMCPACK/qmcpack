#if COMPILATION_INSTRUCTIONS
time mpicxx -O3 -std=c++14 -Wfatal-errors -Wall $0 -o $0x.x -lboost_system && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/shm/managed_shared_memory.hpp" // there is a bug in boost 1.64, this needs to be included first
#include "alf/boost/mpi3/main.hpp"

#include<atomic>

namespace mpi3 = boost::mpi3;
using std::cout; 

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::communicator node = world.split_shared(0);
	int N = 1000;

	mpi3::shm::managed_shared_memory msm(node, 100000);
	std::atomic<int>& i = *msm.construct<std::atomic<int>>(0);
	node.barrier();
	for(int j = 0; j != N; ++j) ++i;
	node.barrier();
	int snapshot = i;
	if(node.rank() == 0) assert(snapshot == node.size()*N);
	msm.destroy<std::atomic<int>>(i);
	node.barrier();

	return 0;
}

