#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wfatal-errors -Wall $0 -o $0x.x && time mpirun -np 10 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/shared_window.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char *argv[], mpi3::communicator& world){
	int rank = world.rank();
	int size = world.size();

	mpi3::shared_window sw = world.make_shared_window<int>(world.rank()?0:size);
	int* shared = (int*)sw.base(0);

	if(rank==0){
	//	sw.lock_exclusive(0);
		for(int i = 0; i < size; i++) shared[i] = -1;
	//	sw.unlock(0);
	}

	world.barrier();

	std::vector<int> local(size);

	sw.lock_shared(0);
	for(int i = 0; i != 10; i++) sw.get_n(&local[i], 1, 0, i);
	cout << "processor " << world.rank() << " (before) : ";
	for(int i = 0; i != size; ++i) cout << local[i] << " "; 
	cout << std::endl;
	sw.unlock(0);

	sw.lock_exclusive(0);
	sw.put_n(&rank, 1, 0, rank);
	sw.unlock(0);

	sw.lock_shared(0);
	for(int i = 0; i != 10; i++) sw.get_n(&local[i], 1, 0, i);
	cout << "processor " << world.rank() << " (after) : ";
	for(int i = 0; i != size; ++i) cout << local[i] << ' '; 
	cout << std::endl;
	sw.unlock(0);

	world.barrier();

	return 0;
}

