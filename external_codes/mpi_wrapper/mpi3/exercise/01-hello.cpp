#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++17 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	if(world.rank() == 0) cout << mpi3::version() << '\n';

	cout << "hello from task " << world.rank() << " in host " << mpi3::processor_name() << std::endl;

	if(world.rank() == 0) cout << "numer of tasks " << world.size() << std::endl;
	return 0;
}

