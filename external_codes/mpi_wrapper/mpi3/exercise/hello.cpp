#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++17 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"
#include "alf/boost/mpi3/processor_name.hpp"

#include<iostream>

using std::cout;
using std::endl;

int boost::mpi3::main(int argc, char* argv[], boost::mpi3::communicator const& world){
	if(world.rank() == 0) cout << boost::mpi3::version() << '\n';

	cout << "hello from task " << world.rank() << " in host " << boost::mpi3::processor_name() << endl;

	if(world.rank() == 0) cout << "numer of tasks " << world.size() << endl;
}

