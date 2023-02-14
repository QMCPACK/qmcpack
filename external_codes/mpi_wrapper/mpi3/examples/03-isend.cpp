#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++17 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	if(world.size() % 2 != 0){
		if(world.root()) cout << "Quitting. Need an even number of tasks: numtasks = " << world.size() << "\n";
		return 1;
	}

	int message = -1;
	cout << "Hello from task " << world.rank() << " on host " << mpi3::processor_name() << "\n";
	if(world.root()) cout << "MASTER: number of mpi tasks is " << world.size() << "\n";

	int partner = world.rank()<world.size()/2?world.rank() + world.size()/2:world.rank()-world.size()/2;

	auto r1 = world.ireceive_value(message, partner, 1);
	auto r2 = world.isend_value(world.rank(), partner, 1);
//	r1.wait();
//	r2.wait();
	mpi3::wait(r1, r2);

	cout << "Task " << world.rank() << " is partner with " << message << std::endl;

}

