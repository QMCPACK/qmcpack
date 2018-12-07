#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/version.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], boost::mpi3::communicator& world){

	if(world.size() % 2 != 0){
		if(world.root()) cout << "Quitting. Need an even number of tasks: numtasks = " << world.size() << "\n";
	}else{
		int message = -1;
		if(world.root()) cout << "MASTER: number of mpi tasks is " << world.size() << "\n";
		if(world.rank() < world.size()/2){
			int partner = world.rank() + world.size()/2;
			int rank = world.rank();
			world.send_value(rank, partner);
			world.receive_value(message, partner);
		}else if(world.rank() >= world.size()/2){
			int partner = world.rank() - world.size()/2;
			world.receive_value(message, partner);
			int rank = world.rank();
			world.send_value(rank, partner);
		}
		cout << "Task " << world.rank() << " is partner with " << message << "\n";
	}

	if(world.rank()==0){
		std::vector<double> code = {1.,2.,3.,4.};
		world.send_n(code.begin(), 4, 1);
	}

	if(world.rank()==1){
		std::vector<double> code(4);
		world.receive_n(code.begin(), 4, 0);
		assert(( code == std::vector<double>{1.,2.,3.,4.} ));
	}

	return 0;
}

