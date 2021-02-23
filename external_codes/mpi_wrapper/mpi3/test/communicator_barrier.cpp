#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	cout <<"Before barrier, I am "<< world.rank() <<" of "<< world.size() << std::endl;
	world.barrier();
	cout <<"After barrier, I am "<< world.rank() <<" of "<< world.size() << std::endl;
	return 0;
}

