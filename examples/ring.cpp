#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	int token = -1;

	if(world.rank() != 0){
		world[world.rank() - 1] << token;
		cout << "proc " << world.rank() << " received token " << token << " from proc " << world.rank() -1 << '\n'; 
	};

	world[(world.rank() + 1)%world.size()] << token;

	if(world.rank() == 0){
		world[world.size() - 1] >> token;
		cout << "proc " << world.rank() << " received token " << token << " from proc " << world.size() - 1 << '\n';
	}

	return 0;
}

