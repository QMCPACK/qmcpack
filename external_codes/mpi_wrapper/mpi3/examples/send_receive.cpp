#if COMPILATION_INSTRUCTIONS
export PATH=/home/correaa/prj/alf/boost/mpi3/fake:$PATH
mpic++ -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 1 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	if( world.size() < 2 ) throw std::runtime_error("needs 2 processes");	

	int number = 0;
	if(world.rank() == 0){
		number = -1;
		world[1] << number;
	}else if(world.rank() == 1){
		world[0] >> number;
		cout << "process 1 received number " << number << " from process 0\n";
	}

	return 0;
}

