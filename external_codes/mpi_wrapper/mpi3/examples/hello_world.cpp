#if COMPILATION_INSTRUCTIONS
export PATH=/home/correaa/prj/alf/boost/mpi3/fake:$PATH
mpic++ -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	cout 
		<< "Hello world from processor '" << mpi3::processor_name() 
		<< "' rank " << world.rank() 
		<< " out of " << world.size() << " processors\n"
	;
//s	assert(world.rank());
	return 0;
}

