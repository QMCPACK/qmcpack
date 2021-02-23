#if COMPILATION_INSTRUCTIONS
#export PATH=/home/correaa/prj/alf/boost/mpi3/fake:$PATH
mpic++ $0 -o $0x&&mpirun -n 4 $0x $@&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	cout 
		<< "Hello world from processor '" << mpi3::processor_name() 
		<< "' rank " << world.rank() 
		<< " out of " << world.size() << " processors\n"
	;

	return 0;

}

