#if COMPILATION_INSTRUCTIONS
mpic++ -std=c++14 -O3 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && mpirun -n 4 $0x.x $@; rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	world.abort(911);

	return 0;
}

