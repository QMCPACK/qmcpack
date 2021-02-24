#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif

#include "../main.hpp"
#include "../communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	world.abort(911);

	return 0;
}

