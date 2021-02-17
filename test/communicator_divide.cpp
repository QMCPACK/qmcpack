#if COMPILATION_INSTRUCTIONS
mpic++  -Wall -Wextra $0 -o $0x &&mpirun -n 6 valgrind --suppressions=communicator_main.cpp.openmpi.supp $0x&&rm $0x;exit
#endif
// © Alfredo Correa 2018-2020

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	assert( world.size() == 6 );

	mpi3::communicator fifth = world/5;

	cout << "I am rank " << world.rank() << " in " << world.name() << ", ";

	if(fifth) cout <<"I am also "   << fifth.rank() << " in " << fifth.name() << '\n';
	else      cout <<"I am not in " << fifth.name() << '\n';

	return 0;
}

