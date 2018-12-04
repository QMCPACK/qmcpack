#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra $0 -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	assert( world.size() == 8 );
	
	mpi3::communicator third(world / 3);
	mpi3::communicator third2 = std::move(third);
	assert(not third);

	cout 
		<< "I am rank " << world.rank() << " in " << world.name() << ", "
		<< "I am also " << third2.rank() << " in " << third2.name() << '\n'
	;
	
	mpi3::communicator third3;
	third3 = std::move(third2);

	cout 
		<< "I am rank " << world.rank() << " in " << world.name() << ", "
		<< "I am also " << third3.rank() << " in " << third3.name() << '\n'
	;

	return 0;
}

