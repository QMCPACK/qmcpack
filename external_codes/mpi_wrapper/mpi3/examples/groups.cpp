#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 15 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator& world){

	mpi3::communicator prime = world.create(world.group().include( {2, 3, 5, 7, 11, 13} ));

  // If this rank isn't in the new communicator, it will be an invalid (null) communicator
  // Using rank() or size() in a null communicator is erroneous
	cout << "world " << world.rank() << "/" << world.size() << " ---> ";
	if(prime) cout << "prime " << prime.rank() << "/" << prime.size() << std::endl;
	else      cout << "not in prime comm" << std::endl;

	return 0;
}

