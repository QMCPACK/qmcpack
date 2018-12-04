#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/process.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	int const limit = 10;
	if( world.size() != 2 ) throw std::runtime_error("needs 2 processes");

	int partner = (world.rank() + 1)%2;

	for(int count = 0; count != limit; ){
		if(world.rank() == count % 2){
			++count;
			world[partner] << count;
			cout << "proc " << world.rank() << " sent and incremented count " << count << " to proc " << partner << '\n';
		}else{
			world[partner] >> count;
			cout << "proc " << world.rank() << " received count " << count << " from proc " << partner << '\n'; 
		}
	}

	return 0;
}

