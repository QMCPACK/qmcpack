#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() >= 2 );

	auto comm = (world < 2);
	cout <<"First, I am rank "<< world.rank() <<"\n";
	if(comm){
		assert(comm.size() == 2);
		cout <<"Second, I am rank "<< world.rank() <<" in world and "<< comm.rank() <<" in comm\n";
	}else{
		assert(comm.size() == 0);
		assert(comm.empty());
		cout <<"Second, I am rank "<< world.rank() <<" but I am not in comm\n";
	}

	return 0;
}

