#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -np 8s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/environment.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/cartesian_communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int main(int argc, char* argv[]){
	mpi3::environment env(argc, argv);
	std::vector<int> dims = {env.world().size()};
	mpi3::cartesian_communicator comm(env.world(), dims);
	{
		std::vector<int> remain(1);
		remain[0] = true;
		mpi3::cartesian_communicator newcomm = comm.sub(remain);
	 	assert( env.self().compare(newcomm) == boost::mpi3::unequal );
	}

	return 0;
}

