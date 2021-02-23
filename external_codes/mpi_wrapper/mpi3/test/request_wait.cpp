#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	assert( world.size() == 4);

	std::vector<int> buffer(10);
	std::vector<int> buffer2(10); std::iota(buffer2.begin(), buffer2.end(), 0);

	int right = (world.rank() + 1)% world.size();
	int left = world.rank() - 1; if(left < 0) left = world.size() - 1;

	mpi3::request r = world.ireceive(buffer.begin(), buffer.end(), left, 123);
	world.send(buffer2.begin(), buffer2.end(), right, 123);
//	r.wait();
	mpi3::wait(r);
	assert( buffer == buffer2 );

}

