#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 4s $0x.x $@ && rm -f $0x.x; exit
#endif

#define BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	mpi3::type t = mpi3::type::int_[100]; // mpi3::type::int_.contiguous(100);
	t.commit_as<std::array<int, 100>>();
	t.commit_as<int[100]>();
//	std::array<int, 100> buffer;
	int buffer[100];

	if(world.rank() == 0) world.send_n(&buffer, 1, 1, 123);
	else if(world.rank() == 1) world.receive_n(&buffer, 1, 0, 123);

	return 0;
}

