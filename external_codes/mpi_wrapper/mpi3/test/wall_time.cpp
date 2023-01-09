#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x && time mpirun -np 4s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include <chrono>
#include <thread>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], communicator& world){
    using namespace std::chrono_literals;

	auto t1 = mpi3::wall_time();
	std::this_thread::sleep_for(2s);
	auto t2 = mpi3::wall_time();
	cout 
		<< "started at " << t1 << " seconds.\n"
		<< "ended   at " << t2 << " seconds.\n"
		<< "duration   " << t2 - t1 << " seconds.\n"
		<< "resolution " << mpi3::wall_tick() << " seconds.\n"
	;
	return 0;
}

