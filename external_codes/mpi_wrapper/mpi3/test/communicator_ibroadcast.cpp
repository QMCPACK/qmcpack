#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert(world.size() > 2);


	mpi3::ostream wout(world);

	std::vector<int> large(10);
	if(world.root())
		iota(large.begin(), large.end(), 0);

	wout << "before:" << std::endl;
	for(auto& e : large) wout << e << " ";
	wout << std::endl;

	{
		auto req = world.ibroadcast(large.begin(), large.end(), 0);
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(5s);
	}

	wout << "after:" << std::endl;
	for(auto& e : large) wout << e << " ";
	wout << std::endl;

	return 0;
}

