#if COMPILATION_INSTRUCTIONS
mpic++ -std=c++14 -O3 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<int> send(3); iota(begin(send), end(send), 0);
	std::vector<int> recv(world.size()*send.size());
	
	world.all_to_all_n(send.begin(), send.size(), recv.begin());
	
	if(world.rank() == 0)
		for(auto& e : recv) std::cout << e << '\n';

	return 0;
}

