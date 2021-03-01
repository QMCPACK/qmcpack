#if COMPILATION_INSTRUCTIONS
mpicxx.mpich -g -std=c++14 -O3 -Wall -Wextra -Wpedantic $0 -o $0x -lboost_serialization&&mpirun.mpich -n 4 valgrind --error-exitcode=1345 $0x&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	using T = double;
	assert( world.size() > 2 );

// initialize local data ///////////////////////////////////////////////////////
	std::vector<T> v_loc;
	switch(world.rank()){
		case 0: v_loc = {0., 0., 0.}        ; break;
		case 1: v_loc = {1., 1., 1., 1.}    ; break;
		case 2: v_loc = {2., 2., 2., 2., 2.}; break;
	}

// gather communication ////////////////////////////////////////////////////////
	std::vector<T> v;
//	v.reserve(v_local.size()*world.size()); // optional! avoid a few allocations
	world.all_gatherv(begin(v_loc), end(v_loc), std::back_inserter(v)); 
//	v.shrink_to_fit();                      // optional! save a few memory

// check communication /////////////////////////////////////////////////////////
	assert((v==std::vector<T>{0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 2.}));
	return 0;
}

