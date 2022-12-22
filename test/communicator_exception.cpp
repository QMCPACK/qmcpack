#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -O3 -Wall -Wextra $0 -o $0x.x -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` && time mpirun -n 3 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() > 1 );
	std::vector<double> v(3);
	switch(world.rank()){
		case 0: v = {3, 4, 5}; try{world.send_n(begin(v), -1, 1);}catch(...){world.send_n(begin(v), 3, 1);} break; // -1 is a bad argument
		case 1:                world.receive(begin(v), 0)   ; break;
	}
	return 0;
}

