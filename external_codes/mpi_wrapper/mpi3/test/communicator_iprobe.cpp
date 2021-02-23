#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	int send_message = 123;
	int receive_message = 0;

	if(world.rank() == 0){
		mpi3::request r = world.isend(&send_message, &send_message + 1, 0, 0);
		while(not world.iprobe(0, 0) ){};
		assert( world.iprobe(0, 0)->count<int>() );
		world.receive(&receive_message, &receive_message + 1, 0, 0);
		assert( receive_message == send_message );
	//	r.wait(); // wait is called on the desctructor now
	}

	return 0;
}

