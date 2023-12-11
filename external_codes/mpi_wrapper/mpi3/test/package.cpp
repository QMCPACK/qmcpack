#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wfatal-errors $0 -o $0x.x -lboost_serialization && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
//#include "../../mpi3/process.hpp"

#include "mpi.h"
#include <stdio.h>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char *[], mpi3::communicator world){

	std::vector<char> buffer(110);
	if (world.rank() == 0){
		mpi3::package p(world);
		int i = 123;
		p.pack_n(&i, 1);
		double d = 5.;
		p.pack_n(&d, 1);
		cout << "about to send" << std::endl;
		cout << "size " << p.buffer_.size() << std::endl;
		p.send(1, 12345);
		cout << "already sent" << std::endl;
	}

	if (world.rank() == 1){
		mpi3::package p(world);
		cout << "about to receive" << std::endl;
		p.receive(0, 12345);
		cout << "already received" << std::endl;
		int i = -1;
		p.unpack_n(&i, 1);
		assert( i == 123 );
		double d = -1;
		p.unpack_n(&d, 1);
		assert( d == 5. );
	}

	return 0;
}

