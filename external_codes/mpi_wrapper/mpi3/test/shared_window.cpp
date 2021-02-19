#if COMPILATION_INSTRUCTIONS// -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-
mpicxx -O3 -std=c++14 $0 -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/shared_window.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char *argv[], mpi3::communicator world){

	auto node = world.split_shared();
	int n = 100;
	auto sw = node.make_shared_window<double>(node.root()?n:0);
	assert( sw.size(0) == n );
	assert( sw.size(1) == 0 );

//	mpi3::shared_window<double> sw(node, node.rank()==0?n:0);


//	int* arr = (int*)sw.base(0);
#if 0
//	std::cout<<"size "<< sw.size(0) <<std::endl;
//	assert( n == sw.size(0)/sizeof(double) );
//	int size = sw.size(0)/sizeof(double);

	cout << "In rank " << world.rank() << " baseptr = " << arr << std::endl;

	for(int i = world.rank(); i < n; i+=world.size())
		arr[i] = world.rank();

	world.barrier();

	if(world.rank() == 0){
		for(int i = 0; i < n; i++) cout << arr[i];
		cout << '\n';
	}
#endif
	return 0;
}

