#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif
// (C) Copyright 2019 Alfredo A. Correa
#include "../../../mpi3/main.hpp"
#include "../../../mpi3/shm/allocator.hpp"

namespace mpi3 = boost::mpi3;

using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();

	mpi3::shm::allocator<double> A1(&node);
	mpi3::shm::pointer<double> data1 = A1.allocate(80);

	using ptr = decltype(data1);
	std::pointer_traits<ptr>::pointer pp = data1;
	double* dp = std::addressof(*data1);
	double* dp2 = mpi3::shm::pointer_traits<ptr>::to_address(data1);

	if(node.root()) data1[3] = 3.4;
	data1.w_->fence();
	node.barrier();
	assert( *dp == *data1 );
	assert( *dp2 == *data1 );
	assert(data1);
	assert(!!data1);
	assert(not (data1 < data1));
	assert(data1[3] == 3.4);

	A1.deallocate(data1, 80);

	return 0;
}
