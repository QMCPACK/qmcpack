#if COMPILATION_INSTRUCTIONS/*-*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-*/
mpic++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#include "../../mpi3/environment.hpp"
#include "../../mpi3/group.hpp"
#include "../../mpi3/communicator.hpp"

#include "../../mpi3/main.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	
	mpi3::group wg{world};
	mpi3::communicator w2 = wg;
	assert( w2.rank() == world.rank() );
	assert( w2.size() == world.size() );


	mpi3::communicator half = world/2;
	mpi3::group hg{half};

	mpi3::communicator h2 = hg;
	assert(half.rank() == h2.rank());

	return 0;
}

