#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include "../../mpi3/main.hpp"
#include "../../mpi3/window.hpp"
#include "../../mpi3/group.hpp"

#include<cassert>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char*[], mpi3::communicator world){
	int buf[10];
	mpi3::window<int> win{buf, 10, world}; 
	mpi3::group wing(win);
	mpi3::group g(world);
	assert( g == wing );
	assert( g.rank() == wing.rank() );
	assert( g.rank() == world.rank() );
	assert( mpi3::group(world) == mpi3::group(win) );
	assert( mpi3::group(world).rank() == mpi3::group(win).rank() );
	return 0;
}

