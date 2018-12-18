#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/environment.hpp"
#include "../../mpi3/group.hpp"
#include "../../mpi3/communicator.hpp"

#include "../../mpi3/main.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert(world.size() == 4);

	mpi3::group base = world;
		
	mpi3::communicator w1 = world;
	mpi3::group g1 = w1;
	
	assert( base.size() == world.size() );
	assert( base.rank() == world.rank() );
	
	assert( g1.rank() == base.rank() );

	mpi3::communicator w2 = w1.create(g1);

	assert( w2.size() == w1.size() );
	assert( w2.rank() == w1.rank() );
	assert( w2.rank() == world.rank() );


	mpi3::group g4 = base.exclude({0});
	mpi3::group g5 = base.exclude({1, 2, 3}); 

	mpi3::group g6 = set_union(g5, g4);
	assert( g6 == base);

	mpi3::group g9 = set_intersection(base, g4);
	assert( g4 == g9 );
	
	mpi3::group g44 = g1.sliced(0, g1.size()-1, 2);
	mpi3::group g55 = g1.sliced(1, g1.size(), 2);
	mpi3::group g4455 = set_union(g44, g55);
	assert( g4455.size() == g1.size() );
	assert( is_permutation(g4455, g1) );
	mpi3::group dg4455 = set_intersection(g44, g55);
	assert( dg4455.size() == 0 );
	assert( dg4455.empty() );
	assert( dg4455 == mpi3::group{} );
	
	mpi3::group g11 = world.reversed();
	assert( g11.translate_rank(3, base) == 0 );
	assert( g11.translate_rank(2, base) == 1 );
	assert( g11.translate_rank(1, base) == 2 );
	assert( g11.translate_rank(0, base) == 3 );

	return 0;
}

