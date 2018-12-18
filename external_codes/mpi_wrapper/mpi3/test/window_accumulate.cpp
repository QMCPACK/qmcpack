#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 2 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include "../../mpi3/main.hpp"
#include "../../mpi3/window.hpp"

#include<cassert>

namespace mpi3 = boost::mpi3;
using std::cout;

#define NROWS 100
#define NCOLS 100

int mpi3::main(int, char*[], mpi3::communicator world){

	double A[NROWS][NCOLS];

	if(world.root()){
		for(int i = 0; i != NROWS; ++i)
			for(int j = 0; j != NCOLS; ++j)
				A[i][j] = i*NCOLS + j;

		mpi3::window<> w = world.make_window();
		w.fence(); // note the two fences here
		w.accumulate_n( (double*)A, NROWS*NCOLS, 1 );
		w.fence();
	}else{
		for(int i = 0; i != NROWS; ++i)
			for(int j = 0; j != NCOLS; ++j)
				A[i][j] = i*NCOLS + j;
		mpi3::window<> w = world.make_window( (double*)A, NROWS*NCOLS );
		w.fence(); // note the two fences here
		w.fence();
		for(int i = 0; i != NROWS; ++i){
			for(int j = 0; j != NCOLS; ++j){
				if(world.rank() == 1) assert( A[i][j] == (i*NCOLS + j)*2 );
				else assert( A[i][j] == i*NCOLS + j );
			}
		}
	}

	return 0;
}

