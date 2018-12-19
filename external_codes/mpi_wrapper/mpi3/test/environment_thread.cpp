#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 8 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#include "../../mpi3/environment.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;

int main(){
	mpi3::environment env(mpi3::multiple);
	assert( env.is_thread_main() );
	assert( env.query_thread() == mpi3::multiple );
}

