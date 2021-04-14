#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x -lboost_serialization&&mpirun -n 3 $0x&&rm $0x;exit
#endif
// © Alfredo Correa 2018-2020

#include "../../mpi3/main.hpp"
#include "../../mpi3/environment.hpp"
#include "../../mpi3/shm/mutex.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator){
	return 0;
}
