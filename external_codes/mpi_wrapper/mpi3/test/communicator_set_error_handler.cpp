#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/error_handler.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	world.set_error_handler(mpi3::error_handler::code); // default, internal function returns codes
	double d = 5.;
	try{
		world.send(&d, &d + 1, 100);
	}catch(...){
		cout << "catched exception" << std::endl;
		return 0;
	}

//	world.set_error_handler(mpi3::error_handler::fatal); // fail immediately 
//	world.send(&d, &d + 1, 100);

	return 1;
}

