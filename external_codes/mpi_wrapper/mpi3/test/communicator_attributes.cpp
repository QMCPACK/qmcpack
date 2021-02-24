#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++17 -Wall -Wextra $0 -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/environment.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	assert(world.get_attribute_as<int>(mpi3::tag_ub) >= 32767);

	int h = world.get_attribute_as<int>(mpi3::host);
	if((h < 0 || h >= world.size()) && h != mpi3::process_null) assert(0);

	int io = world.get_attribute_as<int>(mpi3::io);
	if((io < 0 || io >= world.size()) && io != mpi3::any_source && io != mpi3::process_null) assert(0);

	bool wtime_is_global = world.get_attribute_as<bool>(mpi3::wtime_is_global);
	assert(!wtime_is_global);

//	int* appnum = world.get_attribute_as<int*>(mpi3::application_number);

//	int us = world.get_attribute_as<int>(mpi3::universe_size);
//	(void)us;
//	cout << us << std::endl;
//	assert(us >= world.size());

//	int luc = world.get_attribute_as<int>(mpi3::last_used_code);
//	(void)luc;

	return 0;
}

