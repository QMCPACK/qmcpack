#if COMPILATION_INSTRUCTIONS
(echo "#include<"$0">" > $0x.cpp) && mpicxx -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_SHM_POOL $0x.cpp -o $0x.x -lboost_system && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_POOL_HPP
#define BOOST_MPI3_SHM_POOL_HPP

#include <boost/pool/simple_segregated_storage.hpp>

#include "../../mpi3/shared_window.hpp"
#include<memory>

namespace boost{
namespace mpi3{
namespace shm{

}}}

#ifdef _TEST_BOOST_MPI3_SHM_POOL


#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::shared_window sw = world.make_shared_window<int>(world.rank()?0:10000);
	int* shared = static_cast<int*>(sw.base(0));

	boost::simple_segregated_storage<std::size_t> storage;
	storage.add_block(sw.base(0), sw.size(0), 256);

	int *i; 
	*i = static_cast<char*>(storage.malloc()) - (char*)&i;

	if(world.rank() == 2) *i = 5;
	world.barrier();
	if(world.rank() == 3) assert(*i == 5);
//	sw.unlock();

	return 0;
}
#endif
#endif

