#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpicxx -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_COMMUNICATOR_OPERATORS $0x.cpp -o $0x.x && time mpirun -np 3 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_COMMUNICATOR_OPERATORS_HPP
#define BOOST_MPI3_COMMUNICATOR_OPERATORS_HPP

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace boost{
namespace mpi3{
//	communicator operator/(communicator const& comm, int n){
//		int color = comm.rank()*n/comm.size();
//		return comm.split(color);
//	}
}}

#ifdef _TEST_BOOST_MPI3_COMMUNICATOR_OPERATORS
#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	auto hemi = world/2;
	assert(hemi);
	cout 
		<< "I am " << world.name() << " " << world.rank() 
		<< " also I am " << hemi.name() << " " << hemi.rank() 
		<< std::endl
	;
}

#endif
#endif

