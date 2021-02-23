#if COMPILATION_INSTRUCTIONS
mpic++ $0 -c -o $0.o;exit
#endif
// Test for separate compilation / library usage.

#include "../communicator.hpp"
#include "../environment.hpp"

namespace mpi3 = boost::mpi3;

void do_broadcast(mpi3::communicator &c){
	int a = 2;
	c.broadcast_value(a);
}

