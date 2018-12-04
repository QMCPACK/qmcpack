#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wfatal-errors -Wall -Wextra -Wunused-result $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() == 2);

	int right = world.rank() + 1; if(right >= world.size()) right -= world.size();
	int left  = world.rank() - 1; if(left  <  0           ) left  += world.size();

	std::vector<int> buffer1(10); 
	std::vector<int> buffer2(10); 
	iota(begin(buffer2), end(buffer2), 0);

	mpi3::request r1 = world.ireceive_n(buffer1.begin(), buffer1.size(), left , 123);
	mpi3::request r2 = world.issend_n  (buffer2.begin(), buffer2.size(), right, 123);
	r1.wait();
	r2.wait();

	assert( buffer1 == buffer2 );

	return 0;
}

