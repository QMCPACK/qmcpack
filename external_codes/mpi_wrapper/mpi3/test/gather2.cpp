#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors $0 -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/detail/iterator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	using dd = std::tuple<double, double>;

	std::vector<dd> v_local(10, dd{world.rank(), world.rank() + 1});
	std::vector<dd> v(world.root()?v_local.size()*world.size():0);
	auto last = world.gather(begin(v_local), end(v_local), begin(v));
	if(world.root()){
		assert(last == end(v)); 
		assert(( v[ 0] == dd{0.,1.} ));
		assert(( v[10] == dd{1.,2.} ));
		assert(( v[20] == dd{2.,3.} ));
	}else{
		assert( last == begin(v) );
	}
	return 0;
}

