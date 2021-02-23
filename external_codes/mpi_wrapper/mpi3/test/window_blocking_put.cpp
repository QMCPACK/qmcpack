#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 `#-Wfatal-errors` $0 -o $0x.x && time mpirun -np 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/window.hpp"

#include<cassert>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	std::vector<double> buf(10);
	for(int i = 0; i != 10; ++i) buf[i] = world.rank()*10.0 + i;
	{
		auto w = world.make_window(buf.data(), buf.size()); 

		if(world.rank() == 0){
			cout << "before the put from proc " << world.rank() << std::endl;
			for(int i = 0; i != 10; ++i) cout << buf[i] << " ";
			cout << std::endl;
		}
		else if(world.rank() == 1){
			w.blocking_put_n(buf.begin(), 10, 0);
		}
		// window destructor (calling Win_free_implicitly syncrhonizes)
		//    alternatively use w.fence();
	}
	if(world.rank() == 0){
		cout << "after the put from proc " << world.rank() << std::endl;
		for(int i = 0; i != 10; ++i){
			cout << buf[i] << " ";
			assert(buf[i] == 10. + i);
		}
		cout << std::endl;
	}

	return 0;
}

