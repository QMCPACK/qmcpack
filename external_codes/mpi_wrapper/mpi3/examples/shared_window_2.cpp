#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors -Wall $0 -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/mutex.hpp"
#include "alf/boost/mpi3/shared_window.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char *argv[], mpi3::communicator& world){

	mpi3::communicator node = world.split_shared(0);

	mpi3::shared_window win = node.make_shared_window<int>(1);
	win.lock_all();
	win.sync();

	node.barrier();

	std::vector<mpi3::size_t> b_size(node.size());
	std::vector<int> b_disp(node.size());
	std::vector<int*> buf(node.size());

//	mpi3::mutex m(node);
	if(node.rank() != 0){
		for(int i = 0; i != node.size(); ++i){
			b_size[i] = win.size(i);
			b_disp[i] = win.disp_unit(i);
			buf[i] = static_cast<int*>(win.base(i));
			if(i == node.rank()){
				*buf[i] = node.rank() + 100;
				cout << "node rank " << node.rank() << ": *buf[" << i << "] = " << *buf[i] << std::endl;
			}
		}
	}else cout << "node rank " << node.rank() << " master just watches." << std::endl;

	win.unlock_all();

	if(node.rank() != 0)
		for(int i=0; i != node.size(); ++i)
			if(buf[i]) cout << "node rank " << node.rank() << ", target = " << i << ", buf = " << *buf[i] << ", b_size = " << b_size[i] << ", b_disp = " << b_disp[i] << std::endl;

	return 0;
}

