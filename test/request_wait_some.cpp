#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	assert( world.size() == 4);

	std::vector<int> buffer(400);
	if(world.root()){
		for(int i = 0; i != world.size(); ++i) buffer[i] = i/100;
		std::vector<mpi3::request> rs(3);
		for(int i = 0; i != world.size() - 1; ++i){
			rs[i] = world.isend(&buffer[i*100], &buffer[i*100] + 100, i + 1, 123);
		}
		int remaining = world.size() - 1;
		while(remaining > 0){
			int count = wait_some(rs.begin(), rs.end()).size(); 
			if(count > 0){
				cout << "sends completed " << count << '\n';
				remaining =- count;
			}
		}
	}else world.receive(buffer.begin(), buffer.begin() + 100, 0, 123);

}

