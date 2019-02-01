#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -np 4s $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/environment.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/request.hpp"

#include<thread> // sleep_for

namespace mpi3 = boost::mpi3;
using std::cout;
int main(){
	mpi3::environment env;
	auto& world = env.world();
	assert( world.size() == 4 );

	if(world.rank() == 0){
		int rem = 3;
		std::vector<mpi3::request> rs(3);
		for(int i = 1; i != world.size(); ++i){
			std::vector<int> buffer(100);
			rs[i - 1] = world.ireceive(buffer.begin() + i, buffer.begin() + i + 1, i);
		}
		while(rem > 0){
			std::vector<int> completed = mpi3::completed_some(rs);
			if(completed.size() > 0){
				cout << "finished " << completed.size() << "\n";
				rem -= completed.size();
			}else std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}else{
		std::vector<int> buffer(100);
		world.send(buffer.begin(), buffer.begin() + 1, 0);
	}
}


