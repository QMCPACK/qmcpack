#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/process.hpp"

#include <random>

namespace mpi3 = boost::mpi3;
using std::cout;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis;

int mpi3::main(int, char*[], mpi3::communicator& world){

	int elements_per_proc = 5000;

	std::vector<double> v;
	if(world.rank() == 0){
		v.resize(elements_per_proc*world.size(), -1.);
		std::generate(v.begin(), v.end(), [&](){return dis(gen);});
		cout << "serial average " << std::accumulate(v.begin(), v.end(), 0.)/v.size() << '\n';
	}

	std::vector<double> partial_v(elements_per_proc, -1.);
	world.scatter_from(partial_v.begin(), partial_v.end(), v.begin(), 0);
	double partial_ave = std::accumulate(partial_v.begin(), partial_v.end(), 0.)/partial_v.size();
	{
		std::vector<double> partial_aves = world.gather_value(partial_ave, 0); //	std::vector<double> subaves = (world[0] += subave);
		if(world.rank() == 0){
			assert(partial_aves.size() == (unsigned)world.size());
			cout << "parallel average: " << std::accumulate(partial_aves.begin(), partial_aves.end(), 0.)/partial_aves.size() << '\n';
		}else{
			assert(partial_aves.empty());
		}
	}
	{
	//	std::vector<double> partial_aves = world.all_gather_value(partial_ave);
		std::vector<double> partial_aves = (world |= partial_ave);
		assert(partial_aves.size() == (unsigned)world.size());
		cout << "parallel average on " << world.rank() << ": " << std::accumulate(partial_aves.begin(), partial_aves.end(), 0.)/partial_aves.size() << std::endl;
	}
	return 0;
}

