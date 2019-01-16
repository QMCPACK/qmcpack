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

	std::vector<double> v(elements_per_proc);
	std::generate(v.begin(), v.end(), [&](){return dis(gen);});

	double local_sum = std::accumulate(v.begin(), v.end(), 0.);
	cout << "local sum for proc " << world.rank() << " is " << local_sum << " average is " << local_sum/v.size() << std::endl;

	{//	reduce
//		double global_sum = world.reduce_value(local_sum, mpi3::sum, 0);
		mpi3::optional<double> global_sum = (world[0] += local_sum);
		if(global_sum){
			assert(world.rank() == 0);
			cout << "total sum is " << *global_sum << ", average is " << *global_sum/(world.size()*elements_per_proc) << std::endl;
		}
	}

	{//	all reduce, all reduce is necessary to compute standard deviations, so everyone knows the mean
		double global_sum = world.all_reduce_value(local_sum, mpi3::sum);
		double mean = global_sum /(elements_per_proc*world.size());
		double local_squares = 0.;
		for(auto& e : v){
			local_squares += (e - mean)*(e - mean);
		}
	//	double global_squares = world.reduce_value(local_squares, mpi3::sum, 0);
		mpi3::optional<double> global_squares = (world[0] += local_squares);
		if(global_squares){
			double stddev = sqrt( *global_squares / elements_per_proc / world.size() );
			cout << "mean = " << mean << ", standard deviation = " << stddev << std::endl;
		}
	}

	return 0;
}

