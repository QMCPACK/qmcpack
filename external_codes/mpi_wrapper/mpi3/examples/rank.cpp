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

	double number = dis(gen);
	std::vector<double> numbers = world.gather_value(number, 0);

	std::vector<int> ranks;
	if(world.rank() == 0){
		assert( not numbers.empty() );
		ranks.resize(numbers.size());
		std::iota(ranks.begin(), ranks.end(), 0);
		std::sort(ranks.begin(), ranks.end(), [&numbers](auto a, auto b){return numbers[a] < numbers[b];});
	} else assert( numbers.empty() );

	int rank = world.scatter_value(ranks);
	cout << "processor " << world.rank() << " has value " << number << " and ranks " << rank << std::endl;
	return 0;
}

