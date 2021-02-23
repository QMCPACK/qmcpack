#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++17 `#-Wall` -Wfatal-errors -I$HOME/prj $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/process.hpp"
#include "alf/boost/mpi3/shm/vector.hpp"
#include "alf/boost/mpi3/mutex.hpp"

#include<random>
#include<thread> //sleep_for
#include<mutex> //lock_guard

#include<boost/container/flat_map.hpp>

int rand(int lower, int upper){
	static std::mt19937 rng(std::random_device{}());
	static std::uniform_int_distribution<int> uni(lower, upper); 
	return uni(rng);
}
int rand(int upper = RAND_MAX){return rand(0, upper);}

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

//	mpi3::shm::vector<double>::allocator_type alloc(world);
	mpi3::shm::vector<double> v(10, world);
	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );

	mpi3::mutex m(world);
	std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
	{
		std::lock_guard<mpi3::mutex> lock(m);
	//	m.lock();
		for(int i = 0; i != 10; ++i){
			v[i] = world.rank();
			std::this_thread::sleep_for(std::chrono::milliseconds(rand(10))); // slow down the writing
		}
	//	m.unlock();
	}

	world.barrier();

	if(world.rank() == 0){
		for(int i = 0; i != 10; ++i)
			cout << v[i] << " ";
		cout << std::endl;
	}
	// check that only one process had exclusive access 
	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );

//	if(world.rank() == 0){
		v.resize(2);
//	}
	world.barrier();
	assert( v.size() == 2 );
	return 0;
	{
		boost::container::flat_map<double, double, std::less<double>, std::allocator<std::pair<double, double>>> fm;
		fm[4.1] = 6.;
		assert( fm[4.1] == 6. );
	}
	{
		using shm_flat_map = boost::container::flat_map<double, double, std::less<double>, mpi3::shm::allocator<std::pair<double, double>>>;
		shm_flat_map fm(world);

//		if(world.rank() == 0) fm[4.1] = 6.;
		world.barrier();

//		assert( fm[4.1] == 6. );
	}
	return 0;
}

