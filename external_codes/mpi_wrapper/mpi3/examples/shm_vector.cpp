#if COMPILATION_INSTRUCTIONS
time mpicxx -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/shm/vector.hpp"

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/mutex.hpp"

#include<chrono>
#include<mutex> // lock_guard
#include<random>
#include<thread> // sleep_for

int rand(int lower, int upper){
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_int_distribution<int> uni(lower, upper); 
	return uni(rng);
}
int rand(int upper = RAND_MAX){return rand(0, upper);}

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator& world){

	mpi3::shm::managed_shared_memory mpi3mshm(world);

	using elem = double;

	cout << "about to construct" << std::endl;
	mpi3::shm::vector<elem> v(10, mpi3mshm.get_allocator<elem>());
	assert(v.data());
	v.resize(55);
	cout << v.data() << std::endl;
	assert(v.data());
	assert(not v.empty() and v.front() == 0 and std::equal(std::next(v.begin()), v.end(), v.begin()) );
	return 0;

	mpi3::mutex m(world);
	std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
	{
		std::lock_guard<mpi3::mutex> lock(m); // m.lock();
		for(int i = 0; i != v.size(); ++i){
			v[i] = world.rank();
			std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
		} //	m.unlock();
	}
	world.barrier();

	if(world.rank() == 0){
		for(int i = 0; i != v.size(); ++i) cout << v[i] << " ";
		cout << std::endl;
	}
	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );
	cout << "about to resize" << std::endl;
	assert(v.data());
	v.resize(20);
	assert(v.data());
	cout << "after resize" << std::endl;
	world.barrier();
	assert(v.size()==20);
	cout << "before assign" << std::endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
	{
		std::lock_guard<mpi3::mutex> lock(m);
		for(int i = 0; i != v.size(); ++i){
			v[i] = world.rank();
			std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
		}
	}
	cout << "after assign" << std::endl;
	world.barrier();
	cout << "after barrier" << std::endl;

	if(world.rank() == 0){
	//	for(int i = 0; i != v.size(); ++i) cout << v[i] << " ";
		cout << std::endl;
	}
//	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );
	cout << "end" << std::endl;
	return 0;
}

