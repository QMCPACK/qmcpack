#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -I${HOME}/prj -fpermissive -std=c++17 -Wall `#-Wfatal-errors` -D_TEST_BOOST_MPI3_SHM_MULTI $0x.cpp -o $0x.x -lblas && time mpirun -n 4 $0x.x $@; exit 
#&& rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_MULTI_HPP
#define BOOST_MPI3_SHM_MULTI_HPP

#include "../shm/allocator.hpp"

#ifdef _TEST_BOOST_MPI3_SHM_MULTI

#include "../../mpi3/main.hpp"

#include<alf/boost/multi/array.hpp>
#include<alf/boost/multi/adaptors/blas/dot.hpp>

namespace mpi3 = boost::mpi3; 
namespace multi = boost::multi;

using std::cout;

namespace boost::mpi3::intranode{

template<class Size, class T>
T dot(Size n, array_ptr<T> x, Size incx, array_ptr<T> y, Size incy){
	auto& c = x.wSP_->get_communicator(); assert(c==y.wSP_->get_communicator());
	assert( n % c.size() == 0 );
	auto local = boost::multi::blas::dot(n/c.size(), to_address(x + c.rank()*n/c.size()*incx), incx, to_address(y + c.rank()*n/c.size()*incy), incy);
	return (c += local);
}

}

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();

	multi::array<double, 2, mpi3::shm::allocator<>> M({10, 100000000}, node);
	if(node.root())	std::iota(to_address(M.data()), to_address(M.data() + M.num_elements()), 0.);
	node.barrier();

	auto d = multi::blas::dot(M[5], M[3]);

	if(node.root()) std::cerr << "d = " << d << std::endl;
	return 0;
}

#if 0
#include<iostream>
#include<algorithm> // generate
#include<random>
#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/mutex.hpp"
#include <mutex>
#include<thread> 
#include<chrono>
#include<cassert>

int rand(int lower, int upper){
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_int_distribution<int> uni(lower, upper); 
	return uni(rng);
}
int rand(int upper = RAND_MAX){return rand(0, upper);}

namespace mpi3 = boost::mpi3;
using std::cout; 

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){

	mpi3::shm::managed_shared_memory mpi3mshm(world);

	using elem = double;

	mpi3::shm::vector<elem> v(10, mpi3mshm.get_allocator<elem>());
	assert(not v.empty() and v.front() == 0 and std::equal(std::next(v.begin()), v.end(), v.begin()) );

	mpi3::mutex m(world);
	std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
	{
		std::lock_guard<mpi3::mutex> lock(m); // m.lock();
		for(int i = 0; i != 10; ++i){
			v[i] = world.rank();
			std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
		} //	m.unlock();
	}
	world.barrier();

	if(world.rank() == 0){
		for(int i = 0; i != 10; ++i)
			cout << v[i] << " ";
		cout << std::endl;
	}
	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );

	v.resize(15);

	return 0;
#if 0

	mpi3::mutex m(world);
	std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
	{
		std::lock_guard<boost::mpi3::mutex> lock(m);
	//	m.lock();
		for(int i = 0; i != 10; ++i){
			v[i] = world.rank();
			std::this_thread::sleep_for(std::chrono::milliseconds(rand(10)));
		}
	//	m.unlock();
	}

	world.barrier();

	if(world.rank() == 0){
		for(int i = 0; i != 10; ++i)
			cout << v[i] << " ";
		cout << std::endl;
	}
	assert( std::equal(std::next(v.begin()), v.end(), v.begin()) );
//	v.resize(2);
#endif
}
#endif

#endif
#endif

