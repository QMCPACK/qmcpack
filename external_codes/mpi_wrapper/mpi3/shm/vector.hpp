#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -fpermissive -std=c++14 -Wall -Wfatal-errors -D_TEST_BOOST_MPI3_SHM_VECTOR $0x.cpp -o $0x.x && time mpirun -np 10 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_VECTOR_HPP
#define BOOST_MPI3_SHM_VECTOR_HPP

//#include "../../mpi3/shm/memory.hpp"
//#include<boost/container/vector.hpp>
#include "../shared_window.hpp"
#include "../shm/allocator.hpp"

namespace boost{
namespace mpi3{
namespace shm{

template<class T>
using allocator = boost::mpi3::intranode::allocator<T>;

template<class T>
struct vector : std::vector<T, boost::mpi3::shm::allocator<T>>{
	using std::vector<T, boost::mpi3::shm::allocator<T>>::vector;
	using alloc_traits = std::allocator_traits<boost::mpi3::shm::allocator<T>>;
	vector(std::size_t n) = delete;
	auto& start(){return this->_M_impl._M_start;}
	auto& finish(){return this->_M_impl._M_finish;}
	auto& get_stored_allocator(){return this->_M_impl;}
	vector(std::size_t n, shm::allocator<T> const& alloc) : vector(alloc){
		start() = alloc_traits::allocate(this->_M_impl, n);
		finish() = start() + n;
		uninitialized_construct();
	}
	vector(std::size_t n, const T& value, boost::mpi3::shm::allocator<T> const& alloc) : vector(alloc){
		start() = alloc_traits::allocate(get_stored_allocator(), n);
		finish() = start() + n;
		uninitialized_construct(value);
	}
	template<class It>
	vector(It first, It last, shm::allocator<T> const& alloc) : vector(alloc){
		start() = alloc_traits::allocate(get_stored_allocator(), std::distance(first, last));
		finish() = start() + n;
		uninitialized_copy(first);
	}
	void resize(std::size_t n){
		vector tmp(this->begin(), this->begin() + n, get_stored_allocator());
		swap(start(), tmp.start());
		swap(finish(), tmp.finish());
	}
	~vector(){destroy();}
	protected:
	template<class... Args>
	auto uninitialized_construct(Args&&... args){
		using std::for_each;
		for_each(
			start(), finish(), 
			[&](auto&& e){
				alloc_traits::construct(
					get_stored_allocator(), std::addressof(e), 
					std::forward<Args>(args)...
				);
			}
		);
	}
	template<class It>
	auto uninitialized_copy(It first){
		using std::for_each;
		for_each(
			start(), finish(), 
			[&](auto&& e){
				alloc_traits::construct(
					get_stored_allocator(), std::addressof(e), 
					*first
				);
				++fist;
			}
		);
	}
	void destroy(){
		using std::for_each;
		for_each(
			start(), finish(), 
			[&](auto&& e){
				alloc_traits::destroy(get_stored_allocator(), std::addressof(e));
			}
		);
	}
};

}}}

#ifdef _TEST_BOOST_MPI3_SHM_VECTOR
#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3; 
using std::cout;


int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();
	mpi3::shm::vector<double> v(100, node);
	for(int i = 0; i != node.size(); ++i)
		assert(v[i] == 0.);
	node.barrier();
	v.resize(90);
	node.barrier();
	v[node.rank()] = node.rank()*10.;
	node.barrier();
	node.barrier();
	for(int i = 0; i != node.size(); ++i)
		assert(v[i] == i*10.);
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

