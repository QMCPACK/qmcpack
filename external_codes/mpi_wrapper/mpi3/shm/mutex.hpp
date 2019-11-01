#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wfatal-errors -D_TEST_BOOST_MPI3_SHM_MUTEX $0x.cpp -o $0x.x -lrt && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_MUTEX_HPP
#define BOOST_MPI3_SHM_MUTEX_HPP

#include "../../mpi3/shm/allocator.hpp"

namespace boost{
namespace mpi3{
namespace shm{

class mutex{
	mpi3::shared_communicator& scomm_;
	mpi3::shm::allocator<std::atomic_flag> alloc_;//(node);
	mpi3::shm::pointer<std::atomic_flag> f_;
	public:
	mutex(mpi3::shared_communicator& scomm) : scomm_(scomm), alloc_(std::addressof(scomm_)), f_(alloc_.allocate(1)){
		if(scomm_.root()) alloc_.construct(&*f_, false);
		scomm_.barrier();
	}
	mutex(mutex const&) = delete;
	mutex(mutex&&) = delete;
	void lock(){while((*f_).test_and_set(std::memory_order_acquire));}
	void unlock(){(*f_).clear(std::memory_order_release);}
	~mutex(){
		if(scomm_.root()) alloc_.destroy(&*f_);
		scomm_.barrier();
		alloc_.deallocate(f_, 1);
	}
};

}}}

#ifdef _TEST_BOOST_MPI3_SHM_MUTEX

#include "../../mpi3/main.hpp"
#include<thread> // sleep_for
#include <mutex> // lock_guard

namespace mpi3 = boost::mpi3;
using std::cout; 

int mpi3::main(int argc, char* argv[], mpi3::communicator world){
	mpi3::shared_communicator node = world.split_shared();
	
	mpi3::shm::mutex m(node);
	using namespace std::chrono_literals;
	{
		std::lock_guard<mpi3::shm::mutex> guard(m);
		cout << "I am rank "; 
		std::this_thread::sleep_for(2s);
		cout << node.rank() << '\n';
	}

	return 0;
}
#endif
#endif

