// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_SHM_MUTEX_HPP
#define BOOST_MPI3_SHM_MUTEX_HPP

#include <mpi3/shm/allocator.hpp>

namespace boost {
namespace mpi3 {
namespace shm {

class mutex {
	mpi3::shared_communicator& scomm_;
	using allocator_type = mpi3::shm::allocator<std::atomic_flag>;
	allocator_type alloc_;
	mpi3::shm::ptr<std::atomic_flag> f_;
	public:
	explicit mutex(mpi3::shared_communicator& scomm) : scomm_(scomm), alloc_{std::addressof(scomm_)}, f_(alloc_.allocate(1)){
		if(scomm_.root()) {std::allocator_traits<mpi3::shm::allocator<std::atomic_flag>>::construct(alloc_, &*f_, false);}
		scomm_.barrier();
	}
	mutex(mutex const&) = delete;
	mutex(mutex&&) = delete;

	mutex& operator=(mutex const&) = delete;
	mutex& operator=(mutex     &&) = delete;
	void lock() {
		// spin
		while(f_->test_and_set(std::memory_order_acquire)) {};  // NOLINT(altera-unroll-loops)
	}
	void unlock(){f_->clear(std::memory_order_release);}
	~mutex() {  // noexcept(false) {
		try {
			if(scomm_.root()) {std::allocator_traits<allocator_type>::destroy(alloc_, &*f_);}
			scomm_.barrier();
			alloc_.deallocate(f_, 1);
		} catch(...) {}
	}
};

}  // end namespace shm
}  // end namespace mpi3
}  // end namespace boost

//#if not __INCLUDE_LEVEL__ // TEST_BELOW
//#include "../../mpi3/main.hpp"
//#include<thread> // sleep_for
//#include <mutex> // lock_guard

//namespace mpi3 = boost::mpi3;
//using std::cout; 

//int mpi3::main(int argc, char* argv[], mpi3::communicator world){
//	mpi3::shared_communicator node = world.split_shared();
//	
//	mpi3::shm::mutex m(node);
//	using namespace std::chrono_literals;
//	{
//		std::lock_guard<mpi3::shm::mutex> guard(m);
//		cout << "I am rank "; 
//		std::this_thread::sleep_for(2s);
//		cout << node.rank() << '\n';
//	}

//	return 0;
//}
//#endif
#endif
