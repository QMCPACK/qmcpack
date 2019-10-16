#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 `#-Wfatal-errors` -D_TEST_BOOST_MPI3_SHM_MEMORY $0x.cpp -o $0x.x -lrt && time mpirun -np 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_MEMORY_HPP
#define BOOST_MPI3_SHM_MEMORY_HPP

#include <boost/pool/simple_segregated_storage.hpp>
#include "../../mpi3/shared_window.hpp"
#include<memory>


namespace boost{
namespace mpi3{
namespace shm{

struct shared_memory_object{
	shared_communicator& comm_;
	std::unique_ptr<mpi3::shared_window> swUP_;
	shared_memory_object(shared_communicator& c) : comm_(c){}
	shared_memory_object(shared_communicator& c, mpi3::size_t n) : comm_(c){
		truncate(n);
	}
	shared_memory_object(shared_memory_object const&) = delete;
	void truncate(mpi3::size_t n){
		swUP_ = std::make_unique<mpi3::shared_window>(comm_.make_shared_window<char>(comm_.rank()==0?n:0));
	}
};

struct managed_shared_memory{
	shared_communicator& comm_;
	shared_window sw_;
	boost::simple_segregated_storage<std::size_t> storage_;
	managed_shared_memory(shared_communicator& comm, mpi3::size_t n) : 
		comm_(comm), 
		sw_(comm.make_shared_window<char>(comm.rank()?0:n))
	{
	//	if(comm_.rank()==0) 
		storage_.add_block(sw_.base(0), sw_.size(0), n);
	}
	void* malloc(){return storage_.malloc();}
#if 0
	array_ptr<void> allocate(mpi3::size_t n){
		array_ptr<void> ret;
		ret.swP_ = new shared_window(comm_.make_shared_window<char>(comm_.rank()==0?n:0)); //(comm_.rank()==0?n:0);
	//	ret.swSP_ = std::make_shared<shared_window>(comm_.make_shared_window<char>(comm_.rank()==0?n:0));
		return ret;
	}
	void deallocate(array_ptr<void> ptr){
		if(not ptr.swP_) return;
		delete ptr.swP_;
	}
	template<class T>
	allocator<T> get_allocator();
#endif
//	template<class T, class... Args>
//	T* construct(Args&&... args){
//		if(comm_.rank()==0){
//			T* p = storage_.construct<T>(std::forward<Args>(args)...);
//		}
//	}
	
};


template<class T>
allocator<T> managed_shared_memory::get_allocator(){
	return allocator<T>(*this);
}

struct mapped_region{
	shared_memory_object& smo_;
	mapped_region(shared_memory_object& smo) : smo_(smo){}
	mapped_region(mapped_region const&) = delete;
	void* get_address(){return smo_.swUP_->base(0);}
	mpi3::size_t get_size(){return smo_.swUP_->size(0);}
//	mpi3::size_t get_free_memory() const;
//	void zero_free_memory() const;
//	bool all_memory_deallocated();
//	bool check_sanity();
};

}}}


namespace boost{
namespace mpi3{
namespace shm{

template<class Alloc, typename T, typename Size, typename TT>
	pointer<T> uninitialized_fill_n(Alloc const& a, pointer<T> f, Size n, TT const& val){
	assert(0);
	using std::uninitialized_fill_n;
	if(mpi3::group(*f.wSP_).root()) uninitialized_fill_n(to_address(f), n, val);
	first.wSP_->fence();
	first.wSP_->fence();
	return first + n;
}

}}}

#ifdef _TEST_BOOST_MPI3_SHM_MEMORY

#include "../../mpi3/main.hpp"
#include "../../mpi3/mutex.hpp"

namespace mpi3 = boost::mpi3;
using std::cout; 

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){




	{
		mpi3::shm::shared_memory_object mpi3shm(world);
		mpi3shm.truncate(100);
		mpi3::shm::mapped_region mpi3region(mpi3shm);
		if(world.rank() == 0){
			std::memset(mpi3region.get_address(), 1, mpi3region.get_size());
		}
		world.barrier();
		if(world.rank() == 1){
			char* mem = static_cast<char*>(mpi3region.get_address());
			for(int i = 0; i != mpi3region.get_size(); ++i) assert(mem[i] == 1);
		}
	}
	{
		mpi3::shm::managed_shared_memory mpi3mshm(world);

		mpi3::shm::array_ptr<void> ptr = mpi3mshm.allocate(100);
		mpi3mshm.deallocate(ptr);

		ptr = mpi3mshm.allocate(200);
		mpi3mshm.deallocate(ptr);

		{
			mpi3::shm::allocator<double> a = mpi3mshm.get_allocator<double>();
			mpi3::shm::array_ptr<double> ptr = a.allocate(10);
			if(world.rank() == 0) std::fill_n(ptr, 10, 5);
			world.barrier();
			if(world.rank() == 1) for(int i = 0; i != 10; ++i) assert(ptr[i] == 5);
			a.deallocate(ptr);
		}
		{
			std::allocator<int> a;
			int* ptr = a.allocate(100);
			a.deallocate(ptr, 100);
		}
	}

	world.barrier();

	{
		mpi3::shm::managed_shared_memory mpi3mshm(world);
		std::atomic<int>& atomic = *mpi3mshm.construct<std::atomic<int>>(0);
	}

	return 0;
}
#endif
#endif

