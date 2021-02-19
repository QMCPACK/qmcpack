#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_SHARED_MUTEX $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHARED_MUTEX_HPP
#define BOOST_MPI3_SHARED_MUTEX_HPP

#include "../mpi3/shared_window.hpp"
//#include "../mpi3/detail/basic_mutex.hpp"
#include<boost/interprocess/offset_ptr.hpp>

namespace boost{
namespace mpi3{

template<class WindowT>
struct shared_mutex{ //https://gist.github.com/aprell/1486197#file-mpi_mutex-c-L61
	static int tag_counter;
	using flag_t = unsigned char;

	shared_communicator& comm_;
	int rank_; //home
	boost::interprocess::offset_ptr<flag_t> addr_;
	WindowT win_;

//	std::vector<flag_t> wait_list; //[comm_.size()];
	boost::interprocess::offset_ptr<flag_t> wait_list;

//	int tag_;

	shared_mutex(shared_mutex&&) = delete;
	shared_mutex(shared_mutex const&) = delete;
	shared_mutex& operator=(shared_mutex const&) = delete;
	shared_mutex& operator=(shared_mutex&&) = delete;

	shared_mutex(shared_communicator& comm, int rank = 0) : 
		comm_(comm),
		rank_(rank), 
		addr_((comm.rank() == rank)?(flag_t*)mpi3::malloc(sizeof(flag_t)*comm.size()):nullptr),
		win_(addr_, addr_?1:0, comm_),
		wait_list((comm.rank() == rank)?(flag_t*)mpi3::malloc(sizeof(flag_t)*comm.size()):nullptr)//comm.size())
	{
		if(addr_) std::memset(addr_.get(), 0, comm.size());
	//	tag_ = tag_counter;
		++tag_counter;
		comm.barrier();
	}

	void lock(){
		std::fill_n(wait_list, comm_.size(), 0); //	flag_t wait_list[comm_.size()];
		flag_t lock = 1;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank()); //	win_.put_n(&lock, 1, rank_, comm_.rank());
		win_.get_n(wait_list.get(), comm_.size(), rank_);
		win_.unlock(rank_);
		for(int i = 0; i != comm_.size(); ++i){
			if(wait_list[i] == 1 and i != comm_.rank()){
				comm_.receive_n(&lock, 0/*, MPI_ANY_SOURCE, tag_*/); //dummy receive
				break;
			}
		}
	}
	bool try_lock(){
		std::fill_n(wait_list, comm_.size(), 0);
		flag_t lock = 1;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank());
		win_.get_n(wait_list.get(), comm_.size(), rank_);
		win_.unlock(rank_);
		for(int i = 0; i != comm_.size(); ++i)
			if(wait_list[i] == 1 and i != comm_.rank()) return false;
		return true;
	}
	void unlock(){
		std::fill_n(wait_list, comm_.size(), 0);// wait_list.assign(unsigned char(0)); //	flag_t wait_list[comm_.size()];
		flag_t lock = 0;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank());
		win_.get_n(wait_list.get(), comm_.size(), rank_);
		win_.unlock(rank_);

		for(int i = 0; i != comm_.size(); ++i){
			int next = (comm_.rank() + i +1) % comm_.size();
			if(wait_list[next] == 1){
				comm_.send_n(&lock, 0, next/*, tag_*/);
				break;
			}
		}
	}
	~shared_mutex(){
	//	comm_.barrier();
		if(addr_) mpi3::free(addr_.get());
	}
};

}}

#ifdef _TEST_BOOST_MPI3_SHARED_MUTEX

#include "../mpi3/main.hpp"

#include<iostream>
#include<mutex> // lock_guard

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

//	mpi3::shared_mutex m(world);
	{
//		mpi3::shared_mutex& m = *segment.construct<mpi3::shared_mutex>(world);
	}
	
/*	{
		std::lock_guard<mpi3::shared_mutex> lock(m);
		cout << "locked from " << world.rank() << '\n';
		cout << "never interleaved " << world.rank() << '\n';
		cout << "forever blocked " << world.rank() << '\n';
		cout << std::endl;
	}*/
	return 0;
}
#endif
#endif

