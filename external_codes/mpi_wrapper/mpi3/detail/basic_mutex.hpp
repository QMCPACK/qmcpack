#ifndef BOOST_MPI3_DETAIL_BASIC_MUTEX_HPP
#define BOOST_MPI3_DETAIL_BASIC_MUTEX_HPP

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/allocator.hpp"

namespace boost{
namespace mpi3{
namespace detail{

template<class WindowT>
struct basic_mutex{ //https://gist.github.com/aprell/1486197#file-mpi_mutex-c-L61
	static int tag_counter;
	using flag_t = unsigned char;

	communicator& comm_;
	int rank_; //home
	flag_t* addr_;
	WindowT win_;

	std::vector<flag_t> wait_list; //[comm_.size()];

//	int tag_;

	basic_mutex(basic_mutex&&) = delete;
	basic_mutex(basic_mutex const&) = delete;
	basic_mutex& operator=(basic_mutex const&) = delete;
	basic_mutex& operator=(basic_mutex&&) = delete;

	basic_mutex(communicator& comm, int rank = 0) : 
		comm_(comm),
		rank_(rank), 
		addr_((comm.rank() == rank)?(flag_t*)mpi3::malloc(sizeof(flag_t)*comm.size()):nullptr),
		win_(addr_, addr_?1:0, comm_),
		wait_list(comm.size())
	{
		if(addr_) std::memset(addr_, 0, comm.size());
	//	tag_ = tag_counter;
		++tag_counter;
		comm.barrier();
	}

	void lock(){
		std::fill(wait_list.begin(), wait_list.end(), 0); //	flag_t wait_list[comm_.size()];
		flag_t lock = 1;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank()); //	win_.put_n(&lock, 1, rank_, comm_.rank());
		win_.get_n(wait_list.data(), comm_.size(), rank_);
		win_.unlock(rank_);
		for(int i = 0; i != comm_.size(); ++i){
			if(wait_list[i] == 1 and i != comm_.rank()){
				comm_.receive_n(&lock, 0/*, MPI_ANY_SOURCE, tag_*/); //dummy receive
				break;
			}
		}
	}
	bool try_lock(){
		std::fill(wait_list.begin(), wait_list.end(), 0);
		flag_t lock = 1;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank());
		win_.get_n(wait_list.data(), comm_.size(), rank_);
		win_.unlock(rank_);
		for(int i = 0; i != comm_.size(); ++i)
			if(wait_list[i] == 1 and i != comm_.rank()) return false;
		return true;
	}
	void unlock(){
		std::fill(wait_list.begin(), wait_list.end(), 0);// wait_list.assign(unsigned char(0)); //	flag_t wait_list[comm_.size()];
		flag_t lock = 0;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank());
		win_.get_n(wait_list.data(), comm_.size(), rank_);
		win_.unlock(rank_);

		for(int i = 0; i != comm_.size(); ++i){
			int next = (comm_.rank() + i +1) % comm_.size();
			if(wait_list[next] == 1){
				comm_.send_n(&lock, 0, next/*, tag_*/);
				break;
			}
		}
	}
	~basic_mutex(){
		comm_.barrier();
		if(addr_) mpi3::free(addr_);
	}
};

template<class WindowT>
int basic_mutex<WindowT>::tag_counter = 11023;

}}}

#endif

