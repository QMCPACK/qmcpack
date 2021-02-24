#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include \""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_MUTEX $0x.cpp -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_MUTEX_HPP
#define BOOST_MPI3_MUTEX_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "../mpi3/window.hpp"
#include "../mpi3/detail/basic_mutex.hpp"

namespace boost{
namespace mpi3{

struct mutex{ //https://gist.github.com/aprell/1486197#file-mpi_mutex-c-L61
	using window_t = mpi3::window<unsigned char>;
	static int tag_counter;
	using flag_t = unsigned char;

	communicator& comm_;
	int rank_; //home
	flag_t* addr_;
	window_t win_;

	std::vector<flag_t> wait_list; //[comm_.size()];

//	int tag_;

	mutex(mutex&&) = delete;
	mutex(mutex const&) = delete;
	mutex& operator=(mutex const&) = delete;
	mutex& operator=(mutex&&) = delete;

	mutex(communicator& comm, int rank = 0) : 
		comm_(comm),
		rank_(rank), 
		addr_((comm.rank() == rank)?(flag_t*)mpi3::malloc(sizeof(flag_t)*comm.size()):nullptr),
		win_(addr_, addr_?comm.size():0, comm_),
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
	~mutex(){
		comm_.barrier();
	//	if(addr_) 
		mpi3::free(addr_);
	}
};

int mutex::tag_counter = 11023;

struct rmutex{ //https://gist.github.com/aprell/1486197#file-mpi_mutex-c-L61
	using flag_t = int;
	using window_t = mpi3::window<flag_t>;
	static int tag_counter;


	communicator& comm_;
	int rank_; //home
	flag_t* addr_;
	window_t win_;

	std::vector<flag_t> wait_list; //[comm_.size()];

//	int tag_;

	rmutex(rmutex&&) = delete;
	rmutex(rmutex const&) = delete;
	rmutex& operator=(rmutex const&) = delete;
	rmutex& operator=(rmutex&&) = delete;

	rmutex(communicator& comm, int rank = 0) : 
		comm_(comm),
		rank_(rank), 
		addr_((comm.rank() == rank)?(flag_t*)mpi3::malloc(sizeof(flag_t)*(comm.size()+2)):nullptr),
		win_(addr_, addr_?(comm.size()+2):0, comm_),
		wait_list(comm_.size())
	{
		win_.fence();
		if(addr_){
			std::memset(addr_ + 2, 0, comm.size());
			addr_[0] = 0;
			addr_[1] = -10000;
			assert(comm_.rank() == 0);
		}
		win_.fence();
	//	tag_ = tag_counter;
		++tag_counter;
		comm_.barrier();
	}

	void lock(){
		win_.lock_exclusive(rank_);
		int count = -1;
		win_.get_n(&count, 1, rank_);
		std::cout << "count " << count << std::endl;
		assert(count == 0);
		{
			int owner = -1;
			win_.get_n(&owner, 1, rank_, 1);
			std::cout << "count is " << count << std::endl;
			std::cout << "onwer is " << owner << std::endl;
			if(count > 0) assert(owner == comm_.rank());
			++count;
			owner = comm_.rank();
		//	win_.put_value(count, rank_, 0);
		//	win_.put_value(owner, rank_, 1);
			if(count > 1) return;
		}
		win_.unlock(rank_);

		std::fill(wait_list.begin(), wait_list.end(), 0); //	flag_t wait_list[comm_.size()];
		flag_t lock = 1;
		win_.lock_exclusive(rank_);
		win_.put_value(lock, rank_, comm_.rank() + 2); //	win_.put_n(&lock, 1, rank_, comm_.rank());
		win_.get_n(wait_list.data(), comm_.size(), rank_, 2);
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
		win_.lock_exclusive(rank_);
		int count;
		win_.get_value(count, rank_, 0);
		if(count > 0){
			int owner;
			win_.get_value(owner, rank_, 1);
			assert(owner == comm_.rank());
			--count;
			win_.put_value(count, rank_, 0);
			if(count > 0) return;
			assert(count == 0);
		}
		win_.unlock(rank_);
		
		std::fill(wait_list.begin(), wait_list.end(), 0);// wait_list.assign(unsigned char(0)); //	flag_t wait_list[comm_.size()];
		flag_t lock = 0;
		win_.lock_exclusive(rank_);
		
		win_.put_value(lock, rank_, comm_.rank() + 2);
		win_.get_n(wait_list.data(), comm_.size(), rank_, 2);
		win_.unlock(rank_);

		for(int i = 0; i != comm_.size(); ++i){
			int next = (comm_.rank() + i +1) % comm_.size();
			if(wait_list[next] == 1){
				comm_.send_n(&lock, 0, next/*, tag_*/);
				break;
			}
		}
	}
	~rmutex(){
		comm_.barrier();
	//	if(addr_) 
		mpi3::free(addr_);
	}
};

int rmutex::tag_counter = 11024;


template<class T>
struct atomic{
	int rank_;
	T* addr_;
	communicator& comm_;
//	T counter = -999;
	window<T> win_;

	atomic(T const& value, communicator& comm, int rank = 0) : 
		comm_(comm),
		rank_(rank),
		addr_(static_cast<T*>(comm.rank()==rank?mpi3::malloc(sizeof(T)):nullptr)),
		win_(addr_, addr_?sizeof(T):0, comm_)
	{
		if(addr_) new (addr_) T(value);
//		if(addr_) *addr_ = value;
	}
	atomic& operator+=(T const& t){
		win_.lock_exclusive(rank_);
		win_.fetch_sum_value(t, *addr_, rank_);
	//	win_.fetch_sum_value(t, counter, rank_);
		win_.unlock(rank_);
		return *this;
	}
	atomic& operator++(){
		win_.lock_exclusive(rank_);
		T t = T(1);
		win_.fetch_sum_value(t, *addr_, rank_);
	//	win_.fetch_sum_value(t, counter, rank_);
		win_.unlock(rank_);
		return *this;
	}
	atomic& operator--(){
		win_.lock_exclusive(rank_);
		T t = -T(1);
		win_.fetch_sum_value(t, *addr_, rank_);
	//	win_.fetch_sum_value(t, counter, rank_);
		win_.unlock(rank_);
		return *this;
	}

	atomic& operator-=(T const& t){return operator+=(-t);}
	atomic& operator*=(T const& t){
		win_.lock_exclusive(rank_);
		win_.fetch_prod_value(t, *addr_, rank_);
	//	win_.fetch_prod_value(t, counter, rank_);
		win_.unlock(rank_);
		return *this;
	}
	atomic& operator/=(T const& t);
	atomic& operator=(T const& t){
		win_.lock_exclusive(rank_);
		win_.fetch_replace_value(t, *addr_, rank_);
	//	win_.fetch_replace_value(t, counter, rank_);
		win_.unlock(rank_);
		return *this;
	}
//	T const& load() const{return counter;}
	operator T(){
		T t;
		win_.lock_exclusive(0);
		win_.put_value(*addr_, 0, comm_.rank());
		win_.get_n(&t, 1, 0);
		win_.unlock(0);
		return t;
	}
	~atomic(){
		comm_.barrier();
		if(addr_) boost::mpi3::free(addr_);
	}
};

template<> atomic<double>& atomic<double>::operator/=(double const& d){return operator*=(1./d);}
template<> atomic<float >& atomic<float >::operator/=(float  const& d){return operator*=(1./d);}

}}

#ifdef _TEST_BOOST_MPI3_MUTEX

#include "../mpi3/main.hpp"

#include<thread>
#include<random>
#include<chrono>

#include<iostream>
#include <mutex>

namespace mpi3 = boost::mpi3; using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	{
		mpi3::mutex m(world);
		{
			m.lock();
			cout << "locked from " << world.rank() << '\n';
			cout << "never interleaved " << world.rank() << '\n';
			cout << "forever blocked " << world.rank() << '\n';
			cout << std::endl;
			m.unlock();
		}
	}
	
	return 0;
	cout << "=================" << std::endl;
	world.barrier();
	{
		mpi3::rmutex m(world);
		{
			m.lock();
		//	m.lock();
			cout << "locked from " << world.rank() << '\n';
			cout << "never interleaved " << world.rank() << '\n';
			cout << "forever blocked " << world.rank() << '\n';
			cout << std::endl;
		//	m.unlock();
			m.unlock();
		}
	}
	cout << "=================" << std::endl;
	world.barrier();
	return 0;


	{
		mpi3::atomic<int> counter(0, world);
		counter += 1;
		if(counter + 1 == world.size()) cout << "Process #" << world.rank() << " did the last updated" << std::endl;
	}
	{
		mpi3::atomic<int> counter(0, world);
		counter += 1;
		world.barrier();
		{
			mpi3::mutex m(world);
			std::lock_guard<mpi3::mutex> lock(m);
			cout << "on process " << world.rank() << " counter = " << (int)counter << std::endl;
		}
	}
	{
		mpi3::mutex m(world);
		std::lock_guard<mpi3::mutex> lock(m);

		cout <<"locked from "<< world.rank() << '\n';
		cout <<"never interleaved" << world.rank() << '\n';
		cout <<"forever blocked "<< world.rank() << '\n';
		cout << std::endl;
	}
	
	cout << "end" << std::endl;
	return 0;
}

#endif
#endif

