#if COMPILATION_INSTRUCTIONS// -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-
mpic++ -D_TEST_MPI3_SHARED_COMMUNICATOR -xc++ $0 -o $0x&&mpirun -n 3 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MPI3_SHARED_COMMUNICATOR_HPP
#define MPI3_SHARED_COMMUNICATOR_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/environment.hpp" // processor_name

//#include "/usr/src/kernels/4.18.16-300.fc29.x86_64/include/linux/getcpu.h"
//#include<sched.h>
//#include<numa.h> // sudo dnf install numactl-devel

#include<boost/uuid/uuid.hpp>
#include<boost/uuid/uuid_generators.hpp>

namespace boost{

namespace mpi3{

template<class T = void>
struct shared_window;

struct shared_communicator : communicator{
	shared_communicator() = default;
	shared_communicator(shared_communicator&&) = default;
	shared_communicator(shared_communicator const&) = default;
	shared_communicator(mpi3::group const& g) : communicator(g){}
	shared_communicator(mpi3::group const& g, int tag) : communicator(g, tag){}
private:
	template<class T> static auto data_(T&& t){
		using detail::data;
		return data(std::forward<T>(t));
	}
	template<class T> friend struct shared_window;
	explicit shared_communicator(communicator&& c) : communicator(std::move(c)){}
	explicit shared_communicator(communicator const& comm, int key = 0){
		auto e = static_cast<enum error>(MPI_Comm_split_type(comm.get(), MPI_COMM_TYPE_SHARED, key, MPI_INFO_NULL, &impl_));
		if(e != mpi3::error::success) throw std::system_error{e, "cannot split"};
		name(comm.name()+":"+mpi3::processor_name());
	}
	shared_communicator(communicator const& comm, mpi3::communicator_type t, int key = 0){
		MPI3_CALL(MPI_Comm_split_type)(comm.get(), static_cast<int>(t), key, MPI_INFO_NULL, &impl_);
		boost::uuids::uuid tag = boost::uuids::random_generator{}(); static_assert(sizeof(unsigned int)<=sizeof(boost::uuids::uuid), "!");
		auto utag = reinterpret_cast<unsigned int const&>(tag);
		this->broadcast_n(&utag, 1, 0);
		auto Tag = std::to_string(utag);
		std::string const& base = comm.name();
		// !!! switch-case don't work here because in some MPI impls there are repeats !!!
		if(communicator_type::shared==t){
			#if __linux__
			set_name(base+":shared/pu" + std::to_string(::sched_getcpu())); //same as ::getcpu() // TODO
			#else
			set_name(base+":shared/pu" + Tag);
			#endif
		}
		else if(communicator_type::core     ==t){
			#if __linux__
			set_name(base+":core/pu" + std::to_string(::sched_getcpu())); //same as ::getcpu() // TODO
			#else
			set_name(base+":core/pu" + Tag);
			#endif
		}
		else if(communicator_type::hw_thread==t) set_name(base+":hw_thread"+Tag);
		else if(communicator_type::l1_cache ==t) set_name(base+":l1_cache" +Tag);
		else if(communicator_type::l2_cache ==t) set_name(base+":l2_cache" +Tag);
		else if(communicator_type::l3_cache ==t) set_name(base+":l3_cache" +Tag);
		else if(communicator_type::socket   ==t) set_name(base+":socket"   +Tag);
		else if(communicator_type::numa     ==t) set_name(base+":numa"     +Tag);
		else if(communicator_type::board    ==t) set_name(base+":board"    +Tag);
		else if(communicator_type::host     ==t) set_name(base+":cu"       +Tag);
		else if(communicator_type::cu       ==t) set_name(base+":cu"       +Tag);
		else if(communicator_type::cluster  ==t) set_name(base+":cluster"  +Tag);
	}
	friend class communicator;
public:
	shared_communicator& operator=(shared_communicator const& other) = default;
	shared_communicator& operator=(shared_communicator&& other) = default;
	inline shared_communicator split(int key) const{return split_shared(key);}
	auto split(int color, int key) const{
		return shared_communicator{communicator::split(color, key)};
	}

	template<class T = char>
	shared_window<T> make_shared_window(mpi3::size_t size);
	template<class T = char>
	shared_window<T> make_shared_window();
};

inline shared_communicator communicator::split_shared(int key /*= 0*/) const{
	return shared_communicator(*this, key);
}

inline shared_communicator communicator::split_shared(communicator_type t, int key /*= 0*/) const{
	return shared_communicator(*this, t, key);
}

}}

#ifdef _TEST_MPI3_SHARED_COMMUNICATOR

#include "../mpi3/main.hpp"
#include "../mpi3/operation.hpp"
#include "../mpi3/shared_window.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	
	auto numa = world.split_shared(communicator_type::numa); // fallback to world.split_shared() if OMPI is not available
	auto win = numa.make_shared_window<int>(numa.rank()?0:1);
	assert(win.base() != nullptr and win.size() == 1);
	win.lock_all();
	if(numa.rank() == 0){
		*win.base() = 42;
		win.sync();
	}
	for(int j=1; j != numa.size(); ++j){
		if(numa.rank()==0) numa.send_n((int*)nullptr, 0, j, 666);
		else if(numa.rank()==j) numa.receive_n((int*)nullptr, 0, 0, 666);
	}
	if(numa.rank() != 0) win.sync();
//	int l = *win.base();
	win.unlock_all();

#if 0
	auto win = node.make_shared_window<int>(node.rank()?0:1);
	assert(win.base() != nullptr and win.size() == 1);
	win.lock_all();
	if(node.rank()==0){
		*win.base() = 42;
		win.sync();
	}
	for(int j=1; j != node.size(); ++j){
		if(node.rank()==0) node.send_n((int*)nullptr, 0, j, 666);
		else if(node.rank()==j) node.receive_n((int*)nullptr, 0, 0, 666);
	}
	if(node.rank() != 0) win.sync();
	int l = *win.base();
	win.unlock_all();
	int minmax[2] = {-l,l};
//	node.reduce_in_place_n(&minmax[0], 2, mpi3::max<>{}, 0);
	node.all_reduce_n(&minmax[0], 2, mpi3::max<>{});
	assert( -minmax[0] == minmax[1] );
#endif
	return 0;
}

#endif
#endif

