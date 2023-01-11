// -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-
// © Alfredo A. Correa 2018-2021

#ifndef MPI3_SHARED_COMMUNICATOR_HPP
#define MPI3_SHARED_COMMUNICATOR_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/environment.hpp" // processor_name

//#include "/usr/src/kernels/4.18.16-300.fc29.x86_64/include/linux/getcpu.h"
//#include<sched.h>
//#include<numa.h> // sudo dnf install numactl-devel

#include<boost/uuid/uuid.hpp>
#include<boost/uuid/uuid_generators.hpp>

namespace boost {
namespace mpi3 {

template<class T = void>
struct shared_window;

struct shared_communicator : communicator {
	shared_communicator() = default;

	shared_communicator(shared_communicator const&) = delete;
	shared_communicator(shared_communicator     &&) = default;
	shared_communicator(shared_communicator      &) = default;

	explicit shared_communicator(mpi3::group const& g) : communicator(g) {}
	shared_communicator(mpi3::group const& g, int tag) : communicator(g, tag) {}

 private:
	template<class T> friend struct shared_window;
	explicit shared_communicator(communicator&& c) : communicator(std::move(c)) {}
	explicit shared_communicator(communicator& comm, int key = 0) {
		MPI_(Comm_split_type)(&comm, MPI_COMM_TYPE_SHARED, key, MPI_INFO_NULL, &impl_);
		set_name(comm.name() +":"+ mpi3::processor_name());
	}
	shared_communicator(communicator& comm, mpi3::communicator_type t, int key = 0) {
		MPI_(Comm_split_type)(&comm, static_cast<int>(t), key, MPI_INFO_NULL, &impl_);
		boost::uuids::uuid const tag = boost::uuids::random_generator{}(); static_assert(sizeof(unsigned int)<=sizeof(boost::uuids::uuid), "!");
		auto utag = reinterpret_cast<unsigned int const&>(tag);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa)
		this->broadcast_n(&utag, 1, 0);
		auto Tag = std::to_string(utag);
		std::string const& base_name = comm.name();

		// switch-case doesn't work here because in some MPI impls there are "repeats" in the cases
		if(communicator_type::shared == t) {
			#if __linux__
			set_name(base_name +":shared/pu"+ std::to_string(::sched_getcpu()));  // same as ::getcpu() TODO(correaa)
			#else
			set_name(base_name +":shared/pu"+ Tag);
			#endif
		} else if(communicator_type::core == t) {
			#if __linux__
			set_name(base_name +":core/pu"+ std::to_string(::sched_getcpu()));  // same as ::getcpu() TODO(correaa)
			#else
			set_name(base_name +":core/pu"+ Tag);
			#endif
		}
		else if(communicator_type::hw_thread==t) {set_name(base_name +":hw_thread"+Tag);}
		else if(communicator_type::l1_cache ==t) {set_name(base_name +":l1_cache" +Tag);}
		else if(communicator_type::l2_cache ==t) {set_name(base_name +":l2_cache" +Tag);}
		else if(communicator_type::l3_cache ==t) {set_name(base_name +":l3_cache" +Tag);}
		else if(communicator_type::socket   ==t) {set_name(base_name +":socket"   +Tag);}
		else if(communicator_type::numa     ==t) {set_name(base_name +":numa"     +Tag);}
		else if(communicator_type::board    ==t) {set_name(base_name +":board"    +Tag);}
		else if(communicator_type::host     ==t) {set_name(base_name +":host"     +Tag);}
		else if(communicator_type::cu       ==t) {set_name(base_name +":cu"       +Tag);}
		else if(communicator_type::cluster  ==t) {set_name(base_name +":cluster"  +Tag);}
	}
	friend class communicator;

 public:
	shared_communicator& operator=(shared_communicator const&) = delete;
	shared_communicator& operator=(shared_communicator     &&) = default;
	shared_communicator& operator=(shared_communicator      &) = default;  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)

	~shared_communicator() = default;

	shared_communicator      * operator&()      & {return this;}  // NOLINT(google-runtime-operator)
	shared_communicator const* operator&() const& {return this;}  // NOLINT(google-runtime-operator)
	shared_communicator      * operator&()     && {return this;}  // NOLINT(google-runtime-operator)

	inline shared_communicator split(int key) {return split_shared(key);}
	auto split(int color, int key) {
		return shared_communicator{communicator::split(color, key)};
	}

	template<class T = char>
	shared_window<T> make_shared_window(mpi3::size_t size);
	template<class T = char>
	shared_window<T> make_shared_window();
};

inline shared_communicator communicator::split_shared(int key /*= 0*/) {
	return shared_communicator{*this, key};
}

inline shared_communicator communicator::split_shared(communicator_type t, int key /*= 0*/) {
	return shared_communicator{*this, t, key};
}

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_SHARED_COMMUNICATOR

//#include "../mpi3/main.hpp"
//#include "../mpi3/operation.hpp"
//#include "../mpi3/shared_window.hpp"

//#include<iostream>

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int, char*[], mpi3::communicator world) {
//	auto numa = world.split_shared(communicator_type::numa); // fallback to world.split_shared() if OMPI is not available
//	auto win = numa.make_shared_window<int>(numa.rank()?0:1);
//	assert(win.base() != nullptr and win.size() == 1);
//	win.lock_all();
//	if(numa.rank() == 0) {
//		*win.base() = 42;
//		win.sync();
//	}
//	for(int j=1; j != numa.size(); ++j) {
//		if     (numa.rank()==0) {numa.send_n((int*)nullptr, 0, j, 666);}
//		else if(numa.rank()==j) {numa.receive_n((int*)nullptr, 0, 0, 666);}
//	}
//	if(numa.rank() != 0) {win.sync();}
//	win.unlock_all();

//	return 0;
//}

//#endif
#endif
