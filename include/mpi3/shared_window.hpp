// -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-
// © Alfredo A. Correa 2018-2020

#ifndef MPI3_SHARED_WINDOW_HPP
#define MPI3_SHARED_WINDOW_HPP

#include "../mpi3/dynamic_window.hpp"
#include "../mpi3/group.hpp"
#include "../mpi3/shared_communicator.hpp"

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

#include<mpi.h>

namespace boost {
namespace mpi3 {

template<class T>
struct shared_window : window<T> {
	shared_window(shared_communicator& comm, mpi3::size_t n, int disp_unit = alignof(T)) //: //sizeof(T)) : // here we assume that disp_unit is used for align
	//	window<T>{}//, comm_{comm}
	{
		void* base_ptr = nullptr;
		MPI_(Win_allocate_shared)(n*static_cast<mpi3::size_t>(sizeof(T)), disp_unit, MPI_INFO_NULL, comm.get(), &base_ptr, &this->impl_);
	}
	explicit shared_window(shared_communicator& comm, int disp_unit = alignof(T)) : 
		shared_window(comm, 0, disp_unit)
	{}


	shared_window(shared_window const&) = default;  // cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3
	shared_window(shared_window &&) noexcept = default;  // cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3

	shared_window& operator=(shared_window const&) = delete;
	shared_window& operator=(shared_window     &&) = delete;
	shared_window& operator=(shared_window      &) = delete;

	~shared_window() = default;

	group get_group() const{group r; MPI_(Win_get_group)(this->impl_, &(r.get())); return r;}
	struct query_t{
		mpi3::size_t size;
		int disp_unit;
		void* base;
	};
	struct query_t query(int rank = MPI_PROC_NULL) const{
		query_t r;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init) delayed init
		MPI3_CALL(MPI_Win_shared_query)(this->impl_, rank, &r.size, &r.disp_unit, &r.base);
		return r;
	}
	template<class TT = T>
	mpi3::size_t size(int rank = 0) const{return query(rank).size/sizeof(TT);}
	int disp_unit(int rank = 0) const{return query(rank).disp_unit;}
	template<class TT = T>
	TT* base(int rank = 0) const{return static_cast<TT*>(query(rank).base);}
};

template<class T /*= char*/> 
shared_window<T> shared_communicator::make_shared_window(mpi3::size_t size){
	return shared_window<T>(*this, size);
}

template<class T /*= char*/>
shared_window<T> shared_communicator::make_shared_window(){
	return shared_window<T>(*this);//, sizeof(T));
}

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_MPI3_SHARED_WINDOW

//#include "../mpi3/main.hpp"
//#include "../mpi3/ostream.hpp"

//namespace mpi3 = boost::mpi3;

//int mpi3::main(int, char*[], mpi3::communicator world){

//	mpi3::ostream wout{world};

//	mpi3::shared_communicator node = world.split_shared();

//	mpi3::shared_window<int> win = node.make_shared_window<int>(node.root()?node.size():0);
//	mpi3::shared_communicator node_cpy{win.get_group()};
////	mpi3::shared_communicator node_cpy{win.get_group(), 0};
//	assert( node_cpy == node );


//	assert( win.base() != nullptr );
//	assert( win.size() == node.size() );

//	win.base()[node.rank()] = node.rank() + 1;
//	node.barrier();
//	for(int i = 0; i != node.size(); ++i) assert(win.base()[i] == i + 1);
//	{
//		mpi3::shared_window<int> win = node.make_shared_window<int>(0);
//	//	assert( mpi3::shared_communicator(win.get_group()) == node );
//	}
//	{
//		mpi3::shared_window<int> win = node.make_shared_window<int>(node.root()?node.size():0);
//		mpi3::shared_window<int> win2 = std::move(win);
//	}
//	return 0;
//}

//#endif
#endif
