#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
OMPI_CXX=$CXX mpic++ $0 -o $0x -lboost_serialization&&mpirun --oversubscribe -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef MPI3_GROUP_HPP
#define MPI3_GROUP_HPP

#include "detail/equality.hpp"

#include "../mpi3/error.hpp"

#include "../mpi3/detail/iterator_traits.hpp"
#include "../mpi3/detail/call.hpp"

#include<cassert>

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

namespace boost{
namespace mpi3{

//class communicator;
//template<class T = void> struct window;

class group{
	MPI_Group impl_ = MPI_GROUP_EMPTY;
public:
	friend class communicator;
	template<class T> friend struct window;
	MPI_Group operator&(){return impl_;}
	MPI_Group& get(){return impl_;}
//	std::pointer_traits<MPI_Group>::element_type const* operator&() const{return impl_;} // this doesn't work because in mpich MPI_Group is not a pointer type
	group() = default;
	group(group&& other) noexcept : impl_{std::exchange(other.impl_, MPI_GROUP_EMPTY)}{}
	group(group const& other){MPI_(Group_excl)(other.impl_, 0, nullptr, &impl_);}
//	explicit group(communicator const& c);//{MPI_Comm_group(c.impl_, &impl_);}
//	explicit group(window<> const& w);
	void swap(group& other) noexcept{std::swap(impl_, other.impl_);}
	group& operator=(group other) noexcept{swap(other); return *this;}
	void clear(){
		if(impl_ != MPI_GROUP_EMPTY) MPI_(Group_free)(&impl_);
		impl_ = MPI_GROUP_EMPTY;
	}
	~group(){if(impl_ != MPI_GROUP_EMPTY) MPI_(Group_free)(&impl_);}
	group include(std::initializer_list<int> il){
		group ret; MPI_(Group_incl)(impl_, il.size(), il.begin(), &ret.impl_); return ret;
	}
	group exclude(std::initializer_list<int> il){
		group ret; MPI_(Group_excl)(impl_, il.size(), il.begin(), &ret.impl_); return ret;
	}
	int rank() const{int rank = -1; MPI_(Group_rank)(impl_, &rank); return rank;}
	bool root() const{assert(not empty()); return rank() == 0;}
	int size() const{int size=-1; MPI_(Group_size)(impl_, &size); return size;}
	group sliced(int first, int last, int stride = 1) const{
		int ranges[][3] = {{first, last - 1, stride}};
		group ret; MPI_(Group_range_incl)(impl_, 1, ranges, &ret.impl_); return ret;
	}
	bool empty() const{return size()==0;}
	friend auto compare(group const& self, group const& other){
		int result; MPI_(Group_compare)(self.impl_, other.impl_, &result);
		return static_cast<mpi3::detail::equality>(result);
	}
	bool operator==(group const& other) const{
		auto e=compare(*this, other); 
		return e == mpi3::detail::identical or e == mpi3::detail::congruent;
	}
	bool operator!=(group const& other) const{return not operator==(other);}
	bool friend is_permutation(group const& self, group const& other){
		return compare(self, other) != mpi3::detail::unequal;
	}
	friend group set_intersection(group const& self, group const& other){
		group ret; MPI_(Group_intersection)(self.impl_, other.impl_, &ret.impl_); return ret;
	}
	friend group set_difference(group const& self, group const& other){
		group ret; MPI_(Group_difference)(self.impl_, other.impl_, &ret.impl_); return ret;
	}
	friend group set_union(group const& self, group const& other){
		group ret; MPI_(Group_union)(self.impl_, other.impl_, &ret.impl_); return ret;
	}
	int translate_rank(int rank, group const& other) const{
		int out; MPI_(Group_translate_ranks)(impl_, 1, &rank, other.impl_, &out); return out;
	}
};

}}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if not __INCLUDE_LEVEL__ // def _TEST_MPI3_GROUP

#include "../mpi3/main.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::communicator w1 = world;
	assert( w1.size() == world.size() );
	assert( w1.rank() == world.rank() );

	mpi3::group g1{w1};
	{
		mpi3::group g2 = g1;
		assert( g1 == g2 );
	}
	assert( g1.rank() == w1.rank() );

	mpi3::communicator w2 = w1.create(g1);
	assert( w2.size() == w1.size() );
	assert( w2.rank() == w1.rank() );
	assert( w2.rank() == world.rank() );

	return 0;
}

#endif
#endif


