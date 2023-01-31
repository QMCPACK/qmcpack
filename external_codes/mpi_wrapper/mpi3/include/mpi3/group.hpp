// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MPI3_GROUP_HPP
#define MPI3_GROUP_HPP

#include "detail/equality.hpp"

#include "../mpi3/error.hpp"

#include "../mpi3/detail/call.hpp"
#include "../mpi3/detail/iterator_traits.hpp"

#include<cassert>

// #define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

namespace boost {
namespace mpi3 {

using ptr = MPI_Group;

class group {
	MPI_Group impl_ = MPI_GROUP_EMPTY;

 public:
	friend class communicator;
	template<class T> friend class window;

	MPI_Group& get() {return impl_;}
	MPI_Group operator&() {return get();}  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions,google-runtime-operator) access implementation as pointer

//	std::pointer_traits<MPI_Group>::element_type const* operator&() const{return impl_;} // this doesn't work because in mpich MPI_Group is not really pointer type

	group() = default;
	group(group&& other) noexcept : impl_{std::exchange(other.impl_, MPI_GROUP_EMPTY)}{}
	group(group const& other){MPI_(Group_excl)(other.impl_, 0, nullptr, &impl_);}

	void swap(group& other) noexcept{std::swap(impl_, other.impl_);}
//	group& operator=(group other) noexcept{swap(other); return *this;}
	group& operator=(group const& other) {group tmp(other); swap(tmp)  ; return *this;}
	group& operator=(group     && other) noexcept {         swap(other); return *this;}

	void clear(){
		if(impl_ != MPI_GROUP_EMPTY) {
			try {
				MPI_(Group_free)(&impl_);
			} catch(...) {}
		}
		impl_ = MPI_GROUP_EMPTY;
	}
	~group(){
		if(impl_ != MPI_GROUP_EMPTY) {
			try {MPI_(Group_free)(&impl_);} catch(...) {}
		}
	}

	group include(std::initializer_list<int> il) const {
		group ret; MPI_(Group_incl)(impl_, static_cast<int>(il.size()), il.begin(), &ret.impl_); return ret;
	}
	group exclude(std::initializer_list<int> il) const {
		group ret; MPI_(Group_excl)(impl_, static_cast<int>(il.size()), il.begin(), &ret.impl_); return ret;  // TODO(correaa) use safe cast
	}
	auto rank() const -> int {int rank = -1; MPI_(Group_rank)(impl_, &rank); return rank;}
	auto root() const -> bool {assert(not empty()); return rank() == 0;}
	auto size() const -> int {int size = -1; MPI_(Group_size)(impl_, &size); return size;}

	group sliced(int first, int last, int stride = 1) const {
		int ranges[][3] = {{first, last - 1, stride}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		group ret; MPI_(Group_range_incl)(impl_, 1, ranges, &ret.impl_); return ret;
	}

	bool empty() const {return size()==0;}
	friend auto compare(group const& self, group const& other){
		int result;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_(Group_compare)(self.impl_, other.impl_, &result);
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
		int out;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_(Group_translate_ranks)(impl_, 1, &rank, other.impl_, &out); 
		return out;
	}
};

}  // end namespace mpi3
}  // end namespace boost

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


