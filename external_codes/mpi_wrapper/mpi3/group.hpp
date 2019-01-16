#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -std=c++14 -O3 -Wall -Wextra -fmax-errors=2 `#-Wfatal-errors` -D_TEST_MPI3_GROUP $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef MPI3_GROUP_HPP
#define MPI3_GROUP_HPP

#include "../mpi3/detail/iterator_traits.hpp"
#include "../mpi3/detail/strided.hpp"
#include "../mpi3/equality.hpp"
#include "../mpi3/communicator.hpp"
#include "../mpi3/window.hpp"

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<algorithm> // copy_n
#include<numeric> // iota
#include<vector>

namespace boost{
namespace mpi3{

class group{
	MPI_Group impl_ = MPI_GROUP_NULL;
public:
	MPI_Group& operator&(){return impl_;}
	MPI_Group const& operator&() const{return impl_;}
	group() : impl_{MPI_GROUP_EMPTY}{}
	group(mpi3::communicator const& c){
		int s = MPI_Comm_group(&c, &impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot construct group"};
	}
	group(mpi3::window<> const& w){ // TODO window<void>
		int s = MPI_Win_get_group(&w, &impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot get group"};
	}
	group(group const& d) = delete;
/*	{
		std::vector<int> v(d.size());
		std::iota(v.begin(), v.end(), 0);
		int s = MPI_Group_incl(d.impl_, v.size(), v.data(), &impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot copy construct"};
	}*/
	group(group&& g) : impl_{std::exchange(g.impl_, MPI_GROUP_EMPTY)}{}
	group& operator=(group const&) = delete;
	group& operator=(group&& other){
		swap(other);
		other.clear();
		return *this;
	}
	void swap(group& other){std::swap(impl_, other.impl_);}
	void clear(){
		if(impl_ != MPI_GROUP_EMPTY){
			int s = MPI_Group_free(&impl_);
			assert( s == MPI_SUCCESS ); // don't want to throw from ctor
		//	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot free group"};
		}
		impl_ = MPI_GROUP_EMPTY;
	}
	~group(){clear();}
	group include(std::initializer_list<int> il){
		group ret;
		int s = MPI_Group_incl(impl_, il.size(), il.begin(), &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"rank not available"};
		return ret;
	}
	group exclude(std::initializer_list<int> il){
		group ret;
		int s = MPI_Group_excl(impl_, il.size(), il.begin(), &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"rank not available"};
		return ret;
	}
	int rank() const{
		int rank = -1; 
		int s = MPI_Group_rank(impl_, &rank);
		if(s != MPI_SUCCESS) throw std::runtime_error("rank not available");
		return rank;
	}
	bool root() const{return rank() == 0;}
	int size() const{
		int size = -1;
		int s = MPI_Group_size(impl_, &size);
		if(s != MPI_SUCCESS) throw std::runtime_error{"rank not available"};
		return size;
	}
	group sliced(int first, int last, int stride = 1) const{
		int ranges[][3] = {{first, last - 1, stride}};
		group ret;
		int s = MPI_Group_range_incl(impl_, 1, ranges, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot slice"};
		return ret;
	}
	bool empty() const{return size() == 0;}
	friend mpi3::equality compare(group const& self, group const& other){
		int result;
		int s = MPI_Group_compare(self.impl_, other.impl_, &result);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot compare"};
		return static_cast<boost::mpi3::equality>(result);
	}
	bool operator==(group const& other) const{
		mpi3::equality e = compare(*this, other); 
		return e == mpi3::identical or e == mpi3::congruent;
	}
	bool operator!=(group const& other) const{return not operator==(other);}
	bool friend is_permutation(group const& self, group const& other){
		return compare(self, other) != mpi3::unequal;
	}
	friend group set_intersection(group const& self, group const& other){
		group ret;
		int s = MPI_Group_intersection(self.impl_, other.impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot difference"};
		return ret;
	}
	friend group set_difference(group const& self, group const& other){
		group ret;
		int s = MPI_Group_difference(self.impl_, other.impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot difference"};
		return ret;
	}
	friend group set_union(group const& self, group const& other){
		group ret;
		int s = MPI_Group_union(self.impl_, other.impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error{"cannot union"};
		return ret;
	}
	int translate_rank(int rank, group const& other) const{
		int out;
		int s = MPI_Group_translate_ranks(impl_, 1, &rank, other.impl_, &out);
		if(s != MPI_SUCCESS) throw std::runtime_error{"error translating"};
		return out;
	}
};

inline communicator communicator::create(group const& g) const{
	communicator ret;
	int s = MPI_Comm_create(impl_, &g, &ret.impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot crate group"};
	return ret;
}

inline communicator communicator::create_group(class group const& g, int tag = 0) const{
	communicator ret;
	int s = MPI_Comm_create_group(impl_, &g, tag, &ret.impl_);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot create_group"};
	return ret;
}

}}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef _TEST_MPI3_GROUP

#include "../mpi3/main.hpp"

#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	mpi3::communicator w1 = world;
	assert( w1.size() == world.size() );
	assert( w1.rank() == world.rank() );
	mpi3::group g1 = w1;
	assert( g1.rank() == w1.rank() );
	mpi3::communicator w2 = w1.create(g1);
	assert( w2.size() == w1.size() );
	assert( w2.rank() == w1.rank() );
	assert( w2.rank() == world.rank() );
	return 0;
}

#endif
#endif


