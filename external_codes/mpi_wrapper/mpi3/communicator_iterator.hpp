#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_COMMUNICATOR_ITERATOR $0x.cpp -o $0x.x && time mpirun -n 5 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_COMMUNICATOR_ITERATOR_HPP
#define BOOST_MPI3_COMMUNICATOR_ITERATOR_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

namespace boost{
namespace mpi3{

struct communicator_iterator {
	int n_;
	communicator* commptr_;
	communicator_iterator() : n_(-1), commptr_(nullptr){}
	communicator_iterator(void*** p) : communicator_iterator(){}
//	communicator_iterator(communicator_iterator const& other) = default;
	process operator*() const{return commptr_->operator[](n_);}
	communicator_iterator& operator++(){n_++; return *this;}
	communicator_iterator& operator--(){n_--; return *this;}
	communicator_iterator& operator+=(int n){n_+=n; return *this;}
	communicator_iterator& operator-=(int n){n_-=n; return *this;}
};

communicator_iterator next_periodic(communicator_iterator it){
	it.n_++;
	it.n_ = it.n_ % it.commptr_->size();
	return it;
}

communicator_iterator prior_periodic(communicator_iterator it){
	it.n_--;
	if(it.n_ < 0) it.n_ += it.commptr_->size();
	return it;
}

}}

#ifdef _TEST_BOOST_MPI3_COMMUNICATOR_ITERATOR

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
//	for(auto&& p : world)
//		std::periodic_next(p) << world.rank()*10;
	return 0;
}

#endif
#endif


