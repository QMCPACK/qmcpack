#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif
// (C) Copyright 2018 Alfredo A. Correa
#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"
#include "../../mpi3/shm/allocator.hpp"
#include "../../mpi3/mutex.hpp"

#include<random>
#include<thread> //sleep_for
#include<mutex> //lock_guard

#include "../../../boost/multi/array.hpp"

namespace boost{
namespace multi{

template<class It1, class Size, class T>
auto copy_n(
	It1 first, Size n,
	array_iterator<T, 1, boost::mpi3::intranode::array_ptr<T>> d_first
){
	base(d_first).wSP_->fence();
	if(mpi3::group(*base(d_first).wSP_).root()) std::copy_n(first, n, d_first);
	base(d_first).wSP_->fence();
	return d_first + n;
}

template<class RandomAccessIt1, class T>
auto copy(
	RandomAccessIt1 first, RandomAccessIt1 second, 
	array_iterator<T, 1, boost::mpi3::intranode::array_ptr<T>> d_first
){
	using std::distance;
	return copy_n(first, distance(first, second), d_first);
}

}}


namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;

using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();

{
	multi::array<double, 2, mpi3::shm::allocator<double> > A({5, 5}, 1, node);
	multi::array<double, 2, mpi3::shm::allocator<double> > B({5, 5}, 2, node);

	assert( A[1][2] == 1 );
	assert( B[2][1] == 2 );
	
	B = A; // uses copy_n(shm::ptr, shm::ptr, shm::ptr) as optimization

	assert( B[2][1] == 1 );
	assert( B[3][1] == 1 );
}
{
	multi::array<double, 2, mpi3::shm::allocator<double> > A({5, 5}, 1, node);
	multi::array<double, 2, mpi3::shm::allocator<double> > B({50, 25}, 2, node);

	assert( A[1][2] == 1 );
	assert( B[2][1] == 2 );
	
	B = A; // uses uninitialized_copy_n(shm::ptr, shm::ptr, shm::ptr) as optimization

	assert( B[2][1] == 1 );
	assert( B[3][1] == 1 );
}
{
	multi::array<double, 2, mpi3::shm::allocator<double> > A({5, 5}, 1, node);
	multi::array<double, 2, mpi3::shm::allocator<double> > B({5, 5}, 2, node);

	assert( A[1][2] == 1 );
	assert( B[2][1] == 2 );

	B[2] = A.rotated()[1]; // uses copy(multi::it, multi::it, multi::it)

	assert( B[2][1] == 1 );
	assert( B[3][1] == 2 );
}

	return 0;
}

