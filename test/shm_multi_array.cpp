#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x -lboost_serialization&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/process.hpp"
#include "../../mpi3/shm/allocator.hpp"

#include<random>
#include<thread> //sleep_for
#include<mutex> //lock_guard

#include "../../../boost/multi/array.hpp"

namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;

using std::cout;

namespace boost{
namespace mpi3{
namespace shm{

template<class T, boost::multi::dimensionality_type D>
using array = multi::array<T, D, mpi3::shm::allocator<T>>;

}}}

template<class T> using shm_vector = multi::array<T, 1, mpi3::shm::allocator<double>>;

int mpi3::main(int, char*[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();
{ 
	multi::array<double, 1, mpi3::shm::allocator<double>> V(100, 99., &node);
	assert( V[13] == 99. );

	multi::array<double, 1, mpi3::shm::allocator<double>> W(100, 88., &node);
	assert( W[13] == 88. );
	W = V;
	assert( V[13] == 99. );
}
{
	shm_vector<double> V(200, 99., &node);
	assert( V.size() == size(V) );
	assert( size(V) == 200 );

	cout << V[10] << std::endl;

	double x = 10.0 * V[10];
	node.barrier();
	assert( x == 10.*99. );
}
{
	auto const V = [&]{
		shm_vector<double> V(1000, 0., &node);
		if(node.root()) std::iota(V.begin(), V.end(), 10.);
		node.barrier();
		return V;
	}();
	std::cout
		<< "accumulation result in rank "<< node.rank() 
		<<' '<< std::accumulate(V.begin(), V.begin() + node.rank(), 0.) << std::endl
	;
	node.barrier();
}
{
	multi::array<double, 2, mpi3::shm::allocator<double>> A({5, 5}, 1.); // implicitly uses self communicator
	multi::array<double, 2, mpi3::shm::allocator<double>> B({5, 5}, 2.);

	assert( A[1][2] == 1 );
	assert( B[2][1] == 2 );
	
	B = A;

	assert( B[2][1] == 1 );
	assert( B[3][1] == 1 );
}
{
	mpi3::shm::array<double, 2> A({5, 5}, 1., &node);
	mpi3::shm::array<double, 2> B({5, 5}, 2., &node);

	assert( A[1][2] == 1 );
	assert( B[1][2] == 2 );

	B() = A();

	assert( A[1][2] == 1 );
	assert( B[1][2] == 1 );
}
	return 0;
}

