#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++17 -I${HOME}/prj/ -Wfatal-errors $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/environment.hpp"
#include "alf/boost/mpi3/shared_window.hpp"

#include<iostream>

using std::cout;
using std::endl;

namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]){
	MPI_Init(&argc, &argv);
	MPI_Comm node;
	int s = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
	int rank = -1;
	MPI_Comm_rank(node, &rank);
	assert(s == MPI_SUCCESS);

#define ALLOC(N) \
	MPI_Win win##N; \
	void* base_ptr##N = nullptr; \
	int s##N = MPI_Win_allocate_shared((rank==0?80:0)*sizeof(char), 1, MPI_INFO_NULL, node, &base_ptr##N, &win##N); \
	MPI_Aint size##N = -1; \
	int i##N = -1; \
	char* ptr##N = nullptr; \
	MPI_Win_shared_query(win##N, 0, &size##N, &i##N, &ptr##N); \

ALLOC(1);
ALLOC(2);
ALLOC(3);
ALLOC(4);
ALLOC(5);
ALLOC(6);
ALLOC(7);
ALLOC(8);
ALLOC(9);

	MPI_Barrier(node);

#define CHECK(N) \
	if(rank == 0) ptr##N[3] = 'z'; \
	MPI_Barrier(node); \
	assert(ptr##N[3] == 'z');
//    mpi3::environment env(argc, argv);
//    auto node = env.world();
//    mpi3::shared_communicator node = world.split_shared();
//    cout<<" rank:  " <<world.rank() <<endl;

CHECK(9);

#define DEALLOC(N) \
	if(win##N != MPI_WIN_NULL) MPI_Win_free(&win##N);

DEALLOC(1);

	MPI_Finalize();
	return 0;
#if 0
    mpi3::intranode::allocator<char>   A1(node);

    auto data1 = A1.allocate(80);
    if(node.root()) data1[3] = 'z';
    node.barrier();
    assert(data1[3] == 'z');
    node.barrier();
    A1.deallocate(data1, 80);
#endif
}

