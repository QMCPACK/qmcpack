#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wfatal-errors $0 -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/environment.hpp"
#include "../../mpi3/shared_window.hpp"
#include "../../mpi3/shm/allocator.hpp"

#include<iostream>

using std::cout;
using std::endl;

namespace mpi3 = boost::mpi3;
namespace shm = mpi3::shm;

int main(int, char*[]){
	mpi3::environment env;
	auto world = env.world();
	mpi3::shared_communicator node = world.split_shared();

	std::vector<double, shm::allocator<double>> v(100, node);
	
	cout << "v = " << v[32] << std::endl;
	return 0;

	cout<<" rank:  " <<world.rank() <<endl;

	shm::allocator<double>   A1(node);
	shm::allocator<float>    A2(node);
	shm::allocator<int>      A3(node);
	shm::allocator<unsigned> A4(node);
	shm::allocator<double>   A5(node);
	shm::allocator<float>    A6(node);
	shm::allocator<int>      A7(node);
	shm::allocator<unsigned> A8(node);
	shm::allocator<double>   A9(node);

	auto data1 = A1.allocate(80);
	auto data2 = A2.allocate(80);
	auto data3 = A3.allocate(80);
	auto data4 = A4.allocate(80);
	auto data5 = A5.allocate(80);
	auto data6 = A6.allocate(80);
	auto data7 = A7.allocate(80);
	auto data8 = A8.allocate(80);
	auto data9 = A9.allocate(80);

//	assert(data9 != nullptr);

	using ptr = decltype(data1);
	std::pointer_traits<ptr>::element_type dd = 5.6;
	std::pointer_traits<ptr>::pointer pp = data1;
//	double* dppp = std::pointer_traits<ptr>::to_address(data1);

	assert(dd == 5.6);

	if(node.root()) data9[3] = 3.4;
	node.barrier();
	assert(data9[3] == 3.4);
	node.barrier();
	A1.deallocate(data1, 80);
	A2.deallocate(data2, 80);
	A3.deallocate(data3, 80);
	A4.deallocate(data4, 80);
	A5.deallocate(data5, 80);
	A6.deallocate(data6, 80);
	A7.deallocate(data7, 80);
	A8.deallocate(data8, 80);
	A9.deallocate(data9, 80);

	

}
