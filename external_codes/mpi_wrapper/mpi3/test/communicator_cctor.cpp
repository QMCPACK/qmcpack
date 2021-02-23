#if COMPILATION_INSTRUCTIONS
mpic++ -g -O3 -Wall -Wextra $0 -o $0x -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY`#-lboost_serialization` && time mpirun -n 3 valgrind $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#if 1
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

#include<complex>
#include<string>
#include<future>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
{
	mpi3::communicator mty;
	assert( mty.size() == 0 );
	mpi3::communicator mty2 = mty;
	assert( mty2.size() == 0 );
}
{
	std::vector<mpi3::communicator> comms;
	comms.emplace_back(world);
	comms.emplace_back(world);
	comms.emplace_back(world);

	std::vector<mpi3::communicator> comms2; for(auto& e:comms) comms2.emplace_back(e);
}
{
	int const NTHREADS = 10;
	std::vector<std::future<int>> fs;
	for(int i=0; i != NTHREADS; ++i){
#if 0 // this is problematic because copy (mpi_comm_dup) is not thread safe
		fs.emplace_back(std::async([&world](){
			auto comm = world; // hangs here
			return 5;
		}));
#endif
		fs.emplace_back(std::async([thread_comm=world](){
		//	auto comm2 = thread_comm; // works, just to test
			return 5;
		}));
#if 0 // more correct
		fs.emplace_back(std::async([world](){
			auto comm2 = world; // works, just to test
			return 5;
		}));
#endif
		std::cout << "created thread" << std::endl;
	}
	for(int i=0; i != NTHREADS; ++i){
		auto five = fs[i].get();
		assert( five == 5 );
		std::cout << "joined thread" << std::endl;
	}
}
	return 0;
}
#endif

