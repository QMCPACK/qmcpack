#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 `#-Wall` `#-Wfatal-errors` -I$HOME/prj $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif
// (C) Copyright 2018 Alfredo A. Correa
#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"
#include "alf/boost/mpi3/process.hpp"
#include "alf/boost/mpi3/shm/allocator.hpp"
#include "alf/boost/mpi3/mutex.hpp"

#include<random>
#include<thread> //sleep_for
#include<mutex> //lock_guard

#include "alf/boost/multi/array.hpp"

namespace mpi3 = boost::mpi3;
namespace multi = boost::multi;

using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	mpi3::shared_communicator node = world.split_shared();
	multi::array<double, 2, mpi3::shm::allocator<double>> A({4, 5}, node);

	multi::array<double, 2, mpi3::shm::allocator<double>> B({1, 1}, node);
	
	multi::array<double, 2, mpi3::shm::allocator<double>> C({0, 0}, node);
	
	{
		multi::array<double, 2, mpi3::shm::allocator<double>> D({4, 5}, node);
		D.clear();
	//	multi::array<double, 2, mpi3::shm::allocator<double>> D(std::move(A));
		std::clog << "hola" << std::endl;	
	}
	{
		multi::array<double, 2> a({4, 5});
		a.clear();
	}
//	multi::array<double, 2> d(std::move(a));
	
	return 0;
}

