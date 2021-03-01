#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
mpicxx -x c++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_MAIN_HPP
#define BOOST_MPI3_MAIN_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "../mpi3/environment.hpp"
#include "../mpi3/communicator.hpp"
#include "../mpi3/exception.hpp"

namespace boost{
namespace mpi3{

static int main(int, char*[], boost::mpi3::communicator);

}}

int main(int argc, char* argv[]){
	boost::mpi3::environment env{argc, argv};
	return boost::mpi3::main(argc, argv, env.get_world_instance());
}

#if 0
//#ifndef _BOOST_MPI3_MAIN_ENVIRONMENT
int main(int argc, char* argv[]){
	boost::mpi3::environment env{argc, argv};
	try{
		return boost::mpi3::main(argc, argv, env.world());
	}catch(std::exception& e){
		int rank = -1;
		int s = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(s == MPI_SUCCESS and rank == 0){
			std::cerr << 
				"terminate called after throwing\n"
				"  what(): " << e.what() << std::endl
			;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		return 0;
	}catch(...){
		int rank = -1;
		int s = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if(s == MPI_SUCCESS and rank == 0){
			std::cerr << "terminate called after throwing an unknown type\n";
		}
		return 1;
	}
}
#endif

//int main(int argc, char* argv[]){
//	boost::mpi3::environment env(argc, argv);
//	return boost::mpi3::main(argc, argv, env);
//}
#if 0
#endif

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_MAIN

#include "../mpi3/version.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int boost::mpi3::main(int argc, char* argv[], mpi3::communicator world){
	if(world.rank() == 0) cout<< mpi3::version() <<'\n';
	mpi3::communicator duo = world < 2;
	if(duo) cout <<"my rank in comm "<< duo.name() <<" is "<< duo.rank() <<'\n';
	return 0;
}

#endif
#endif


