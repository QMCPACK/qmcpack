#if COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
mpicxx -x c++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_MAIN_ENVIRONMENT_HPP
#define BOOST_MPI3_MAIN_ENVIRONMENT_HPP

#ifdef BOOST_MPI3_MAIN_HPP
#error Include either "mpi3/main.hpp" or "mpi3/main_environment.hpp"
#endif

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "../mpi3/environment.hpp"
#include "../mpi3/communicator.hpp"
#include "../mpi3/exception.hpp"

namespace boost{
namespace mpi3{

static int main(int, char*[], boost::mpi3::environment&);

}}

int main(int argc, char* argv[]){
	boost::mpi3::environment env{argc, argv};
	return boost::mpi3::main(argc, argv, env);
}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_MAIN

#include "../mpi3/version.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int boost::mpi3::main(int argc, char* argv[], mpi3::environment& env){
	auto world = env.world();
	if(world.rank() == 0) cout<< mpi3::version() <<'\n';
	mpi3::communicator duo = world < 2;
	if(duo) cout <<"my rank in comm "<< duo.name() <<" is "<< duo.rank() <<'\n';
	return 0;
}

#endif
#endif


