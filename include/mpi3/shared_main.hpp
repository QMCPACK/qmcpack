#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include<"$0">" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_SHARED_MAIN $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHARED_MAIN_HPP
#define BOOST_MPI3_SHARED_MAIN_HPP

#include<mpi.h>

#include "../mpi3/environment.hpp"
#include "../mpi3/shared_communicator.hpp"

namespace boost{
namespace mpi3{

int main(int argc, char* argv[], boost::mpi3::shared_communicator node);

}}

int main(int argc, char* argv[]){
	boost::mpi3::environment env(argc, argv);
	return boost::mpi3::main(argc, argv, env.world().split_shared());
}

#ifdef _TEST_BOOST_MPI3_SHARED_MAIN

#include "../mpi3/version.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::shared_communicator& node){
	if(node.rank() == 0) cout << mpi3::version() << '\n';
	return 0;
}

#endif
#endif

