// Copyright 2019-2025 Alfredo A. Correa

#ifndef BOOST_MPI3_SHARED_MAIN_HPP
#define BOOST_MPI3_SHARED_MAIN_HPP

#include <mpi3/detail/mpi_impl.h>

#include <mpi3/environment.hpp>
#include <mpi3/shared_communicator.hpp>

namespace boost{
namespace mpi3{

int main(int argc, char* argv[], boost::mpi3::shared_communicator node);

}}

int ssmain(int argc, char* argv[]){
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

