// -*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4; -*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_CORE_HPP
#define BOOST_MPI3_CORE_HPP

#include<mpi3/detail/call.hpp>

#include<mpi.h> // if you get a compilation error here it means that 1) you need to compile or defined your CXX as mpic++ or 2) have not setup the compilation flags to find MPI headers, or 3) not installed an MPI implementation (e.g. `apt install libopenmpi-dev openmpi-bin`)

namespace boost {
namespace mpi3 {

inline bool initialized() { 
	return MPI_(Initialized)();
}

inline bool finalized() {
	return MPI_(Finalized)();
}

}  // end namespace mpi3
}  // end namespace boost

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_ENVIRONMENT

#include<cassert>

namespace mpi3 = boost::mpi3;

int main(){
	assert(not mpi3::initialized() );
}

#endif
#endif