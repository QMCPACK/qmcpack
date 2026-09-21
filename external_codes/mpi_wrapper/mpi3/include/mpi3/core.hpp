// Copyright 2018-2025 Alfredo A. Correa

#ifndef BOOST_MPI3_CORE_HPP
#define BOOST_MPI3_CORE_HPP

#include<mpi3/detail/call.hpp>

#include <mpi3/detail/mpi_impl.h>

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
#endif
