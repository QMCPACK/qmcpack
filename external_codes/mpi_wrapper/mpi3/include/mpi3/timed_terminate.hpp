// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_TIMED_TERMINATE_HPP
#define BOOST_MPI3_TIMED_TERMINATE_HPP

#include "../mpi3/environment.hpp"

namespace boost {
namespace mpi3 {

template<class Duration>
[[noreturn]] void timed_terminate(Duration d, mpi3::communicator& comm = mpi3::environment::get_world_instance()) {
	auto rbarrier = comm.ibarrier();
	auto const t0 = mpi3::wall_time();
	// now spin  
	while(not rbarrier.completed() and (mpi3::wall_time() - t0) < d) {}  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops) investigate (w/clang-tidy 14)

	if(rbarrier.completed()) {
		if(comm.root()) {comm.abort(911);}
	} else {
		std::cout<<"MPI program terminated from rank "<< comm.rank() <<" after timeout of "<< std::chrono::duration_cast<std::chrono::seconds>(d).count() <<" seconds, not all processes failed within that time."<<std::endl;
	}
	comm.abort(911);
	// never call std::terminate from here (it will recurse)
	std::abort();  // necessary to avoid error for returning in a [[noreturn]] function
}

}  // end namespace mpi3
}  // end namespace boost

#endif
