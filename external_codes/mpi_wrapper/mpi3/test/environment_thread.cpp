// Copyright 2018-2021 Alfredo A. Correa
#include "../../mpi3/environment.hpp"

namespace mpi3 = boost::mpi3;

mpi3::environment const mpienv{mpi3::thread::serialized};  // NOLINT(fuchsia-statically-constructed-objects,cert-err58-cpp)

int main() try {

	assert( mpienv.thread_support() == mpi3::thread::single or mpienv.thread_support() == mpi3::thread::funneled or mpienv.thread_support() == mpi3::thread::serialized );
	assert( mpienv.thread_support() <= mpi3::thread::serialized );
	assert( mpienv.thread_support() <  mpi3::thread::multiple );

	assert( mpienv.is_thread_main() );
} catch(...) {return 0;}
