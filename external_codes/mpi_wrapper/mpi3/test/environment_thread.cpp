// Copyright 2018-2021 Alfredo A. Correa
#include "../../mpi3/environment.hpp"

namespace mpi3 = boost::mpi3;

int main() try {
	mpi3::environment env{mpi3::thread::serialized};

	assert( env.thread_support() == mpi3::thread::single or env.thread_support() == mpi3::thread::funneled or env.thread_support() == mpi3::thread::serialized );
	assert( env.thread_support() <= mpi3::thread::serialized );
	assert( env.thread_support() <  mpi3::thread::multiple );

	assert( env.is_thread_main() );
} catch(...) {return 0;}
