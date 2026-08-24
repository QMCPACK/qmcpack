// Copyright 2018-2023 Alfredo A. Correa
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>

namespace mpi3 = boost::mpi3;

mpi3::environment const mpienv{mpi3::thread::serialized};  // NOLINT(fuchsia-statically-constructed-objects,cert-err58-cpp)

int main() try {

#ifndef EXAMPI
	BOOST_TEST( mpienv.thread_support() == mpi3::thread::single or mpienv.thread_support() == mpi3::thread::funneled or mpienv.thread_support() == mpi3::thread::serialized );
	BOOST_TEST( mpienv.thread_support() <= mpi3::thread::serialized );
	BOOST_TEST( mpienv.thread_support() <  mpi3::thread::multiple );
#endif

	BOOST_TEST( mpienv.is_thread_main() );

	auto const cuda = mpienv.cuda_support();  // NOLINT(readability-static-accessed-through-instance)
	BOOST_TEST( not cuda );

	return boost::report_errors();
} catch(...) {return 1;}
