// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION

#include "../main.hpp"
#include "../communicator.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {
	{
		mpi3::type t = mpi3::int_;
		assert( t.size() == sizeof(int) );
		assert( t.extent() == sizeof(int) );
		assert( t.lower_bound() == 0 );
		assert( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type t = mpi3::int_[3];
		assert( t.size() == sizeof(int)*3 );
		assert( t.extent() == sizeof(int)*3 );
		assert( t.lower_bound() == 0 );
		assert( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type t = mpi3::make_type<int>()[3];
		assert( t.size() == sizeof(int)*3 );
		assert( t.extent() == sizeof(int)*3 );
		assert( t.lower_bound() == 0 );
		assert( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type t = mpi3::make_type<double>()[3];
		assert( t.size() == sizeof(double)*3 );
		assert( t.extent() == sizeof(double)*3 );
	}
	{
		mpi3::type t = mpi3::make_type<std::complex<double>>()[3];
		assert( t.size() == sizeof(std::complex<double>)*3 );
		assert( t.extent() == sizeof(std::complex<double>)*3 );
	}
	{
		mpi3::type t = mpi3::float_int;
		assert( t.size() == sizeof(std::pair<float, int>) );
		struct foo_t {
			float a;
			int b;
		};
		[[maybe_unused]]  // nvcc false warning
		foo_t foo{};

		foo.a = 1.2F;
		foo.b = 5;

		assert( t.size() == sizeof(float) + sizeof(int) );
		assert( t.extent() == sizeof(foo) );
		assert( t.lower_bound() == 0 );
		assert( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		using T = std::complex<double>;
		T buffer[100];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

		switch(world.rank()) {
			break; case 0: {
				buffer[10] = 42.;
				MPI_Send(buffer, 1, &mpi3::make_type<T>()[100], 1, 0          , world.handle());
			}
			break; case 1: {
				MPI_Status status;
			    MPI_Recv(buffer, 1, &mpi3::make_type<T>()[100], 0, MPI_ANY_TAG, world.handle(), &status);
				assert( buffer[10] == 42. );
			}
		}
	}

	return 0;
} catch(...) {return 1;}

