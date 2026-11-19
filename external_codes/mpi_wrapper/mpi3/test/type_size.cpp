// Copyright 2018-2025 Alfredo A. Correa

#define BOOST_MPI3_DISALLOW_AUTOMATIC_POD_COMMUNICATION

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/type.hpp>

#include <boost/core/lightweight_test.hpp>

#include <complex>
#include <utility>

#include <mpi3/detail/mpi_impl.h>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	{
		mpi3::type const& t = mpi3::int_;
		BOOST_TEST( t.size() == sizeof(int) );
		BOOST_TEST( t.extent() == sizeof(int) );
		BOOST_TEST( t.lower_bound() == 0 );
		BOOST_TEST( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type const t = mpi3::int_[3];
		BOOST_TEST( t.size() == sizeof(int)*3 );
		BOOST_TEST( t.extent() == sizeof(int)*3 );
		BOOST_TEST( t.lower_bound() == 0 );
		BOOST_TEST( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type const t = mpi3::make_type<int>()[3];
		BOOST_TEST( t.size() == sizeof(int)*3 );
		BOOST_TEST( t.extent() == sizeof(int)*3 );
		BOOST_TEST( t.lower_bound() == 0 );
		BOOST_TEST( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		mpi3::type const t = mpi3::make_type<double>()[3];
		BOOST_TEST( t.size() == sizeof(double)*3 );
		BOOST_TEST( t.extent() == sizeof(double)*3 );
	}
	{
		mpi3::type const t = mpi3::make_type<std::complex<double>>()[3];
		BOOST_TEST( t.size() == sizeof(std::complex<double>)*3 );
		BOOST_TEST( t.extent() == sizeof(std::complex<double>)*3 );
	}
	{
		mpi3::type const& t = mpi3::float_int;
		BOOST_TEST( t.size() == sizeof(std::pair<float, int>) );
		struct foo_t {
			float a;
			int b;
		};
		[[maybe_unused]]  // nvcc false warning
		foo_t foo{};

		foo.a = 1.2F;
		foo.b = 5;

		BOOST_TEST( t.size() == sizeof(float) + sizeof(int) );
		BOOST_TEST( t.extent() == sizeof(foo) );
		BOOST_TEST( t.lower_bound() == 0 );
		BOOST_TEST( t.upper_bound() == t.extent() - t.lower_bound() );
	}
	{
		using T = std::complex<double>;
		T buffer[100];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

		auto array_type = mpi3::make_type<T>()[100];  // named to avoid freeing the handle before MPI uses it

		switch(world.rank()) {
			case 0: {
				buffer[10] = 42.;
				MPI_Send(buffer, 1, &array_type, 1, 0          , world.handle());
			}
			break; case 1: {
				MPI_Status status;
			    MPI_Recv(buffer, 1, &array_type, 0, MPI_ANY_TAG, world.handle(), &status);
				BOOST_TEST( buffer[10] == 42. );
			}
			break; default: ;
		}
	}
	return boost::report_errors();
}
