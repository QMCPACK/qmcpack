// Copyright 2021-2025 Alfredo A. Correa

#include <mpi3/cartesian_communicator.hpp>
#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <mpi3/detail/mpi_impl.h>

#include <boost/core/lightweight_test.hpp>

#include <iostream>
#include <stdexcept>
#include <tuple>

namespace mpi3 = boost::mpi3;

namespace {
void division_tests1() {
	{
		auto div = mpi3::cartesian_communicator<2>::division(6);
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 3);
		BOOST_TEST(div[1] == 2);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 3);
		BOOST_TEST(div[1] == 2);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 3);
		BOOST_TEST(div[1] == 2);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill, mpi3::fill});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 3);
		BOOST_TEST(div[1] == 2);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {2});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 3);
	}
}

void division_tests2() {
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {2, mpi3::fill});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 3);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::fill, 3});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 3);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(6, {mpi3::_, 3});
		BOOST_TEST(div[0] * div[1] == 6);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 3);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(7);
		BOOST_TEST(div[0] * div[1] == 7);
		BOOST_TEST(div[0] == 7);
		BOOST_TEST(div[1] == 1);
	}
	{
		auto div = mpi3::cartesian_communicator<2>::division(7, {mpi3::fill, mpi3::fill});
		BOOST_TEST(div[0] * div[1] == 7);
		BOOST_TEST(div[0] == 7);
		BOOST_TEST(div[1] == 1);
	}

	try {  // this is an error in MPICH and openMPI
		auto const div = mpi3::cartesian_communicator<2>::division(7, {2, mpi3::fill});
		BOOST_TEST(div[0] * div[1] == 4);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 2);
	} catch(std::runtime_error&) {}  // NOLINT(bugprone-empty-catch)

	try {  // this is an error in MPICH
		auto const div = mpi3::cartesian_communicator<2>::division(6, {2, 2});
		BOOST_TEST(div[0] * div[1] == 4);
		BOOST_TEST(div[0] == 2);
		BOOST_TEST(div[1] == 2);
	} catch(std::runtime_error&) {  // NOLINT(bugprone-empty-catch)
	}
}
} // end namespace

auto main(int argc, char** argv) -> int try {  // NOLINT(readability-function-cognitive-complexity)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() == 6);

	division_tests1();
	division_tests2();

	{
		mpi3::cartesian_communicator<2> cart_comm(world, {3, 2});
		BOOST_TEST(cart_comm.dimensions()[0] == 3);
		BOOST_TEST(cart_comm.dimensions()[1] == 2);

		auto row = cart_comm.axis(0);
		auto col = cart_comm.axis(1);
		BOOST_TEST(row.size() == 3);
		BOOST_TEST(col.size() == 2);
	}
	{
		mpi3::cartesian_communicator<2> cart_comm(world, {mpi3::fill, 2});
		BOOST_TEST(cart_comm.dimensions()[0] == 3);
		BOOST_TEST(cart_comm.dimensions()[1] == 2);

		auto row = cart_comm.axis(0);
		auto col = cart_comm.axis(1);
		BOOST_TEST(row.size() == 3);
		BOOST_TEST(col.size() == 2);
	}
	{
		mpi3::cartesian_communicator<2> cart_comm(world, {3, mpi3::fill});
		BOOST_TEST(cart_comm.dimensions()[0] == 3);
		BOOST_TEST(cart_comm.dimensions()[1] == 2);

		auto row = cart_comm.axis(0);
		auto col = cart_comm.axis(1);
		BOOST_TEST(row.size() == 3);
		BOOST_TEST(col.size() == 2);
	}

	{
		mpi3::cartesian_communicator<2> cart_comm(world, {3, 2});
		BOOST_TEST(cart_comm.dimensions()[0] == 3);
		BOOST_TEST(cart_comm.dimensions()[1] == 2);

		auto row = cart_comm.axis(0);
		auto col = cart_comm.axis(1);
		BOOST_TEST(row.size() == 3);
		BOOST_TEST(col.size() == 2);

		{
			auto comm_sub0 = cart_comm.axis(0);
			BOOST_TEST(comm_sub0.shape()[0] == 3);
			BOOST_TEST(comm_sub0.size() == 3);
		}
		{
			auto comm_sub1 = cart_comm.axis(1);
			BOOST_TEST(comm_sub1.shape()[0] == 2);
			BOOST_TEST(comm_sub1.size() == 2);
		}
		{
			auto plane0 = cart_comm.hyperplane(0);
			static_assert(decltype(plane0)::dimensionality == 1, "!");
			BOOST_TEST(plane0.size() == 2);
		}
		{
			auto plane1 = cart_comm.hyperplane(1);
			static_assert(decltype(plane1)::dimensionality == 1, "!");
#if !defined(MPICH_VERSION) && !defined(MSMPI_VER)
			{  // cartesian communicators of dimension 0 do not work in MPICH or MS-MPI
				BOOST_TEST(plane1.size() == 3);
				auto point = plane1.hyperplane(0);
				BOOST_TEST(point.num_elements() == 1);
			}
#endif
		}
		{
			auto comm_sub0 = cart_comm.axis<0>();
			BOOST_TEST(comm_sub0.shape()[0] == 3);
			BOOST_TEST(comm_sub0.size() == 3);
		}
		{
			auto plane = cart_comm.plane<0, 1>();
			BOOST_TEST(plane.shape() == cart_comm.shape());
		}
	}
	{
		mpi3::cartesian_communicator<3> cart_comm3d(world, {3, 1, 2});
		BOOST_TEST(cart_comm3d.dimensions()[0] == 3);
		BOOST_TEST(cart_comm3d.dimensions()[1] == 1);
		BOOST_TEST(cart_comm3d.dimensions()[2] == 2);

		auto cart_comm2d(cart_comm3d.flatten());
		BOOST_TEST(cart_comm2d.dimensions()[0] == 3);
		BOOST_TEST(cart_comm2d.dimensions()[1] == 2);

		auto cart_comm1d(cart_comm2d.flatten());
		BOOST_TEST(cart_comm1d.dimensions()[0] == 6);

		BOOST_TEST( world.compare(cart_comm1d) == mpi3::detail::equality::congruent );

		auto cart_comm2d_last(cart_comm3d.flatten_last());
		BOOST_TEST(cart_comm2d_last.dimensions()[0] == 3);
		BOOST_TEST(cart_comm2d_last.dimensions()[1] == 2);

		auto cart_comm1d_last(cart_comm2d.flatten_last());
		BOOST_TEST(cart_comm1d_last.dimensions()[0] == 6);

		BOOST_TEST( world.compare(cart_comm1d_last) == mpi3::detail::equality::congruent );
	}
	{
		mpi3::cartesian_communicator<3> cart_comm(world, {3, 2, 1});
		BOOST_TEST(cart_comm.dimensions()[0] == 3);
		BOOST_TEST(cart_comm.dimensions()[1] == 2);
		BOOST_TEST(cart_comm.dimensions()[2] == 1);

		{
			mpi3::cartesian_communicator<2> const plane = cart_comm.plane<0, 1>();
			BOOST_TEST(plane.shape()[0] == cart_comm.shape()[0]);
			BOOST_TEST(plane.shape()[1] == cart_comm.shape()[1]);
		}
		{
			auto plane = cart_comm.plane<0, 1>();
			BOOST_TEST(plane.shape()[0] == cart_comm.shape()[0]);
			BOOST_TEST(plane.shape()[1] == cart_comm.shape()[1]);
		}
		{
			mpi3::cartesian_communicator<2> const plane = cart_comm.plane<0, 2>();
			BOOST_TEST(plane.shape()[0] == cart_comm.shape()[0]);
			BOOST_TEST(plane.shape()[1] == cart_comm.shape()[2]);
		}
		{
			auto plane_comm = cart_comm.plane(0, 1);
			BOOST_TEST(plane_comm.shape()[0] == cart_comm.shape()[0]);
			BOOST_TEST(plane_comm.shape()[1] == cart_comm.shape()[1]);
		}
		{
			auto plane_comm = cart_comm.plane();
			BOOST_TEST(plane_comm.shape()[0] == cart_comm.shape()[0]);
			BOOST_TEST(plane_comm.shape()[1] == cart_comm.shape()[1]);
		}
	}
	{
		mpi3::cartesian_communicator<2> cart_comm(world, {3, 2});

		BOOST_TEST(cart_comm.rank() == cart_comm.rank(cart_comm.coordinates()));
		BOOST_TEST(cart_comm.coordinates() == cart_comm.coordinates(cart_comm.rank()));

		BOOST_TEST(cart_comm(2, 1).rank() == cart_comm.rank({2, 1}));
		BOOST_TEST(std::apply(cart_comm, cart_comm.coordinates()).rank() == cart_comm.rank());

		BOOST_TEST(cart_comm.rank(cart_comm.coordinates(4)) == 4);

		std::cout << cart_comm.rank({5, 3}) << '\n';

		//  switch(world.rank()) {
		//      case 1: std::cout<< world.rank() <<" "<< cart_comm.coordinates()[0] <<", "<< cart_comm.coordinates()[1] <<std::endl;
		//      case 5: std::cout<< world.rank() <<" "<< cart_comm.coordinates()[0] <<", "<< cart_comm.coordinates()[1] <<std::endl;
		//  }
	}
	{
		mpi3::cartesian_communicator<1> const comm_1D(world, {world.size()});
		BOOST_TEST(comm_1D.periods()[0] == true);
		BOOST_TEST(comm_1D.rank({-1}) == comm_1D.size() - 1);

		auto next_rank = comm_1D.rank({comm_1D.coordinates()[0] + 1});
		auto prev_rank = comm_1D.rank({comm_1D.coordinates()[0] - 1});

		BOOST_TEST(comm_1D.rank() != next_rank);
		BOOST_TEST(comm_1D.rank() != prev_rank);
	}
	{
		auto const periodic = mpi3::cartesian_communicator<1>{world};  // implivit , {}, {true}};
		{
			auto const [source, dest] = periodic.shift<0>(+1);
			auto const periodic_size = periodic.size();
			if(periodic_size != 0) {
				BOOST_TEST(source == ((periodic.rank() - 1 + periodic_size) % periodic_size));
				BOOST_TEST(dest == ((periodic.rank() + 1) % periodic_size));
			}
		}
		{
			auto const [source, dest] = periodic.shift<0>(0);
			BOOST_TEST(source == periodic.rank());
			BOOST_TEST(dest == periodic.rank());
		}
		{
			auto const [source, dest] = periodic.shift<0>(-1);
			auto const periodic_size = periodic.size();
			if(periodic_size != 0) {
				BOOST_TEST(source == ((periodic.rank() + 1) % periodic_size));
				BOOST_TEST(dest == ((periodic.rank() - 1 + periodic_size) % periodic_size));
			}
		}
	}
	{
		auto const nonperiodic = mpi3::cartesian_communicator<1>{world, {}, {false}};
		{
			auto const [source, dest] = nonperiodic.shift<0>(+1);

			if(nonperiodic.rank() > 0) { BOOST_TEST(source == nonperiodic.rank() - 1); }
			if(nonperiodic.rank() == 0) { BOOST_TEST(source == MPI_PROC_NULL); }

			if(nonperiodic.rank() < nonperiodic.size() - 1) { BOOST_TEST(dest == nonperiodic.rank() + 1); }
			if(nonperiodic.rank() == nonperiodic.size() - 1) { BOOST_TEST(dest == MPI_PROC_NULL); }
		}
		{
			auto const [source, dest] = nonperiodic.shift<0>(0);
			BOOST_TEST(source == nonperiodic.rank());
			BOOST_TEST(dest == nonperiodic.rank());
		}
		{
			auto const [source, dest] = nonperiodic.shift<0>(-1);
			if(nonperiodic.rank() < nonperiodic.size() - 1) { BOOST_TEST(source == nonperiodic.rank() + 1); }
			if(nonperiodic.rank() == nonperiodic.size() - 1) { BOOST_TEST(source == MPI_PROC_NULL); }

			if(nonperiodic.rank() > 0) { BOOST_TEST(dest == nonperiodic.rank() - 1); }
			if(nonperiodic.rank() == 0) { BOOST_TEST(dest == MPI_PROC_NULL); }
		}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
