// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <algorithm>
#include <boost/core/lightweight_test.hpp>

#include <cstddef>
#include <future>
#include <iostream>
#include <iterator>
#include <vector>

#include <mpi3/detail/mpi_impl.h>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	// static_assert(sizeof(MPI_Comm) == sizeof(mpi3::communicator));

	{
		mpi3::communicator const& w2 = mpi3::grip_communicator(world.handle());
		BOOST_TEST(  w2 ==  world );
		BOOST_TEST( &world == &w2 );

		BOOST_TEST(  mpi3::grip_communicator(MPI_COMM_WORLD) ==  world );
		BOOST_TEST( &mpi3::grip_communicator(MPI_COMM_WORLD) != &world );
	}

	//  BOOST_TEST( reinterpret_cast<mpi3::communicator&>(MPI_COMM_WORLD) == world );

	{
		mpi3::communicator mty;
		BOOST_TEST( mty.empty() );
		//  BOOST_TEST( mty.size() == 0 );
		mpi3::communicator const mty2 = mty;
		BOOST_TEST( mty2.empty() );
		//  BOOST_TEST( mty2.size() == 0 );
	}
	{
		std::vector<mpi3::communicator> comms;
		comms.emplace_back(world);
		comms.emplace_back(world);
		comms.emplace_back(world);

		std::vector<mpi3::communicator> comms2;
		comms2.reserve(3);
		//  for(auto& e:comms) {comms2.emplace_back(e);}  // ok, but old style
		//  std::copy(comms.begin(), comms.end(), std::back_inserter(comms2));  // doesn't work because it calls cctor
		std::transform(comms.begin(), comms.end(), std::back_inserter(comms2), [](auto&& e) { return e; });  // calls dup ctor    // NOLINT(boost-use-ranges) for C++20, use std::ranges::transform
	}
	{
		std::size_t const             NTHREADS = 10;
		std::vector<std::future<int>> fs;
		for(int i = 0; i != NTHREADS; ++i) {  // NOLINT(altera-unroll-loops)
			// #if 0 // this is problematic because copy (mpi_comm_dup) is not thread safe
			//      fs.emplace_back(std::async([&world](){
			//          auto comm = world; // hangs here
			//          return 5;
			//      }));
			// #endif
			fs.emplace_back(std::async([/*thread_comm = world*/]() {
				//  auto comm2 = thread_comm; // works, just to test
				return 5;
			}));
// #if 0  // more correct
//      fs.emplace_back(std::async([world](){
//          auto comm2 = world; // works, just to test
//          return 5;
//      }));
// #endif
			std::cout << "created thread\n";
		}
		for(std::size_t i = 0; i != NTHREADS; ++i) {  // NOLINT(altera-unroll-loops)
			auto five = fs[i].get();
			BOOST_TEST( five == 5 );
			std::cout << "joined thread\n";
		}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
