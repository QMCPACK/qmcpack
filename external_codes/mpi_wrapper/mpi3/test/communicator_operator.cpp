// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <mpi3/detail/datatype.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <iterator>
#include <numeric>
#include <thread>  // sleep_for
#include <vector>

#include <mpi3/detail/mpi_impl.h>

namespace mpi3 = boost::mpi3;

using std::cout;
using namespace std::chrono_literals;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	std::vector<double> inbuf(100);

	switch(world.rank()) {
		case 0: {
			std::vector<double> outbuf(100);
			iota(begin(outbuf), end(outbuf), 0.0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
			std::this_thread::sleep_for(2s);  // NOLINT(misc-include-cleaner) bug in clang-tidy 18.1
			cout <<"world["<< world.rank() <<"] about to isent\n";
			mpi3::request const r = world.isend(begin(outbuf), end(outbuf), 1);
			cout <<"comm["<< world.rank() <<"] isent\n";
		//  r.wait();
		}; break;
		case 1: {
			cout <<"comm["<< world.rank() <<"] about to ireceive\n";
			mpi3::request r;//= world.ireceive_n(inbuf.begin(), inbuf.size(), 0);
			MPI_Irecv(
				inbuf.data(), static_cast<int>(inbuf.size()),
				mpi3::detail::basic_datatype<double>{},
				MPI_ANY_SOURCE, MPI_ANY_TAG, world.get(), &r.impl_
			);
			cout <<"comm["<< world.rank() <<"] ireceived\n";
			// cppcheck-suppress[cstyleCast,intToPointerCast] 
			MPI_Wait(&r.impl_, MPI_STATUS_IGNORE);  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast) for macro
			//  r.wait();
		}; break;
		default: break;
	}
	cout <<"comm["<< world.rank() <<"] completed op\n";

	if(world.rank() == 1) {BOOST_TEST( inbuf[9] == 9. );}

	return boost::report_errors();
} catch(...) {
	return 1;
}
