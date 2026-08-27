// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/window.hpp>

#include <boost/core/lightweight_test.hpp>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	mpi3::communicator const comm = (world < 2);
	(void)comm;
	// #if 0
	if(comm) {
		// std::vector<double> inbuf(100);
		// std::vector<double> outbuf(100);

		// std::iota(outbuf.begin(), outbuf.end(), 0.0);
		// mpi3::window<double> win{
		// 	comm, comm.rank()==1?inbuf.data():nullptr,
		// 	comm.rank()==1?inbuf.size():0
		// };
		// auto win = comm.make_window<double>(comm.rank()==1?inbuf.size():0);
		// comm.make_window<double>(comm.rank()==1?inbuf.data():nullptr, comm.rank()==1?inbuf.size():0);
		// comm.make_window<double>(comm.rank()==1?inbuf.data():nullptr, comm.rank()==1?inbuf.size():0);

		// win.fence();
		// if(world.rank() == 0) win.put_n(outbuf.data(), outbuf.size(), 1);
		// win.fence();

		// std::cout << "rank " << comm.rank() << " " << inbuf[7] << '\n';
		// if(world.rank() == 1) BOOST_TEST( inbuf[7] == 7.0 );
	}
	// #endif

	return boost::report_errors();
} catch(...) {
	return 1;
}
