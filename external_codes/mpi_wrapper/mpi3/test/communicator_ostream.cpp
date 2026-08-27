// Copyright 2022-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/ostream.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <limits>
#include <random>

namespace mpi3 = boost::mpi3;


auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST( world.size() > 2 );

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis{-1.0, +1.0};

	world.set_name("world");

	mpi3::ostream wout{world, std::cout};

//  wout << mpi3::set_communicator_logging;  // TODO

	wout << "Program starts\n" << std::flush;
	wout << "Hello! for world using "<< world.size() <<" processes\n" << std::flush;

	wout << "Hello! I am rank "<< world.rank()<< " in " << world.name() << '\n' << std::flush;

	wout << (world.root()?"this precess is root\n":"this process is NOT root\n") << std::flush;

	wout << "rank "         << world.rank()                                << '\t' << std::flush;
	wout << "small random " << dis(gen)                                    << '\t' << std::flush;
	wout << "large random " << dis(gen)*std::numeric_limits<double>::max() << '\t' << std::flush;

	wout << "-------------------\n" << std::flush;

	wout << "raw_stuff " << world.rank() << std::flush;
	wout << "\nsomething random"         << std::flush;

	wout << "Program Ends\n" << std::flush;

	if(mpi3::communicator firsttwo = (world < 2) ) {
		firsttwo.set_name("firsttwo");
		mpi3::ostream fout(firsttwo);
		fout
			<<"Hola! I am rank "<< firsttwo.rank() <<" in "<< firsttwo.name()
			<<" and also rank "<< world.rank() <<" in "<< world.name()
			<<'\n'
		;
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
