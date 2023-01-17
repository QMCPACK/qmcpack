// Copyright 2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"
#include "../../mpi3/ostream.hpp"

#include <fstream>
#include <random>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argv*/, char** /*argc*/, mpi3::communicator world) -> int try {
	assert( world.size() > 2 );

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis{-1.0, +1.0};

	world.set_name("world");

	mpi3::ostream wout{world, std::cout};

//  wout << mpi3::set_communicator_logging;  // TODO

	wout << "Program starts" << std::endl;
	wout << "Hello! for world using "<< world.size() <<" processes" << std::endl;

	wout << "Hello! I am rank "<< world.rank()<< " in " << world.name() << std::endl;

	wout << (world.root()?"this precess is root":"this process is NOT root") << std::endl;

	wout << "rank "         << world.rank()                                << '\t' << std::flush;
	wout << "small random " << dis(gen)                                    << '\t' << std::flush;
	wout << "large random " << dis(gen)*std::numeric_limits<double>::max() << '\t' << std::flush;

	wout << "-------------------" << std::endl;

	wout << "raw_stuff " << world.rank() << std::flush;
	wout << "\nsomething random"         << std::flush;

	wout << "Program Ends" << std::endl;

	if(mpi3::communicator firsttwo = (world < 2) ) {
		firsttwo.set_name("firsttwo");
		mpi3::ostream fout(firsttwo);
		fout
			<<"Hola! I am rank "<< firsttwo.rank() <<" in "<< firsttwo.name()
			<<" and also rank "<< world.rank() <<" in "<< world.name()
			<<std::endl
		;
	}

	return 0;
} catch(...) {return 1;}
