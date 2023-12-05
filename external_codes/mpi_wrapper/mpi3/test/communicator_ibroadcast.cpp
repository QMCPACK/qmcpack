#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>
#if not defined(EXAMPI)
#include <mpi3/ostream.hpp>
#endif

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	assert(world.size() > 2);

	std::vector<int> large(10);
	if(world.root()) {
		iota(large.begin(), large.end(), 0);
	}

#if not defined(EXAMPI)
	mpi3::ostream wout(world);
	wout << "before:" << std::endl;
	std::copy(large.begin(), large.end(), std::ostream_iterator<int>(wout, " "));

	wout << std::endl;

	{
		auto req = world.ibroadcast(large.begin(), large.end(), 0);
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(5s);
	}

	wout << "after:" << std::endl;
	std::copy(large.begin(), large.end(), std::ostream_iterator<int>(wout, " "));
	wout << std::endl;
#endif

	return 0;
} catch(...) {return 1;}
