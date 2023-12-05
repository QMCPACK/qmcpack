// Copyright 2023 Alfredo A. Correa

// based on http://mpi.deino.net/mpi_functions/MPI_Alltoall.html

#include <mpi3/main.hpp>
#if not defined(EXAMPI)
#include <mpi3/ostream.hpp>
#endif

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {
	std::size_t chunk = 5;

	auto sb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);
    std::iota(sb.begin(), sb.end(), 40000 + world.rank()*100);

	auto rb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);

	auto sz = static_cast<std::size_t>(world.size()); assert( sz != 0 );
	assert( sb.size() % sz == 0);

	world.all_to_all_n(sb.data(), sb.size()/sz, rb.data());

#if not defined(EXAMPI)
    mpi3::ostream wout(world);
    std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout<<"sb = ", ", ")); wout<<std::endl;
    std::copy(rb.begin(), rb.end(), std::ostream_iterator<int>(wout<<"rb = ", ", ")); wout<<std::endl;

	world.all_to_all_inplace_n(sb.data(), sb.size()/sz);
	// do_all_to_all_n(world, sb.data(), sb.size(), sb.data());
	// world.all_to_all_n(sb.data(), sb.size()); //  , sb.data());
	// std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout<<"sb (inplace) = ", ", ")); wout<<std::endl;
	assert(sb == rb);
#endif

	return 0;
} catch(...) {return 0;}
