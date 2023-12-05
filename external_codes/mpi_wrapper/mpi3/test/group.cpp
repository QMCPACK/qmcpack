/*-*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-*/
// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/environment.hpp"
#include "../../mpi3/group.hpp"

#include "../../mpi3/main.hpp"

#include <iostream>

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {
	mpi3::group wg{world};
	mpi3::communicator w2{wg};

	assert(w2.rank() == world.rank());
	assert(w2.size() == world.size());

	mpi3::communicator half = world / 2;
	mpi3::group hg{half};

	mpi3::communicator h2{hg};
	assert(half.rank() == h2.rank());

	static_assert(std::is_same<decltype(&wg), MPI_Group>{}, "!");

	[[maybe_unused]] mpi3::group const& wgc = wg;
	static_assert(std::is_same<decltype(&wgc), mpi3::group const*>{}, "!");

	// static_assert( std::is_same<decltype(*&wg), mpi3::group&>{}, "!" );

	return 0;
} catch(...) {
	return 1;
}
