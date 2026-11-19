// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <mpi3/detail/mpi_impl.h>
#include <mpi3/group.hpp>

#include <boost/core/lightweight_test.hpp>
#include <type_traits>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	mpi3::group const        wg{world};
	mpi3::communicator const w2{wg};

	BOOST_TEST(w2.rank() == world.rank());
	BOOST_TEST(w2.size() == world.size());

	mpi3::communicator const half = world / 2;
	mpi3::group const        hg{half};

	mpi3::communicator const h2{hg};
	BOOST_TEST(half.rank() == h2.rank());

	static_assert(std::is_same_v<decltype(&const_cast<mpi3::group&>(wg)), MPI_Group>);  // NOLINT(cppcoreguidelines-pro-type-const-cast) test const_cast

	[[maybe_unused]] mpi3::group const& wgc = wg;
	static_assert(std::is_same_v<decltype(&wgc), mpi3::group const*>);

	// static_assert( std::is_same<decltype(*&wg), mpi3::group&>{}, "!" );

	return boost::report_errors();
} catch(...) {
	return 1;
}
