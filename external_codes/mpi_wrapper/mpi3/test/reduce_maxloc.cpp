// Copyright 2020-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>  // NOLINT(misc-include-cleaner) used indirectly

#include <boost/core/lightweight_test.hpp>

#include <cstdlib>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	int const data = [&] {
		switch(world.rank()) {
		case 0 : return 12;
		case 1 : return 34;
		case 2 : return 56;
		case 3 : return 78;
		default: world.abort(EXIT_FAILURE);
		}
		return 0;
	}();

	{
		auto reduction = world.max_loc(data);

		BOOST_TEST(reduction.value == 78);
		BOOST_TEST(reduction.location == 3);
	}
#if __cpp_structured_bindings >= 201606
	{
		auto [value, location] = world.max_loc(data);
		BOOST_TEST(value == 78);
		BOOST_TEST(location == 3);
	}
	{
		auto&& [value, procs] = world.max_location(data);
		BOOST_TEST(value == 78);
		BOOST_TEST(procs.rank() == 3);
	}
#endif
	{
		int const* max_ptr = world.max_element(data);
		if(world.rank() == 3) {
			BOOST_TEST(max_ptr and max_ptr == &data);
		} else {
			BOOST_TEST(not max_ptr);
		}
	}
	return boost::report_errors();
}
