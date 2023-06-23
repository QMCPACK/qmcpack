#include "../../mpi3/main.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {

	int const data = [&] {switch(world.rank()) {
	    case 0: return 12;
	    case 1: return 34;
	    case 2: return 56;
	    case 3: return 78;  // default: world.abort(EXIT_FAILURE);
	} return 0;}();

	{
		auto reduction = world.max_loc(data);

		assert( reduction.value == 78 );
		assert( reduction.location == 3 );
	}
#if __cpp_structured_bindings >= 201606
	{
		auto [value, location] = world.max_loc(data);
		assert( value == 78 );
		assert( location == 3 );
	}
	{
		auto&& [value, procs] = world.max_location(data);
		assert( value == 78 );
		assert( procs.rank() == 3 );
	}
#endif
	{
		int const* max_ptr = world.max_element(data);
		if(world.rank() == 3) {
			assert( max_ptr and max_ptr == &data );
		} else {
			assert( not max_ptr );
		}
	}
	return 0;
} catch(...) {
	return 1;
}
