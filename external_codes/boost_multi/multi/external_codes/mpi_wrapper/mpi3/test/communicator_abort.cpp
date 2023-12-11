#include "../main.hpp"
#include "../communicator.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int {
	world.abort(911);

	return 0;
}
