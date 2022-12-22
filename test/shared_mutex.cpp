// © Alfredo Correa 2018-2020

#include "../../mpi3/main.hpp"
#include "../../mpi3/environment.hpp"
#include "../../mpi3/shm/mutex.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator/*world*/) -> int {
	return 0;
}
