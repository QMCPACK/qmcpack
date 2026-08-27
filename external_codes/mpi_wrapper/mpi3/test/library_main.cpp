//#if COMPILATION_INSTRUCTIONS
//mpic++ library_check.cpp library_main.cpp
//#mpic++ -c library_check.cpp
//#ar rcs liblibrary_check.a library_check.o
//#mpic++ library_main.cpp -o library_main.x -L. -llibrary_check
//mpirun -n 4 ./library_main.x;exit
//#endif
// Compile-time test that all functions are appropriately declared 'inline' and
// will not give multiple definition errors

// Copyright 2018-2025 Alfredo A. Correa

// For a simple check on multi-file compilation, use
// mpic++ library_check.cpp library_main.cpp

// To check when one file is in a library, use
// mpic++ -c library_check.cpp
// ar rcs liblibrary_check.a library_check.o
// mpic++ library_main.cpp -L. -llibrary_check

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

namespace mpi3 = boost::mpi3;

void do_broadcast(mpi3::communicator& c);  // NOLINT(misc-use-internal-linkage)  // defined in library_check.cpp

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	int b = 1;
	world.broadcast_value(b);
	do_broadcast(world);
}
