
// Compile-time test that all functions are appropriately declared 'inline' and
// will not give multiple definition errors

// For a simple check on multi-file compilation, use
// mpic++ library_check.cpp library_main.cpp

// To check when one file is in a library, use
// mpic++ -c library_check.cpp
// ar rcs liblibrary_check.a library_check.o
// mpic++ library_main.cpp -L. -llibrary_check

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;

// defined in library_check.cpp
void do_broadcast(mpi3::communicator &c);

int mpi3::main(int, char*[], mpi3::communicator world)
{
  int b = 1;
  world.broadcast_value(b);
  do_broadcast(world);

  return 0;
}
