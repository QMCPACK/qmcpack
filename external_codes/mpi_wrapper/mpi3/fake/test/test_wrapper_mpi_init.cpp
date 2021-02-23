
#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

namespace mpi3 = boost::mpi3;

int main(int argc, char **argv)
{
  mpi3::environment(argc, argv);
  return 0;
}
