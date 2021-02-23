
#include <mpi.h>

// Test MPI initialization in C++

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Finalize();
  return 0;
}
