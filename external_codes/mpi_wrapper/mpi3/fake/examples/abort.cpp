#include <mpi.h>
int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    MPI_Abort(MPI_COMM_WORLD, 911);
    /* No further code will execute */
	assert(0);
    MPI_Finalize();
    return 0;
}
