#include "mpi.h"
#include <stdio.h>

int main( int argc, char *argv[] )
{
    int errs = 0;
    int size, dims[2], periods[2], remain[2];
    int result;
    MPI_Comm comm, newcomm;

    MPI_Init( &argc, &argv );

    /* First, create a 1-dim cartesian communicator */
    periods[0] = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    dims[0] = size;
    MPI_Cart_create( MPI_COMM_WORLD, 1, dims, periods, 0, &comm );

    /* Now, extract a communicator with no dimensions */
    remain[0] = 0;
    MPI_Cart_sub( comm, remain, &newcomm );

    /* This should be congruent to MPI_COMM_SELF */
    MPI_Comm_compare( MPI_COMM_SELF, newcomm, &result );
    if (result != MPI_CONGRUENT) {
        errs++;
        printf( "cart sub to size 0 did not give self\n" );fflush(stdout);
    }

    /* Free the new communicator */
    MPI_Comm_free( &newcomm );
    MPI_Comm_free( &comm );

    MPI_Finalize();
    return 0;
}
