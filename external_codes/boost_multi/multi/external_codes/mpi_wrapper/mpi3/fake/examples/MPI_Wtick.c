#include "mpi.h"
#include <stdio.h>

int main( int argc, char *argv[] )
{
    double tick;

    MPI_Init( 0, 0 );
    tick = MPI_Wtick();
    printf("A single MPI tick is %0.9f seconds\n", tick);
    fflush(stdout);
    MPI_Finalize( );
    return 0;
} 
