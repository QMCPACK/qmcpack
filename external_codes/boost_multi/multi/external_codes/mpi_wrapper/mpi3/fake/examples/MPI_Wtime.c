#include "mpi.h"
//#include <windows.h>
#include <unistd.h>
#include <stdio.h>

int main( int argc, char *argv[] )
{
    double t1, t2;

    MPI_Init( 0, 0 );
    t1 = MPI_Wtime();
    sleep(5);
    t2 = MPI_Wtime();
    printf("MPI_Wtime measured a 1 second sleep to be: %1.2f\n", t2-t1);fflush(stdout);
    MPI_Finalize( );
    return 0;
} 
