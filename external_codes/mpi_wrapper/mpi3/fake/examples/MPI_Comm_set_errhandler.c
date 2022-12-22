#include "mpi.h"
#include <stdio.h>
#include<assert.h>
static int calls = 0;
static int errs = 0;
static MPI_Comm mycomm;

void eh( MPI_Comm *comm, int *err, ... )
{
    printf( "eh called\n" );fflush(stdout);
    if (*err != MPI_ERR_OTHER) {
        errs++;
        printf( "Unexpected error code\n" );fflush(stdout);
    }
    if (*comm != mycomm) {
        errs++;
        printf( "Unexpected communicator\n" );fflush(stdout);
    }
    calls++;
    return;
}

int main( int argc, char *argv[] )
{
    MPI_Comm comm;
    MPI_Errhandler newerr;

    MPI_Init( &argc, &argv );
    comm = MPI_COMM_WORLD;
    mycomm = comm;
    MPI_Comm_create_errhandler( eh, &newerr );
    int s = -1;
    s = MPI_Comm_set_errhandler( comm, newerr ); // fatal error in OpenMPI
    assert(s == MPI_SUCCESS);
    int rank = -1;
    MPI_Comm_rank( comm, &rank ); // error but handler not called
    MPI_Comm_call_errhandler( comm, MPI_ERR_OTHER ); // ok, it seems
    MPI_Errhandler_free( &newerr );
    if (calls != 1) {
        errs++;
        printf( "Error handler not called\n" );fflush(stdout);
    }
    MPI_Finalize();
    return 0;
}
