#include "mpi.h"
#include <stdio.h>
 
int main(int argc, char *argv[])
{
    int myrank;
    MPI_Status status;
    MPI_Datatype type;
    int buffer[100];
 
    MPI_Init(&argc, &argv);
 
    MPI_Type_contiguous( 100, MPI_CHAR, &type );
    MPI_Type_commit(&type);
 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    if (myrank == 0)
    {
        MPI_Send(buffer, 1, type, 1, 123, MPI_COMM_WORLD);
    }
    else if (myrank == 1)
    {
        MPI_Recv(buffer, 1, type, 0, 123, MPI_COMM_WORLD, &status);
    }
 
    MPI_Finalize();
    return 0;
}
