#include "mpi.h"
#include <stdio.h>
 
int main(int argc, char *argv[])
{
//    int myrank;
//    MPI_Status status;
    MPI_Datatype type, type2;
//    int buffer[100];
 
    MPI_Init(&argc, &argv);
 
    MPI_Type_contiguous( 100, MPI_CHAR, &type );
    MPI_Type_commit(&type);
    MPI_Type_dup(type, &type2);
    MPI_Type_free(&type);
 
/*    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    if (myrank == 0)
    {
        MPI_Send(buffer, 1, type2, 1, 123, MPI_COMM_WORLD);
    }
    else if (myrank == 1)
    {
        MPI_Recv(buffer, 1, type2, 0, 123, MPI_COMM_WORLD, &status);
    }
*/
    MPI_Finalize();
    return 0;
}
