#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    int buffer[6];
    int rank, size, i;
    int receive_counts[4] = { 0, 1, 2, 3 };
    int receive_displacements[4] = { 0, 0, 1, 3 };

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size != 4)
    {
        if (rank == 0)
        {
            printf("Please run with 4 processes\n");fflush(stdout);
        }
        MPI_Finalize();
        return 0;
    }
    for (i=0; i<rank; i++)
    {
        buffer[i] = rank;
    }
    MPI_Gatherv(buffer, rank, MPI_INT, buffer, receive_counts, receive_displacements, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        for (i=0; i<6; i++)
        {
            printf("[%d]", buffer[i]);
        }
        printf("\n");
        fflush(stdout);
    }
    MPI_Finalize();
    return 0;
}
