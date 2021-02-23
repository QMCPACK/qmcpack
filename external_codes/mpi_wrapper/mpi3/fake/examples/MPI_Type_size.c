#include "mpi.h"
#include <stdio.h>
#include <string.h>
 
int main(int argc, char **argv)
{
    int /*mpi_err,*/ errs = 0, size;
    MPI_Aint lb, ub, extent;
    MPI_Datatype type;
    struct { float a; int b; } foo;

    MPI_Init(&argc, &argv);
 
    type = MPI_INT;
    MPI_Type_size(type, &size);
    if (size != sizeof(int)) {
        fprintf(stderr, "MPI_Type_size of MPI_INT incorrect size (%d); should be %d.\n", size, (int) sizeof(int));fflush(stderr);
        errs++;
    }
 
    MPI_Type_get_extent(type, &lb, &extent);
    if (extent != sizeof(int)) {
        fprintf(stderr, "MPI_Type_get_extent of MPI_INT returned incorrect extent (%d); should be %d.\n", (int) extent, (int) sizeof(int));fflush(stderr);
        errs++;
    }
    if (lb != 0) {
        fprintf(stderr, "MPI_Type_get_extent of MPI_INT returned incorrect lb (%d); should be 0.\n", (int) lb);fflush(stderr);
        errs++;
    }
 
    MPI_Type_ub(type, &ub);
    if (ub != extent - lb) {
        fprintf(stderr, "MPI_Type_ub of MPI_INT returned incorrect ub (%d); should be %d.\n", (int) ub, (int) (extent - lb));fflush(stderr);
        errs++;
    }
 
    type = MPI_FLOAT_INT;
    MPI_Type_size(type, &size);
    if (size != sizeof(float) + sizeof(int)) {
        fprintf(stderr, "MPI_Type_size of MPI_FLOAT_INT returned incorrect size (%d); should be %d.\n", size, (int) (sizeof(float) + sizeof(int)));fflush(stderr);
        errs++;
    }
 
    MPI_Type_get_extent(type, &lb, &extent);
    if (extent != sizeof(foo)) {
        fprintf(stderr, "MPI_Type_get_extent of MPI_FLOAT_INT returned incorrect extent (%d); should be %d.\n", (int) extent, (int) sizeof(foo));fflush(stderr);
        errs++;
    }
    if (lb != 0) {
        fprintf(stderr, "MPI_Type_get_extent of MPI_FLOAT_INT returned incorrect lb (%d); should be 0.\n", (int) lb);fflush(stderr);
        errs++;
    }
 
    MPI_Type_ub(type, &ub);
    if (ub != extent - lb) {
        fprintf(stderr, "MPI_Type_ub of MPI_FLOAT_INT returned incorrect ub (%d); should be %d.\n", (int) ub, (int) (extent - lb));fflush(stderr);
        errs++;
    }
 
    MPI_Finalize();
    return errs;
}
