#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <mpi.h>



int main(int argc, char **argv){ 
    // argc is argument count and argv is a list of arguments..
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Finalize();
    exit(0);
}