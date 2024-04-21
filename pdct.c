#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>
#include <mpi.h>


/*
A program to compute the parallel DCT and IDCT using the FFT and IFFT

The program takes (at the moment) and extra agrument (more than standard and MPI arguments)
this argument is the exponent of two which defines the length of the list
 
 TODO:: 
 1. implement scheme for DFT and invDFT 
 2. parallelize embarrasingly parallel post processing for the inverse transform. 
 3. figure out scheme for pre processing DCT (and implement)
 4. figure out scheme for pre processing inverse DCT (and implement)

 (
    5. start with multidim schemes 
    6. try algorithm on picture files (black and white)
    7. Possibly try scheme on color picture as just three separate schemes
 )
 */


void pfft(complex double *vec, int len, int localLen, int rank, int size){
    int tag1 = 1, tag2 = 2;
    printf("got here 3");
    complex double recvd, dummy;
    for (int i = log2(len)-1; i >= 0; i--){
        int fac = pow(2, i);
        for (int j = 0; j < localLen; j++){
            printf("got here 4");
            if (fac > localLen){ // communicate
                if (j % 2*fac < fac) { // forward communication 
                    MPI_Send(&vec[j], 1, MPI_C_DOUBLE_COMPLEX, rank+fac, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd, 1, MPI_C_DOUBLE_COMPLEX, rank + fac, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    vec[j] = vec[j] + cexp(-2.0*I*M_PI*j/len) * recvd;
                } else { // backwards communication
                    MPI_Recv(&recvd, 1, MPI_C_DOUBLE_COMPLEX, rank - fac, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Send(&vec[j], 1, MPI_C_DOUBLE_COMPLEX, rank- fac, tag2, MPI_COMM_WORLD);
                    vec[j] = recvd - cexp(-2.0*I*M_PI*(j-fac)/len) * vec[j];
                }
            } else{//no communication
                dummy = vec[j] + cexp(-2.0*I*M_PI*(j)/len) * vec[j+ fac]; 
                vec[j+ fac] = vec[j] - cexp(-2.0*I*M_PI*(j)/len) * vec[j+ fac];
                vec[j] = dummy;
            }
        }
    }
}


void pIfft(complex double *vec, int len, int localLen, int rank, int size){

}




int main(int argc, char **argv){ 
    // argc is argument count and argv is a list of arguments..
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if (argc < 2){
        printf("Not enough input arguments");
        exit(0);
    }
    int N = pow(2, atoi(argv[2]));
    srand(time(NULL));
    printf("got here 1");
    int J = (rank < N%size) ? N/size + 1: N/size;
    double complex* randomVec = malloc(J*sizeof(double complex));
    for (int i=0; i < J; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
    }
    printf("got here 2");
    pfft(randomVec, N, J, rank, size);
    free(randomVec);
    MPI_Finalize();
    exit(0); 
}