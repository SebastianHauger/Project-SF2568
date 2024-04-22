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

void arrayPlus(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int i = 0; i < localLength; i++){
        vec1[i] += cexp(-(2.0*I*M_PI*i)/length) * vec2[i];
    }
}


void arrayMinus(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int i=0; i < localLength; i++){
        vec1[i] = vec2[i] - cexp(-(2.0*I*M_PI*i)/length) * vec1[i];
    }
}


void pFft(complex double *vec, int len, int localLen, int rank, int size){
    int tag1 = 1, tag2 = 2;
    complex double dummy;
    double complex* recvd = malloc(localLen*sizeof(complex double));
    for (int i = 0; i < log2(len); i++){
        int fac = pow(2, i);
        if (fac >= localLen){ // communicate
            if ((rank*localLen) % (2*fac) < fac) { // forward communication 
                    printf("SEND locallen = %d  rank = %d, expression = %d\n",localLen, rank, (rank*localLen) % (2*fac) );
                    MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX, rank+fac/localLen, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, rank+fac/localLen, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    arrayPlus(vec, recvd, localLen, fac, len);
                } else { // backwards communication
                    printf("Recieve locallen = %d  rank = %d, expression = %d\n",localLen, rank, (rank*localLen) % (2*fac) );
                    MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, rank - fac/localLen, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX, rank- fac/localLen, tag2, MPI_COMM_WORLD);
                    arrayMinus(vec, recvd, localLen, fac, len);
            }
        } else{
            for (int j = 0; j < localLen/fac; j+=2){
                for (int i = 0; i < fac; i++){
                    dummy = vec[j*fac + i] + cexp(-2.0*I*M_PI*i/len) * vec[(j+1)*fac+i];
                    vec[(j+1)*fac+i] = vec[j*fac + i] - cexp(-2.0*I*M_PI*i/len) * vec[(j+1)*fac+i];
                    vec[j*fac + i] = dummy;   
                }
            } 
        }
        if (rank == 0){
            for (int i = 0; i < localLen; i++){
                printf("%f + i%f\n", creal(vec[i]), cimag(vec[i]));
            }
        }
    }
    free(recvd);
}



void pIfft(complex double *vec, int len, int localLen, int rank, int size){
    int tag1 = 1, tag2 = 2;
    printf("got here 3");
    complex double recvd, dummy;
    for (int i = log2(len)-1; i >= 0; i--){
        int fac = pow(2, i);
        for (int j = 0; j < localLen; j++){
            printf("got here 4");
            if (fac > localLen){ // communicate
                if (j % 2*fac < fac) { // forward communication 
                    MPI_Send(&vec[j], 1, MPI_C_DOUBLE_COMPLEX, rank+fac/localLen, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd, 1, MPI_C_DOUBLE_COMPLEX, rank + fac/localLen, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    vec[j] = vec[j] + cexp(2.0*I*M_PI*j/len) * recvd;
                } else { // backwards communication
                    MPI_Recv(&recvd, 1, MPI_C_DOUBLE_COMPLEX, rank - fac/localLen, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Send(&vec[j], 1, MPI_C_DOUBLE_COMPLEX, rank- fac/localLen, tag2, MPI_COMM_WORLD);
                    vec[j] = recvd - cexp(2.0*I*M_PI*(j-fac)/len) * vec[j];
                }
            } else if (j % 2*fac < fac && j + fac < localLen){//no communication
                dummy = vec[j] + cexp(2.0*I*M_PI*(j)/len) * vec[j+ fac]; 
                vec[j+ fac] = vec[j] - cexp(2.0*I*M_PI*(j)/len) * vec[j+ fac];
                vec[j] = dummy;
            }
        }
    } 
    for (int i = 0; i < localLen; i++){
        vec[i] = 1.0/len * vec[i];
    }

}


void writeFile(complex double *vec, int localLen, int rank, int size, char* name){
    int sig = 1;
    int killtag = 100; // 'signal' for process to write it's computed piece and then terminate
    if (rank == 0){
        FILE* file = fopen(name, "w"); 
        for (int i=0;i<localLen;i++){
            // write also first element of file
            fprintf(file, "%f\n",creal(vec[i]));
        }
        fclose(file); 
        MPI_Send(&sig, 1, MPI_INT, 1, killtag, MPI_COMM_WORLD);
    } else{
        MPI_Recv(&sig, 1, MPI_INT, rank-1, killtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        FILE* file = fopen(name, "a");
        for (int i=0;i<localLen;i++){
            fprintf(file, "%f\n",creal(vec[i]));
        }
        fclose(file);
        if (rank<size-1){
            MPI_Send(&sig, 1, MPI_INT, rank+1, killtag, MPI_COMM_WORLD);
        } 
    }
}

void shiftArray(complex double *vec, int localLen, int len, int rank, int size){
    int tag = 1, tag2 = 2;
    complex double dummy;
    if (size == 1){
        return;
    }
    for (int i=0; i<localLen; i++){
        if (rank*localLen < len/2){
            if (i % 2 == 1){
                MPI_Send(&vec[i], 1, MPI_C_DOUBLE_COMPLEX, rank+ len/(2*localLen), tag, MPI_COMM_WORLD);
                MPI_Recv(&vec[i], 1, MPI_C_DOUBLE_COMPLEX, rank + len/(2*localLen), tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } 
        } else {
            if (i%2 ==0 ){
                MPI_Recv(&dummy, 1, MPI_C_DOUBLE_COMPLEX, rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&vec[i], 1, MPI_C_DOUBLE_COMPLEX, rank-1, tag2, MPI_COMM_WORLD);
                vec[i] = dummy;
            }
        }
    }
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
    int N = pow(2, atoi(argv[1]));
    srand(rank+123);
    printf("got here 1 %d\n", atoi(argv[1]));
    int J = (rank < N%size) ? N/size + 1: N/size;
    double complex* randomVec = malloc(J*sizeof(double complex));
    for (int i=0; i < J; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
    }
    printf("got here 2 %d\n", J);
    writeFile(randomVec, J, rank, size, "untouched.txt\0");
    shiftArray(randomVec, J, N, rank, size);
    printf("got here 2 %d\n", J);
    pFft(randomVec, N, J, rank, size);
    writeFile(randomVec, J, rank, size, "after_transform.txt\0");
    // pIfft(randomVec, N, J, rank, size);
    // writeFile(randomVec, J, rank, size, "transed_back.txt\0");
    free(randomVec);
    MPI_Finalize();
    exit(0); 
}