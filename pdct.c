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
    complex double* even = malloc(len*sizeof(complex double));
    complex double* odd = malloc(len*sizeof(complex double));


    for (int i = log2(len)-1; i >= 0; i--){
        int fac = pow(2, i);
        if (fac > localLen){ // communicate
            if (rank % fac < fac) { // forward communication 

            } else { // backwards communication 

            }
        }


    }

}


void fft(complex double  *vec, int len, int localLen){
    // very basic implementation with several faults, that could damage performance for large inputs..
    if (len <= 1) return;
    complex double* even = malloc(len*sizeof(complex double));
    complex double* odd = malloc(len*sizeof(complex double));
    for (int i=0; i<len/2;i++){
        even[i] = vec[2*i];
        odd[i] = vec[i*2+1]; // ok because a factor of two:))
    }
    fft(even, len/2, localLen);
    fft(odd, len/2, localLen);
    for (int i = 0; i < len/2; i++) {
        complex double w = cexp(-2.0*I*M_PI*i/len) * odd[i];
        vec[i] = even[i] + w;
        vec[i + len/2] = even[i] - w;
    }
    free(even);
    free(odd);
}


void invFft(complex double *vec, int len, int localLen){
    if (len <= 1) return;
    double complex even[len/2];
    double complex odd[len/2];
    for (int i = 0; i < len/2; i++ ){
        even[i] = vec[2*i];
        odd[i] = vec[2*i + 1];
    }
    invFft(even, len/2, localLen);
    invFft(odd, len/2, localLen);
    for (int i = 0; i < len/2; i++){
        complex double w = cexp(2.0*I*M_PI*i/len) * odd[i];
        vec[i] = even[i] + w; 
        vec[i + len/2] = even[i] - w; 
    }
}

void invFftHelper(complex double *vec, int len, int localLen){
    invFft(vec, len, localLen);
    for (int i = 0; i < localLen; i++){
        vec[i] = 1.0/len * vec[i];
    }
}


void dct(complex double *vec, int len, int localLen){
    /* the discrete cosine transformation, where depending on if vec is in the frequency 
    dimension or in the real dimension we get either frequencies or real numbers respectively*/ 
    complex double* ftilde = malloc(localLen*sizeof(double complex));
    for (int i = 0; i < len/2; i++){
        ftilde[i] = vec[2*i];
        ftilde[len-i-1] = vec[2*i+1];    
    }
    fft(ftilde, len, localLen);
    for (int i=0; i < len;i++){
        vec[i] = 2.0*creal(cexp(-I*M_PI*i/(2.0*len))*ftilde[i]);
    }
    vec[0] /= sqrt(2);
    free(ftilde);
}


void invDct(complex double *vec, int len, int localLen){
    complex double* ftilde = malloc(len*sizeof(double complex));
    ftilde[0] = (1.0/sqrt(2))  * vec[0];
    ftilde[len/2] = (1.0/sqrt(2)) * vec[len/2];
    for (int i = 1; i < len; i++){
        if (i != len/2){
            ftilde[i] = 0.5*((vec[i]*cos(M_PI*i/(2*len)) + vec[len-i]*sin(M_PI*i/(2*len)))+I*(vec[i]*sin(M_PI*i/(2*len)) - vec[len-i]*cos(M_PI*i/(2*len))));
        }
    }
    invFftHelper(ftilde, len, localLen);
    for (int i=0; i < len/2;i++){
        vec[2*i] = ftilde[i];
        vec[2*i+1] = ftilde[len-i-1];
    }
    free(ftilde);
}


int main(int argc, char **argv){ 
    // argc is argument count and argv is a list of arguments..
    if (argc < 3){
        printf("Not enough input arguments");
        exit(0);
    }
    int rank, size;
    int N = pow(2, atoi(argv[2]));
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL));
    int J = (rank < N%size) ? N/size + 1: N/size;
    double complex* randomVec = malloc(J*sizeof(double complex));
    for (int i=0; i < N; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
    }
    fft(randomVec, N, J);
    MPI_Finalize();
    exit(0);
    free(randomVec);
}