#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>



void fft(complex double  *vec, int len){
    // very basic implementation with several faults, that could damage performance for large inputs..
    if (len <= 1) return;
    complex double even[len/2]; 
    complex double odd[len/2];
    for (int i=0; i<len/2;i++){
        even[i] = vec[2*i];
        odd[i] = vec[i*2+1]; // ok because a factor of two:))
    }
    fft(even, len/2);
    fft(odd, len/2);
    for (int i = 0; i < len/2; i++) {
        complex double w = cexp(-2.0*I*M_PI*i/len) * odd[i];
        vec[i] = even[i] + w;
        vec[i + len/2] = even[i] - w;
    }
}

void invFft(complex double *vec, int len){
    if (len <= 1) return;
    complex double even[len/2];
    complex double even[len/2];
    for (int i = 0; i < len/2; i++ )

}

void invFft_helper(complex double *vec, int len){
    invFft(vec, len);
    for (int i = 0; i < len; i++){
        vec[i] = 1/(2*M_PI) * vec[i];
    }
}


int main(int argc, char **argv){ 
    // argc is argument count and argv is a list of arguments..
    if (argc < 2){
        printf("No input arguments");
        exit(0);
    }
    srand(time(NULL));
    int N = pow(2, atoi(argv[1]));
    double randomVec[N];
    for (int i=0; i < N; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
        // printf("original %f\n", randomVec[i]);
    }
    fft(randomVec, N);
    for (int i = 0; i < N; i++){
        // printf("new %f\n", randomVec[i]);
    } 
    exit(0);
}