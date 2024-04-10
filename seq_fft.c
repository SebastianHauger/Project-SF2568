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
    complex double odd[len/2];
    for (int i = 0; i < len/2; i++ ){
        even[i] = vec[2*i];
        odd[i] = vec[2*i + 1];
    }
    invFft(even, len/2);
    invFft(odd, len/2);
    for (int i = 0; i < len/2; i++){
        complex double w = cexp(2.0*I*M_PI*i/len) * odd[i];
        vec[i] = even[i] + w; 
        vec[i + len/2] = even[i] - w; 
    }
}

void invFft_helper(complex double *vec, int len){
    invFft(vec, len);
    for (int i = 0; i < len; i++){
        vec[i] = 1.0/len * vec[i];
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
    double complex* randomVec = malloc(N*sizeof(double complex));
    FILE* file1 = fopen("untouched.txt", "w");
    for (int i=0; i < N; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
        fprintf(file1, "%f\n", creal(randomVec[i]));
    }
    fclose(file1);
    fft(randomVec, N);
    FILE * file2 = fopen("after_trans.txt", "w");
    for (int i = 0; i < N; i++){
        fprintf(file2, "%f\n", creal(randomVec[i]));
    } 
    fclose(file2);
    invFft_helper(randomVec, N);
    FILE* file3 = fopen("transed_back.txt", "w");
    for (int i = 0; i < N; i++){
        fprintf(file3, "%f\n", creal(randomVec[i]));
    } 
    fclose(file3);
    exit(0);
    free(randomVec);
}