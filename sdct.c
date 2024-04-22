#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>


void fft(complex double  *vec, int len){
    // very basic implementation with several faults, that could damage performance for large inputs..
    if (len <= 1) return;
    complex double* even = malloc(len*sizeof(complex double));
    complex double* odd = malloc(len*sizeof(complex double));
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
    free(even);
    free(odd);
}


void invFft(complex double *vec, int len){
    if (len <= 1) return;
    double complex even[len/2];
    double complex odd[len/2];
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

void invFftHelper(complex double *vec, int len){
    invFft(vec, len);
    for (int i = 0; i < len; i++){
        vec[i] = 1.0/len * vec[i];
    }
}


void dct(complex double *vec, int len){
    /* the discrete cosine transformation, where depending on if vec is in the frequency 
    dimension or in the real dimension we get either frequencies or real numbers respectively*/ 
    complex double* ftilde = malloc(len*sizeof(double complex));
    for (int i = 0; i < len/2; i++){
        ftilde[i] = vec[2*i];
        ftilde[len-i-1] = vec[2*i+1];    
    }
    fft(ftilde, len);
    for (int i=0; i < len;i++){
        vec[i] = 2.0*creal(cexp(-I*M_PI*i/(2.0*len))*ftilde[i]);
    }
    vec[0] /= sqrt(2);
    free(ftilde);
}


void invDct(complex double *vec, int len){
    complex double* ftilde = malloc(len*sizeof(double complex));
    ftilde[0] = (1.0/sqrt(2))  * vec[0];
    ftilde[len/2] = (1.0/sqrt(2)) * vec[len/2];
    for (int i = 1; i < len; i++){
        if (i != len/2){
            ftilde[i] = 0.5*((vec[i]*cos(M_PI*i/(2*len)) + vec[len-i]*sin(M_PI*i/(2*len)))+I*(vec[i]*sin(M_PI*i/(2*len)) - vec[len-i]*cos(M_PI*i/(2*len))));
        }
    }
    invFftHelper(ftilde, len);
    for (int i=0; i < len/2;i++){
        vec[2*i] = ftilde[i];
        vec[2*i+1] = ftilde[len-i-1];
    }
    free(ftilde);
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
    double temp;
    FILE* file1 = fopen("untouched.txt", "r");
    for (int i=0; i < N; i++){
        if (fscanf(file1, "%lf", &temp) != 1){
            printf("failed");
            exit(1);
        }
        randomVec[i] = temp;
    }
    fclose(file1);
    fft(randomVec, N);
    // dct(randomVec, N);
    FILE * file2 = fopen("after_trans.txt", "w");
    for (int i = 0; i < N; i++){
        fprintf(file2, "%f\n", creal(randomVec[i]));
    } 
    fclose(file2);
    // invDct(randomVec,N);
    invFftHelper(randomVec, N);
    FILE* file3 = fopen("transed_back.txt", "w");
    for (int i = 0; i < N; i++){
        fprintf(file3, "%f\n", creal(randomVec[i]));
    } 
    fclose(file3);
    exit(0);
    free(randomVec);
}