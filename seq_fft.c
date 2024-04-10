#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <complex.h>

int reverseInt(int num, int nBits){
    int reversed = 0;
    for (int j = 0; j < nBits; j++){
        if (num & (1 << j)){
            reversed |= 1 << ((nBits -1)-j); 
        }
    }
    return reversed;
}


double fft_from_book(double  vec[], int max, int min, int k, double complex w){
    double even, odd;
    if (max-min>1){
        odd = fft_from_book(vec, (max-min/2), max, k, w);
        even = fft_from_book(vec, min, (max-min)/2, k, w);
    } else {
        return vec[min];
    }
    return even + pow(w, k) * odd;
}


int fft(double vec[], int len){
    // shift elements in the list::
    int k;
    int nBits = log2(len); // there are eight bits in one byte...
    double dummy;
    for (int i = 0; i < len/2; i++){
        k = reverseInt(i, nBits);
        if (k!=i){
            dummy = vec[i];
            vec[i] = vec[k];
            vec[k] = dummy;
        }
    } 
    double new_vec[len];
    const double complex W = cos(2*M_PI/len)-I*sin(2*M_PI/len); // using eulers identity 
    for (int i = 0; i < nBits; i++){
        for (int j = 0; j < len; j++){
            new_vec[i] = vec[i] + pow(W, len); // TODO:: implement solution...
        }
    }
    return 0;
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