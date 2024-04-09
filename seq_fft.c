#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

int reverseInt(int num, int nBits){
    int reversed = 0;
    for (int j = 0; j < nBits; j++){
        if (num & (1 << j)){
            reversed |= 1 << ((nBits -1)-j);
        }
    }

}


int fft(double vec[], int len){
    // shift elements in the list::
    for (int i = 0; i < len/2; i++){
        for 
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
    int N = atoi(argv[1]);
    double randomVec[N];
    for (int i=0; i < N; i++){
        randomVec[i] = (double)rand() / RAND_MAX;
    }
    fft(randomVec);
    
    

    exit(0);
}