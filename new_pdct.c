
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <complex.h>
#include <mpi.h>


/*
A program to compute the parallel DCT and IDCT using the FFT and IFFT

The program takes (at the moment) and extra agrument (more than standard and MPI arguments)
this argument is the exponent of two which defines the length of the list
 
 */

int timeval_subtract (double *result, struct timeval *x, struct timeval *y){
    struct timeval result0;
    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }
    /* Compute the time remaining to wait.
    tv_usec is certainly positive. */
    result0.tv_sec = x->tv_sec - y->tv_sec;
    result0.tv_usec = x->tv_usec - y->tv_usec;
    *result = ((double)result0.tv_usec)/1e6 + (double)result0.tv_sec;
    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}


void arrayPlusForward(complex double *vec1, complex double *vec2, int localLength, int factor, int length, int lI){
    for (int k = 0; k < localLength; k++){
        vec1[k] = vec1[k] + cexp(-2.0*I*M_PI*((k+lI)%factor)/(2*factor)) * vec2[k];
    }
}


void arrayMinusForward(complex double *vec1, complex double *vec2, int localLength, int factor, int length, int lI){
    for (int k=0; k < localLength; k++){
        vec1[k] = vec2[k] - cexp(-2.0*I*M_PI*((k+lI)%factor)/(2*factor)) * vec1[k];
    }
}


void arrayPlusBackward(complex double *vec1, complex double *vec2, int localLength, int factor, int length, int lI){
    for (int k = 0; k < localLength; k++){
        vec1[k] = vec1[k] + cexp(2.0*I*M_PI*((k+lI)%factor)/(2*factor)) * vec2[k];
    }
}


void arrayMinusBackward(complex double *vec1, complex double *vec2, int localLength, int factor, int length, int lI){
    for (int k=0; k < localLength; k++){
        vec1[k] = vec2[k] - cexp(2.0*I*M_PI*((k+lI)%factor)/(2*factor)) * vec1[k];
    }
}


void pFft(complex double *vec, int len, int localLen, int rank, int size,int localInd, complex double *copy){
    int tag1 = 3, tag2 = 4, sendLen, sendRank, fac, i1, i2, sendInd;
    complex double dummy;
    for (int i = 0; i < log2(len); i++){
        fac = pow(2, i);
        if (rank % 2 == 0){
            // iterate forwards through the array
            i1 = -1, i2 = -1;
            for (int j = 0; j < localLen; j++){
                if (((localInd +j) % (2*fac)) != ((localInd + j) % fac)){
                    //backwards
                    sendInd = localInd + j - fac;
                    if (sendInd < localInd){
                        if (i1 == -1){
                            i1 = j;
                        }
                        if ((sendInd+1) % fac == 0){
                            MPI_Recv(copy, sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            MPI_Send(&vec[i1], sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag1, MPI_COMM_WORLD);
                            arrayMinusForward(&vec[i1], copy, sendLen, fac, len, localInd);
                        }
                    }
                } else {
                    // forwards
                    sendInd = localInd + j + fac;
                    if (sendInd > (localInd + localLen)){
                        if (i1 == -1){
                            i1 == j;
                        }
                        if ((sendInd+1) % fac == 0){
                            sendLen = j-i1;
                            MPI_Send(&vec[i1], sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag1, MPI_COMM_WORLD);
                            MPI_Recv(copy, sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            arrayPlusForward(&vec[i1], copy, sendLen, fac, len, localInd);
                            i1 == -1;
                        }
                    } else {
                        dummy = vec[j] + cexp(-2.0*I*M_PI*(j%fac)/(2*fac)) * vec[sendInd];
                        vec[sendInd] = vec[j] - cexp(-2.0*I*M_PI*(j%fac)/(2*fac)) * vec[sendInd];
                        vec[j] = dummy; 
                    }
                } 
            }
        } else {
            // iterate backwards through the array
            i1 = 0; i2 = 0;
            for (int j = localLen-1; j >= 0; j--){
                if (((localInd + j) % (2*fac)) == ((localInd + j) % fac)){
                    // forwards
                    sendInd = localInd + j + fac;
                    if (sendInd > (localInd + localLen)){
                        if (i1 == -1){
                            i1 == j;
                        }
                        if ((sendInd+1) % fac == 0){
                            sendLen = i1-j;
                            MPI_Send(&vec[j], sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag1, MPI_COMM_WORLD);
                            MPI_Recv(copy, sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            arrayPlusForward(&vec[j], copy, sendLen, fac, len, localInd);
                            i1 == -1;
                        }
                    } else {
                        dummy = vec[j] + cexp(-2.0*I*M_PI*(j%fac)/(2*fac)) * vec[sendInd];
                        vec[sendInd] = vec[j] - cexp(-2.0*I*M_PI*(j%fac)/(2*fac)) * vec[sendInd];
                        vec[j] = dummy; 
                    }
                    
                } else {
                    //backwards
                    sendInd = localInd + j - fac;
                    if (sendInd < localInd){
                        if (i1 == -1){
                            i1 = j;
                        }
                        if ((sendInd+1) % fac == 0){
                            sendLen = i1 -j;
                            MPI_Recv(copy, sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            MPI_Send(&vec[j], sendLen, MPI_C_DOUBLE_COMPLEX, rank+1, tag1, MPI_COMM_WORLD);
                            arrayMinusForward(&vec[j], copy, sendLen, fac, len, localInd);
                        }
                    }
                }
            } 

        }
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
        if (size > 1){
            MPI_Send(&sig, 1, MPI_INT, 1, killtag, MPI_COMM_WORLD);
        }
    } else {
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



void get_forward(complex double *vec, complex double *copy, int lL, int block, int chunk, int lI, int *i1, int *i2){
    // get elements that should be sent forward
    *i2 = 0;
    for (int j = 0; j < lL; j++){
        if (((j+ lI)%chunk < chunk/2) && ((j+lI)%(2*block) >= block) && (j+chunk/2 >= lL)){
            // printf("change %d\n", j+lI);
            if (*i2 == 0)
                *i1 = j;
            copy[*i2] = vec[j];
            (*i2)++;
        }
    }
}

void get_backward(complex double *vec, complex double *copy, int lL, int block, int chunk, int lI){
    // exhange elements of the array in order to free space before communication  
    int i2 = 0;
    complex double dummy; // in order 
    for (int j = 0; j < lL; j++){
        if (((j+ lI)%chunk >= chunk/2) && (j-chunk/2 < 0) && ((j+lI)%(2*block) < block)){
            // printf("here we get the wrong indexes %i \n", j+lI );
            // printf("received %f, sending %f\n", creal(copy[i2]), creal(vec[j]));
            dummy = vec[j];
            vec[j] = copy[i2];
            copy[i2] = dummy;
            i2++;
        }
    }
}


void shiftArray(complex double *vec, int localLen, int len, int rank, int size, int localIndex, int rest, int maxLen, complex double *copy){
    // perform local sorting step:: we perform a smaller sort and then merge..
    int chunk, block, tag1=1, tag2=2, sendRank, i1, i2;
    complex double dummy;
    for (int i = 0; i < (int)log2(len)/2; i++){
        chunk = pow(2, log2(len)-i);
        block = pow(2, i);
        // for (int j=0; j<localLen; j++){
        //     // printf("%f\n", creal(vec[j]));
        // }
        for (int j = 0; j < localLen; j++){
            if (((j+ localIndex)%chunk < chunk/2) && ((j+localIndex)%(2*block) >= block) && (j+chunk/2 < localLen)){
                // printf("change %d with %d\n", j, j + chunk/2-block);
                dummy = vec[j];
                vec[j] = vec[j + chunk/2-block];
                vec[j+ chunk/2-block] = dummy;
            }
        }
        if (rank%2 == 0){
            get_forward(vec, copy, localLen, block, chunk, localIndex, &i1, &i2);
            if (i2 > 0){ // communicate forwards
                sendRank = localIndex + i1 + chunk/2 -block;
                sendRank = (sendRank < ((size-rest)*maxLen)) ? sendRank/maxLen: size-rest + (2*(sendRank-(size-rest)*maxLen))/maxLen; 
                if ((rank < size-rest) && (sendRank >= size-rest) && (localLen - i1)>maxLen/2){
                    // printf("send 2 from %d to %d with size %d\n", rank, sendRank, i2/2);
                    // printf("send 2 from %d to %d with size %d\n", rank, sendRank+1, i2/2);
                    // communicate with several ranks
                    MPI_Send(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    // printf("send 1 from %d to %d with size %d\n", rank, sendRank, i2);
                    MPI_Send(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                i2 = 0;
                for (int j = 0; j < localLen; j++){
                    if (((j+ localIndex)%chunk < chunk/2) && ((j+localIndex)%(2*block) >= block) && (j+chunk/2 >= localLen)){
                        // printf("recvd %f have %f, rank = %d\n", creal(copy[i2]), creal(vec[j]), rank); 
                        // printf("change %d\n", j+localIndex);
                        vec[j] = copy[i2];
                        i2++;
                    }
                }
            }
            i2 = 0;
            for (int j = 0; j<localLen; j++){
                if (((j+ localIndex)%chunk >= chunk/2) && (j-chunk/2 < 0) && ((j+localIndex)%(2*block) < block)){
                    // printf("change %d\n", j+localIndex);
                    if (i2 == 0)
                        i1 = j;
                    i2++;
                }
            } 
            if (i2 > 0){ // communicate backwards
                sendRank = localIndex + i1 - chunk/2 + block;
                sendRank = (sendRank < ((size-rest)*maxLen)) ? sendRank/maxLen: size-rest + (2*(sendRank-(size-rest)*maxLen))/maxLen; 
                // printf("first %d\n", (localIndex - chunk)/(2*maxLen));
                // printf("second %d\n", (localIndex- chunk/2)/maxLen);
                // printf("third %d\n", rank-chunk/(2*localLen));
                // printf("Receive from %d to %d with size %d\n", rank, sendRank, i2);
                MPI_Recv(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                get_backward(vec, copy, localLen, block, chunk, localIndex); 
                MPI_Send(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD);
            }
        } else{
            i2 = 0;
            for (int j = 0; j<localLen; j++){
                if (((j+ localIndex)%chunk >= chunk/2) && (j-chunk/2 < 0) && ((j+localIndex)%(2*block) < block)){
                    // printf("change %d\n", j+localIndex);
                    if (i2 == 0)
                        i1 = j;
                    i2++;
                }
            } 
            if (i2 > 0){ // communicate backwards
                sendRank = localIndex + i1 - chunk/2 + block;
                sendRank = (sendRank < ((size-rest)*maxLen)) ? sendRank/maxLen: size-rest + (2*(sendRank-(size-rest)*maxLen))/maxLen; 
                // printf("first %d\n", (localIndex - chunk)/(2*maxLen));
                // printf("second %d\n", (localIndex- chunk/2)/maxLen);
                // printf("third %d\n", rank-chunk/(2*localLen));
                // printf("Receive from %d to %d with size %d\n", rank, sendRank, i2);
                MPI_Recv(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                get_backward(vec, copy, localLen, block, chunk, localIndex);    
                MPI_Send(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD);
            }
            get_forward(vec, copy, localLen, block, chunk, localIndex, &i1, &i2);
            if (i2 > 0){ // communicate forwards
                sendRank = localIndex + i1 + chunk/2 -block;
                sendRank = (sendRank < ((size-rest)*maxLen)) ? sendRank/maxLen: size-rest + (2*(sendRank-(size-rest)*maxLen))/maxLen; 
                if ((rank < size-rest) && (sendRank >= size-rest) && (localLen-i1)>maxLen/2){
                    // printf("send 2 from %d to %d with size %d\n", rank, sendRank, i2/2);
                    // printf("send 2 from %d to %d with size %d\n", rank, sendRank+1, i2/2);
                    // communicate with several ranks
                    MPI_Send(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    // printf("send 1 from %d to %d with size %d\n", rank, sendRank, i2);
                    MPI_Send(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                i2 = 0;
                for (int j = 0; j < localLen; j++){
                    // printf("recvd %f have %f, rank = %d\n", creal(copy[j]), creal(vec[j]), rank);
                    if (((j+ localIndex)%chunk < chunk/2) && ((j+localIndex)%(2*block) >= block) && (j+chunk/2 >= localLen)){
                        // printf("change %d\n", j+localIndex);
                        vec[j] = copy[i2];
                        i2++;
                    }
                }
            }
        }
    }
}


int main(int argc, char **argv){ 
    // argc is argument count and argv is a list of arguments..
    int rank, size;
    struct timeval tdr0, tdr1, tdr2, tdr3, tdr4, tdr5;
    double restime1, restime2, restime3, restime4, restime5;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if (argc < 2){
        printf("Not enough input arguments");
        exit(0);
    }
    // start timer 
    if (rank == 0)
        gettimeofday(&tdr0, NULL);

    // set up environment:: this might cause errors still although it is quite properly checked 
    int log2N = atoi(argv[1]);
    int N = pow(2, log2N);
    int J = (rank < N%size) ? N/size + 1: N/size; 
    int lI = (rank < N%size) ? J*rank: J*rank + N%size;
    double complex* randomVec = malloc(J*sizeof(double complex));

    if (randomVec == NULL) {
        fprintf(stderr, "Failed to allocate memory for randomVec\n");
        exit(1);
    }
    double complex* copy = malloc(J*sizeof(double complex));
    if (copy == NULL) {
        fprintf(stderr, "Failed to allocate memory for copy\n");
        exit(1);
    }
    for (int i=0; i < J; i++){
        randomVec[i] = rank *  + i; 
    }
    gettimeofday(&tdr1, NULL);
    // writeFile(randomVec, J, rank, size, "untouched.txt\0");
    // shiftArray(randomVec, J, N, rank, size, lI, rest, largest, copy);
    // writeFile(randomVec, J, rank, size, "shifted.txt\0");
    if (rank == 0)
        gettimeofday(&tdr2, NULL);
    pFft(randomVec, N, J, rank, size, lI, copy);
    if (rank == 0)
        gettimeofday(&tdr3, NULL);
    // writeFile(randomVec, J, rank, size, "after_transform.txt\0");
    // shiftArray(randomVec, J, N, rank, size, lI, copy);
    if (rank == 0)
        gettimeofday(&tdr4, NULL);
    // pIfft(randomVec, N, J, rank, size, lI,  copy);
    if (rank ==0)
        gettimeofday(&tdr5, NULL);
    // writeFile(randomVec, J, rank, size, "transed_back.txt\0");
    free(randomVec);
    free(copy);
    MPI_Finalize();
    
    if (rank == 0){
        FILE* file = fopen("timevals.txt", "a");
        timeval_subtract(&restime1, &tdr1, &tdr0);
        timeval_subtract(&restime2, &tdr2, &tdr1);
        timeval_subtract(&restime3, &tdr3, &tdr2);
        timeval_subtract(&restime4, &tdr4, &tdr3);
        timeval_subtract(&restime5, &tdr5, &tdr4);
        fprintf(file, "n_processes: %d \n 1: %f \n 2: %f \n 3: %f \n 4: %f \n 5: %f \n", size, restime1, restime2, restime3, restime4, restime5);
        fclose(file);
    }
    exit(0); 
}