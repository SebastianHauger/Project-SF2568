
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
 
 TODO:: 

 3. figure out scheme for pre processing DCT (and implement)
 4. figure out scheme for pre processing inverse DCT (and implement)

 (
    5. start with multidim schemes 
    6. try algorithm on picture files (black and white)
    7. Possibly try scheme on color picture as just three separate schemes
 )
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


void pFft(complex double *vec, int len, int localLen, int rank, int size,int localInd, int rest, int maxLen, complex double *recvd){
    int tag1 = 3, tag2 = 4, sendLen, sendRank;
    complex double dummy;
    // double complex* recvd = malloc(localLen*sizeof(complex double));
    for (int i = 0; i < log2(len); i++){
        int fac = pow(2, i);
        if (fac >= localLen){ // communicate
            if (localInd % (2*fac) < fac) { // forward communication 
                sendLen = (localLen == maxLen) ? ((rank + fac/maxLen < size -rest) ? maxLen: maxLen/2) : maxLen/2;
                if (sendLen == localLen){
                    sendRank = rank + fac/localLen;
                    // printf("fac %d rank %d send rank %d\n\n", fac, rank, sendRank);
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    sendRank = size-rest + (2*(fac+(rank+rest-size)*maxLen))/maxLen;
                    // printf("rank %d send rank %d, \n ", rank, sendRank);
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vec[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                }
                arrayPlusForward(vec, recvd, localLen, fac, len, localInd);
            } else { // backwards communication 
                sendRank = (localLen==maxLen) ? rank - fac/maxLen : ((fac > (rank+ rest-size)*localLen) ? (localInd-fac)/maxLen: rank-fac/localLen); // this is still a bit wrong 
                // printf("rank %d send rank %d, \n ", rank, sendRank);
                MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD);
                arrayMinusForward(vec, recvd, localLen, fac, len, localInd);
            }
        } else{
            for (int j = 0; j < localLen/fac; j+=2){
                for (int k = 0; k < fac; k++){
                    dummy = vec[j*fac + k] + cexp(-2.0*I*M_PI*k/(2*fac)) * vec[(j+1)*fac+k];
                    vec[(j+1)*fac+k] = vec[j*fac + k] - cexp(-2.0*I*M_PI*k/(2*fac)) * vec[(j+1)*fac+k];
                    vec[j*fac + k] = dummy;   
                }
            } 
        }
    }
}


void pIfft(complex double *vec, int len, int localLen, int rank, int size, int localInd, int rest, int maxLen, complex double *recvd){
    int tag1 = 3, tag2 = 4, sendRank, sendLen;
    complex double dummy;
    // printf("Inverse transform \n \n");
    // double complex* recvd = malloc(localLen*sizeof(complex double));
    for (int i = 0; i < log2(len); i++){
        int fac = pow(2, i);
        if (fac >= localLen){ // communicate
            if (localInd % (2*fac) < fac) { // forward communication 
                sendLen = (localLen == maxLen) ? ((rank + fac/maxLen < size -rest) ? maxLen: maxLen/2) : maxLen/2;
                if (sendLen == localLen){
                    sendRank = rank + fac/localLen;
                    // printf("fac %d rank %d send rank %d\n\n", fac, rank, sendRank);
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    sendRank = size-rest + (2*(fac+(rank+rest-size)*maxLen))/maxLen;
                    // printf("rank %d send rank %d, \n ", rank, sendRank);
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vec[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                }
                arrayPlusBackward(vec, recvd, localLen, fac, len, localInd);
            } else { // backwards communication 
                sendRank = (localLen==maxLen) ? rank - fac/maxLen : ((fac > (rank+ rest-size)*localLen) ? (localInd-fac)/maxLen: rank-fac/localLen);
                // printf("rank %d send rank %d \n", rank, sendRank);
                MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX,sendRank, tag2, MPI_COMM_WORLD);
                arrayMinusBackward(vec, recvd, localLen, fac, len, localInd);
            }
        } else{
            for (int j = 0; j < localLen/fac; j+=2){
                for (int k = 0; k < fac; k++){
                    dummy = vec[j*fac + k] + cexp(2.0*I*M_PI*k/(2*fac)) * vec[(j+1)*fac+k];
                    vec[(j+1)*fac+k] = vec[j*fac + k] - cexp(2.0*I*M_PI*k/(2*fac)) * vec[(j+1)*fac+k];
                    vec[j*fac + k] = dummy;   
                }
            } 
        }
    }
    for (int i=0; i < localLen; i++){
        vec[i] *= 1.0/len;
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
                    printf("send 2 from %d to %d with size %d\n", rank, sendRank, i2/2);
                    printf("send 2 from %d to %d with size %d\n", rank, sendRank+1, i2/2);
                    // communicate with several ranks
                    MPI_Send(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    printf("send 1 from %d to %d with size %d\n", rank, sendRank, i2);
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
                printf("Receive from %d to %d with size %d\n", rank, sendRank, i2);
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
                printf("Receive from %d to %d with size %d\n", rank, sendRank, i2);
                MPI_Recv(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                get_backward(vec, copy, localLen, block, chunk, localIndex);    
                MPI_Send(copy, i2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD);
            }
            get_forward(vec, copy, localLen, block, chunk, localIndex, &i1, &i2);
            if (i2 > 0){ // communicate forwards
                sendRank = localIndex + i1 + chunk/2 -block;
                sendRank = (sendRank < ((size-rest)*maxLen)) ? sendRank/maxLen: size-rest + (2*(sendRank-(size-rest)*maxLen))/maxLen; 
                if ((rank < size-rest) && (sendRank >= size-rest) && (localLen-i1)>maxLen/2){
                    printf("send 2 from %d to %d with size %d\n", rank, sendRank, i2/2);
                    printf("send 2 from %d to %d with size %d\n", rank, sendRank+1, i2/2);
                    // communicate with several ranks
                    MPI_Send(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(copy, i2/2, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&copy[i2/2], i2/2, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    printf("send 1 from %d to %d with size %d\n", rank, sendRank, i2);
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
    int logSize = log2(size); 
    int largest = pow(2, log2N-logSize); // largest local size 
    int rest = (size%(N/largest))*2; // number of elements that are not the largest 
    int J = ((size-rank) > rest) ? largest: largest/2;
    int locind = ((size-rank) > rest) ? rank *largest: (size-rest)*largest + (rank+rest-size)*largest/2;
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
        randomVec[i] = rank *largest + i; 
    }
    gettimeofday(&tdr1, NULL);
    writeFile(randomVec, J, rank, size, "untouched.txt\0");
    shiftArray(randomVec, J, N, rank, size, locind, rest, largest, copy);
    // writeFile(randomVec, J, rank, size, "shifted.txt\0");
    if (rank == 0)
        gettimeofday(&tdr2, NULL);
    pFft(randomVec, N, J, rank, size, locind, rest, largest, copy);
    if (rank == 0)
        gettimeofday(&tdr3, NULL);
    writeFile(randomVec, J, rank, size, "after_transform.txt\0");
    shiftArray(randomVec, J, N, rank, size, locind, rest, largest, copy);
    if (rank == 0)
        gettimeofday(&tdr4, NULL);
    pIfft(randomVec, N, J, rank, size, locind, rest, largest, copy);
    if (rank ==0)
        gettimeofday(&tdr5, NULL);
    writeFile(randomVec, J, rank, size, "transed_back.txt\0");
    free(randomVec);
    free(copy);
    MPI_Finalize();
    
    if (rank == 0){
        FILE* file = fopen("timevals.txt", "w");
        timeval_subtract(&restime1, &tdr1, &tdr0);
        timeval_subtract(&restime2, &tdr2, &tdr1);
        timeval_subtract(&restime3, &tdr3, &tdr2);
        timeval_subtract(&restime4, &tdr4, &tdr3);
        timeval_subtract(&restime5, &tdr5, &tdr4);
        fprintf(file, "1: %f \n 2: %f \n 3: %f \n 4: %f \n 5: %f \n", restime1, restime2, restime3, restime4, restime5);
        fclose(file);
    }
    exit(0); 
}