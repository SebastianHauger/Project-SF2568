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

 3. figure out scheme for pre processing DCT (and implement)
 4. figure out scheme for pre processing inverse DCT (and implement)

 (
    5. start with multidim schemes 
    6. try algorithm on picture files (black and white)
    7. Possibly try scheme on color picture as just three separate schemes
 )
 */

void arrayPlusForward(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int k = 0; k < localLength; k++){
        vec1[k] = vec1[k] + cexp(-2.0*I*M_PI*(k%factor)/(2*factor)) * vec2[k];
    }
}


void arrayMinusForward(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int k=0; k < localLength; k++){
        vec1[k] = vec2[k] - cexp(-2.0*I*M_PI*(k%factor)/(2*factor)) * vec1[k];
    }
}


void arrayPlusBackward(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int k = 0; k < localLength; k++){
        vec1[k] = vec1[k] + cexp(2.0*I*M_PI*(k%factor)/(2*factor)) * vec2[k];
    }
}


void arrayMinusBackward(complex double *vec1, complex double *vec2, int localLength, int factor, int length){
    for (int k=0; k < localLength; k++){
        vec1[k] = vec2[k] - cexp(2.0*I*M_PI*(k%factor)/(2*factor)) * vec1[k];
    }
}


void pFft(complex double *vec, int len, int localLen, int rank, int size,int localInd, int rest, int maxLen){
    int tag1 = 1, tag2 = 2, sendLen, sendRank;
    complex double dummy;
    double complex* recvd = malloc(localLen*sizeof(complex double));
    for (int i = 0; i < log2(len); i++){
        int fac = pow(2, i);
        if (fac >= localLen){ // communicate
            if (localInd % (2*fac) < fac) { // forward communication 
                sendLen = (localLen == maxLen) ? ((rank + fac/maxLen < size -rest) ? maxLen: maxLen/2) : maxLen/2;
                if (sendLen == localLen){
                    sendRank = rank + fac/localLen;
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    sendRank = size-rest + 2*(fac+(rank+rest-size)*maxLen)/maxLen;
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vec[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                }
                sendRank = (localLen == maxLen) ? ((rank + fac/maxLen < size-rest) ? rank + fac/maxLen: size-rest + 2*(fac - (rank+rest-size)*maxLen)/maxLen): rank + fac/localLen;
                arrayPlusForward(vec, recvd, localLen, fac, len);
            } else { // backwards communication
                printf("Recieve locallen = %d  rank = %d, expression = %d\n",localLen, rank, (rank*localLen) % (2*fac) );
                MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, rank - fac/localLen, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX, rank- fac/localLen, tag2, MPI_COMM_WORLD);
                arrayMinusForward(vec, recvd, localLen, fac, len);
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
    free(recvd);
}


void pIfft(complex double *vec, int len, int localLen, int rank, int size, int localInd, int rest, int maxLen){
    int tag1 = 1, tag2 = 2, sendRank, sendLen;
    complex double dummy;
    double complex* recvd = malloc(localLen*sizeof(complex double));
    for (int i = 0; i < log2(len); i++){
        int fac = pow(2, i);
        if (fac >= localLen){ // communicate
            if (localInd % (2*fac) < fac) { // forward communication 
                sendLen = (localLen == maxLen) ? ((rank + fac/maxLen < size -rest) ? maxLen: maxLen/2) : maxLen/2;
                if (sendLen == localLen){
                    sendRank = rank + fac/localLen;
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    sendRank = size-rest + 2*(fac+(rank+rest-size)*maxLen)/maxLen;
                    MPI_Send(vec, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag1, MPI_COMM_WORLD);
                    MPI_Recv(recvd, sendLen, MPI_C_DOUBLE_COMPLEX, sendRank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vec[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag1, MPI_COMM_WORLD);
                    MPI_Recv(&recvd[localLen/2], sendLen, MPI_C_DOUBLE_COMPLEX, sendRank+1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                }
                sendRank = (localLen == maxLen) ? ((rank + fac/maxLen < size-rest) ? rank + fac/maxLen: size-rest + 2*(fac - (rank+rest-size)*maxLen)/maxLen): rank + fac/localLen;
                arrayPlusForward(vec, recvd, localLen, fac, len);
            } else { // backwards communication
                printf("Recieve locallen = %d  rank = %d, expression = %d\n",localLen, rank, (rank*localLen) % (2*fac) );
                MPI_Recv(recvd, localLen, MPI_C_DOUBLE_COMPLEX, rank - fac/localLen, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Send(vec, localLen, MPI_C_DOUBLE_COMPLEX, rank- fac/localLen, tag2, MPI_COMM_WORLD);
                arrayMinusForward(vec, recvd, localLen, fac, len);
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
    free(recvd);
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


unsigned int reverseBits(unsigned int num, int len){
    unsigned int new=0;
    for(int i =0; i<len; i++){
        if ((num & 1<<i)){
            new |= 1<<((len-1)-i);
        } 
    }
    return new;
}


void shiftArray(complex double *vec, int localLen, int len, int rank, int size, int localIndex, int rest, int maxLen){
    MPI_Status status;
    int sender_rank;
    int tag = 1, tag2 = 2;
    complex double dummy;
    unsigned int index;
    int process;
    for (int i=0; i<localLen; i++){
        index = reverseBits(i+localIndex, log2(len));
        if (index >= localIndex && index < localIndex + localLen && i < index%localLen){
            // no communication. just swap local indexes.. i<index ensures no swap back
            dummy = vec[i];
            vec[i] = vec[index%localLen];
            vec[index%localLen] = dummy;
        }
        else if (localIndex + localLen <= index){
            // controling rank, from here one can send and receive specific elements
            int corrRank = (index < (size-rest)*maxLen) ? index/maxLen: size-rest + 2*(index-(size-rest)*maxLen)/maxLen; // might be a issue 
            MPI_Send(&index, 1, MPI_INT, corrRank, tag, MPI_COMM_WORLD);
            MPI_Send(&vec[i], 1, MPI_C_DOUBLE_COMPLEX, corrRank, tag2, MPI_COMM_WORLD);
            MPI_Recv(&vec[i], 1, MPI_C_DOUBLE_COMPLEX, corrRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (index < localIndex) {
            MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            sender_rank = status.MPI_SOURCE;
            index = index % localLen;
            MPI_Recv(&dummy, 1, MPI_C_DOUBLE_COMPLEX, sender_rank, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&vec[index], 1, MPI_C_DOUBLE_COMPLEX, sender_rank, tag, MPI_COMM_WORLD);
            vec[index] = dummy;
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
    // set up environment:: this might cause errors still although it is quite properly checked 
    int log2N = atoi(argv[1]);
    int N = pow(2, log2N);
    printf("got here 1 %d\n", atoi(argv[1]));
    int logSize = log2(size); 
    int largest = pow(2, log2N-logSize); // largest local size 
    int rest = (size%(N/largest))*2; // number of elements that are not the largest 
    int J = ((size-rank) > rest) ? largest: largest/2;
    int locind = ((size-rank) > rest) ? rank *largest: (size-rest)*largest + (rank+rest-size)*largest/2;
    double complex* randomVec = malloc(J*sizeof(double complex));
    
    for (int i=0; i < J; i++){
        randomVec[i] = i; 
    }
    writeFile(randomVec, J, rank, size, "untouched.txt\0");
    printf("got here 2");
    shiftArray(randomVec, J, N, rank, size, locind, rest, largest);
    pFft(randomVec, N, J, rank, size, locind, rest, largest);
    writeFile(randomVec, J, rank, size, "after_transform.txt\0");
    shiftArray(randomVec, J, N, rank, size, locind, rest, largest);
    pIfft(randomVec, N, J, rank, size, locind, rest, largest);
    writeFile(randomVec, J, rank, size, "transed_back.txt\0");
    free(randomVec);
    MPI_Finalize();
    exit(0); 
}