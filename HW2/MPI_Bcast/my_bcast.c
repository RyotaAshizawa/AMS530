#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    // Variables
    int rank, size, i;
    double start_time, end_time;


    // Get MPI variables
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;
    MPI_Request request;

    // initialize the array
    int N = pow(2, 14); //datapoints
    float arr[N];
    if (rank==0) {
        for (i = 0; i < N; i++){
            arr[i] = i;
        }
    }

    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();


    //Broad casting and sum
    int source = 0;
    int dest = 1;
    int tag = 0;
    float sum = 0;
    if (rank == 0){
        MPI_Isend(&arr, N, MPI_FLOAT, dest, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        for (i = 0; i < N; i++){
            sum += arr[i];
        }
        printf("my rank %d, sum = %f\n", rank, sum);
    }
    else {
        int prev_rank = rank - 1;
        int next_rank = rank + 1;
        float recv_arr[N];

        MPI_Irecv(&recv_arr, N, MPI_FLOAT, prev_rank, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        for (i = 0; i < N; i++){
            sum += recv_arr[i];
        }
        printf("my rank %d, sum = %f\n", rank, sum);
        if (rank < size - 1) {
            MPI_Isend(&recv_arr, N, MPI_INT, next_rank, tag, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
        }
    }

    //printf("my rank %d, sum = %f\n", rank, sum);

    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    // show the time on the master node
    if (rank==0) {
        printf("%f\n", end_time - start_time);
    }

    MPI_Finalize();
    return(0);
}
