//
// Created by Ryota Ashizawa on 10/2/22.
//

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    // Variables
    int rank, size, i;
    int tag = 0;
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
    // Initialize the array only for the process whose rank == 0
    if (rank==0) {
        for (i = 0; i < N; i++){
            arr[i] = i;
        }
    }

    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    //Broadcasting by using tree structure
    int max_timestep = ceil(log2(size)) - 1; // # of timesteps required to perform broadcasting. -1 for 0-indexing.

    // Loop for each timestep
    int t;
    for (t = 0; t <= max_timestep; t++) { // t: timestep
        int tot_nodes = pow(2, t); // The total number of nodes that the data should be copied in at this time step
        int rank_start = pow(2, t);  // The smallest process id that should be filled in at this time step

        // send data
        int src_rank;
        for (src_rank = 0; src_rank <= tot_nodes - 1; src_rank++){
            int dst_rank = src_rank + tot_nodes;
            if (rank == src_rank && dst_rank < size){
                //printf("Send data from %d to %d\n", src_rank, dst_rank);
                MPI_Isend(&arr, N, MPI_FLOAT, dst_rank, tag, MPI_COMM_WORLD, &request);
            }
        }

        // receive data
        int dst_rank;
        for (dst_rank = rank_start; dst_rank <= (rank_start + tot_nodes - 1); dst_rank++){
            int src_rank = dst_rank - tot_nodes;
            if (rank == dst_rank && dst_rank < size){
                //printf("Recv from %d to %d\n", src_rank, dst_rank);
                MPI_Irecv(&arr, N, MPI_FLOAT, src_rank, tag, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
        }
    }

    // Check the sum
    float sum = 0;
    for (i = 0; i < N; i++){
        sum += arr[i];
    }
    printf("my rank %d, sum = %f\n", rank, sum);

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

