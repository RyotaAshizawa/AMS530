#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    int rank, size, i;
    double start_time, end_time;

    // Get MPI variables
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;


    // initialize the array
    int N = pow(2, 10); //datapoints
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
    MPI_Bcast(arr, N, MPI_FLOAT, 0, MPI_COMM_WORLD);
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
        printf("original broad casting: %f\n", start_time - end_time);
    }
    MPI_Finalize();
    return(0);
}
