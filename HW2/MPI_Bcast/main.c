#include <stdio.h>
#include <mpi.h>
#include <string.h>

int main(int argc, char **argv) {
    int rank, size, i;
    int arr[10];
    int name_length = 10;
    char name[10];

    // Get MPI variables
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;
    // Get the host name
    MPI_Get_processor_name(name, &name_length);

    if (rank==0) {
        for (i = 0; i < 10; i++){
            arr[i] = i;
            printf("%d", arr[i]);
        }
    }

    MPI_Finalize();
    return(0);
}
