#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


void init_given_2Darray_randomly(int N, float *array){
    int i, j;
    for (i=0; i < N; i++){
        for (j=0; j < N; j++) {
            array[i * N + j] = (rand()/ (float)RAND_MAX) * 2 - 1;
        }
    }
}

void printf_given_1Darray(int N, float *array){
    for (int i = 0; i < N; i++){
        printf("%10f, ", array[i]);
    }
    printf("\n");
}

void printf_given_2Darray(int N, float *array){
    printf("[\n");
    for (int i=0; i < N; i++){
        printf("[");
        for (int j=0; j < N; j++) {
            printf("%6.3f, ", array[i * N + j]);
        }
        printf("],\n");
    }
    printf("]\n");
}

void get_transpose_given2Darray(int N, float *original_array, float *transposed_array){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            transposed_array[j * N + i] = original_array[i * N + j];
        }
    }
}

float get_product_of_row_and_col_vectors(int N, float *row, float *col){
    float total = 0;
    for (int i = 0; i < N; i++){
        total += row[i] * col[i];
    }
    return total;
}

int main(int argc, char **argv) {
    // Variables
    int N = atoi(argv[1]);
    int print_option = atoi(argv[2]);
    int rank, size, i;
    int tag = 0;
    int t;
    double start_time, end_time;
    float *row_of_products;
    float *row_of_a;
    float *col_of_b;
    float *a, *b, *c, *transposed_c, *transposed_b;

    // Get MPI variables
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request request;

    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Assign heap
    row_of_products = (float*)malloc(sizeof(float) * N);
    row_of_a = (float*)malloc(sizeof(float) * N);
    col_of_b = (float*)malloc(sizeof(float) * N);
    a = (float*)malloc(sizeof(float) * N * N);
    b = (float*)malloc(sizeof(float) * N * N);
    c = (float*)malloc(sizeof(float) * N * N);
    transposed_c = (float*)malloc(sizeof(float) * N * N);
    transposed_b = (float*)malloc(sizeof(float) * N * N);

    // Get max number of repeat of mpi.
    // e.g. If the matrix size is N = 8 and np = 2, this will be 8/2 = 4.
    int max_n_blocks = ceil((float)N/(float)size);

    // initialize arrays only when the rank == 0
    if (rank == 0){
        init_given_2Darray_randomly(N, a);
        init_given_2Darray_randomly(N, b);
        get_transpose_given2Darray(N, b, transposed_b);  // get the transpose of b
    }


    // Loop till all matrix elements are calculated.
    for (int block = 0; block < max_n_blocks; block++){
        int block_index = block * size * N; // the index of the block dealt in this repeat

        // Scatter columns from the matrix b. Only needed cols will be copied.
        MPI_Scatter(&transposed_b[block_index], N, MPI_FLOAT, col_of_b, N, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Perform matrix multiplication on each processor
        int row_index = block_index;
        for (int step = 0; step < N; step++) {
            // Get index After rollup. Roll up is performed implicitly because explicit roll up requires more time
            if (step != 0 && rank == 0) {
                row_index = row_index + N;
            }

            // Scatter rows from the matrix a. Only needed rows will be copied.
            MPI_Scatter(&a[row_index], N, MPI_FLOAT, row_of_a, N, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // Perform multiplication on the given row and col
            // Get the col index where the product should be saved
            int col_to_fill = block * size + rank + step;
            if (col_to_fill >= N) {
                col_to_fill = col_to_fill - N;
            }
            row_of_products[col_to_fill] = get_product_of_row_and_col_vectors(N, row_of_a, col_of_b);
        }
        // Gather data from processors to the array on rank=0
        MPI_Gather(row_of_products, N, MPI_FLOAT, &c[block_index], N, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    // Get the final result by taking transpose
    if (rank == 0){
        get_transpose_given2Darray(N, c, transposed_c);
    }

    //You can print if you need
    if (rank == 0 && print_option == 1){
        printf("Matrix a:\n");
        printf_given_2Darray(N, a);
        printf("Matrix b:\n");
        printf_given_2Darray(N, b);
        printf("Product matrix, c:\n");
        printf_given_2Darray(N, transposed_c);
    }

    //free the heap
    free(a);
    free(b);
    free(c);
    free(transposed_b);
    free(transposed_c);


    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    // show the time on the master node
    if (rank == 0) {
        printf("%f\n", end_time - start_time);
    }

    MPI_Finalize();
    return(0);

}

/** for bug check to print row_of_products
if (rank == 1 && repeat == 1) {
    printf("Col to fill: %d\n", col_to_fill);
    printf("Row of a\n");
    printf_given_1Darray(N, row_of_a);
    printf("Col of b\n");
    printf_given_1Darray(N, col_of_b);
    printf("product\n");
    printf("%f\n", get_product_of_row_and_col_vectors(N, row_of_a, col_of_b));
    printf_given_1Darray(N, row_of_products);
}
**/
