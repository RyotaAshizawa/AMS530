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

void get_block_from_original_matrix(int first_row, int first_col, int N, int n, float *original_matrix, float *block_matrix){
    // n is the number of elements in the block matrix
    // N is the number of elements in the original matrix
    for (int row = first_row; row < first_row + n; row++){
        for (int col = first_col; col < first_col + n; col++){
            block_matrix[n * (row - first_row) + (col - first_col)] = original_matrix[N * row + col];
        }
    }
}

void send_a_to_each_proc(int blocks_per_dir, int step, int N, int n, float *a, float *block_a_to_send, int tag, MPI_Comm comm, MPI_Request *request) {
    int dst_rank;
    for (int row_i = 0; row_i < blocks_per_dir; row_i++) {
        // get target col to broad cast
        int target_col = step + row_i;
        if (target_col >= blocks_per_dir) {
            target_col = target_col - blocks_per_dir;
        }
        get_block_from_original_matrix(row_i * n, target_col * n, N, n, a, block_a_to_send);
        //printf("Step:%d, row_i:%d, target_col:%d\n", step, row_i, target_col);
        //printf_given_2Darray(n, block_a_to_send);

        for (int col_i = 0; col_i < blocks_per_dir; col_i++) {
            dst_rank = row_i * blocks_per_dir + col_i;
            MPI_Isend(block_a_to_send, n * n, MPI_FLOAT, dst_rank, tag, comm, request);
        }
    }
}

void send_rolledup_b_to_each_proc(int n_of_rollup, int blocks_per_dir, int N, int n, float *b, float *block_b_to_send, int tag, MPI_Comm comm, MPI_Request *request){
    int dst_rank;
    for (int block_row = 0; block_row < blocks_per_dir; block_row++) {
        for (int block_col = 0; block_col < blocks_per_dir; block_col++) {
            // get dst_rank
            dst_rank = blocks_per_dir * block_row + block_col;
            // get [Rolled up block matrix data] for sending.
            int rolledup_block_row = block_row + n_of_rollup;
            if (rolledup_block_row >= blocks_per_dir){
                rolledup_block_row = rolledup_block_row - blocks_per_dir;
            }
            get_block_from_original_matrix(rolledup_block_row * n, block_col * n, N, n, b, block_b_to_send);
            // send data
            MPI_Isend(block_b_to_send, n * n, MPI_FLOAT, dst_rank, tag, comm, request);
        }
    }
}

void multiply_matrices(int N, float *a, float *b, float *c){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++) {
                c[i * N + j] += a[i * N + k] * b[k * N + j];
            }
        }
    }
}

void sum_matrices(int N, float *a, float *b){
     // Add a to b
    for (int row = 0; row < N; row++){
        for (int col = 0; col < N; col++){
            b[row * N + col] = a[row * N + col] + b[row * N + col];
        }
    }
}

void subtract_matrices(int N, float *a, float *b){
    //subtract a from b
    for (int row = 0; row < N; row++){
        for (int col = 0; col < N; col++){
            b[row * N + col] = -a[row * N + col] + b[row * N + col];
        }
    }
}

void init_matrix(int N, float *a){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            a[i * N + j] = 0.0;
        }
    }
}

void rearrange_the_gathered_data(int blocks_per_dir, int n, int N, float *gathered_matrix, float *rearranged_matrix) {
    for (int gathered_block = 0; gathered_block < blocks_per_dir * blocks_per_dir; gathered_block++) { // this is linear gathered_block
        int row_of_block = floor((float) gathered_block / (float) blocks_per_dir);
        int col_of_block = gathered_block - row_of_block * blocks_per_dir;
        for (int ele_in_block = 0; ele_in_block < n * n; ele_in_block++) {
            int row_in_block = floor((float) ele_in_block / (float) n);
            int col_in_block = ele_in_block - row_in_block * n;
            int from_this = gathered_block * n * n + ele_in_block;
            int to_this = (row_of_block * N * n + col_of_block * n) + (row_in_block * N + col_in_block);
            //if (from_this == 201) {
            //printf("linear gathered_block:%d, row_of_block:%d, col_of_block:%d\n", gathered_block, row_of_block, col_of_block);
            //printf("ele_in_block:%d, row_in_block:%d, col_in_block:%d\n", ele_in_block, row_in_block, col_in_block);
            //printf("To this:%d\n", to_this);
            //printf("From this:%d\n", from_this);
            //}
            rearranged_matrix[to_this] = gathered_matrix[from_this];
        }
    }
}


int main(int argc, char **argv) {
    // Variables
    int N = atoi(argv[1]);
    int print_option = atoi(argv[2]);
    int rank, size, i;
    int tag = 0;
    int t;
    double start_time, end_time;
    float *a, *b, *c, *c_temp;
    float *block_a_to_send, *block_a_to_recv; //block matrix
    float *block_b_to_send, *block_b_to_recv; //block matrix
    float *block_c, *block_c_store_partial;

    // Get MPI variables
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request request;

    // Time detection
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Get information about block partitioning
    int blocks_per_dir = sqrt(size); //number of blocks in each direction
    int n = N / blocks_per_dir; // number of N in each block
    int max_size = size; // the number of cpu
    // if the divided block is smaller than 2 x 2, use the 2 x 2 block
    if (n < 2){
        blocks_per_dir= N / 2;
        n = N / blocks_per_dir;
        max_size = blocks_per_dir * blocks_per_dir;
    }

    // Assign heap
    a = (float*)malloc(sizeof(float) * N * N);
    b = (float*)malloc(sizeof(float) * N * N);
    c = (float*)malloc(sizeof(float) * N * N);
    c_temp = (float*)malloc(sizeof(float) * N * N);
    block_a_to_send = (float*)malloc(sizeof(float) * n * n);
    block_a_to_recv = (float*)malloc(sizeof(float) * n * n);
    block_b_to_send = (float*)malloc(sizeof(float) * n * n);
    block_b_to_recv = (float*)malloc(sizeof(float) * n * n);
    block_c = (float*)malloc(sizeof(float) * n * n);
    block_c_store_partial = (float*)malloc(sizeof(float) * n * n);

    // Init matrix if initialization is needed
    init_matrix(n, block_c);

    // initialize arrays only when the rank == 0
    if (rank == 0){
        // init
        init_given_2Darray_randomly(N, a);
        init_given_2Darray_randomly(N, b);
    }

    // do matrix multiplication
    // loop for step
    for (int step = 0; step < blocks_per_dir; step++) {
        /**
         *  About matrix A
         */
        // Copy the diagonal block of A to each processor from rank 0
        if (rank == 0) {
            send_a_to_each_proc(blocks_per_dir, step, N, n, a, block_a_to_send, tag, MPI_COMM_WORLD, &request);
        }
        // recv the block matrix of a by each processor
        MPI_Irecv(block_a_to_recv, n * n, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);

        /**
         *  About matrix B
         */
        // Send rolled up matrix B to all processors from rank 0.
        // [!] Roll up is implicitly calculated in the function to save time.
        if (rank == 0) {
            send_rolledup_b_to_each_proc(step, blocks_per_dir, N, n, b, block_b_to_send, tag, MPI_COMM_WORLD, &request);
        }
        // recv the block matrix of b by each processor
        MPI_Irecv(block_b_to_recv, n * n, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);

        /**
         *  Matrix multiplication on block A and block B
         */
        init_matrix(n, block_c_store_partial);
        multiply_matrices(n, block_a_to_recv, block_b_to_recv, block_c_store_partial);
        sum_matrices(n, block_c_store_partial, block_c);

    }
    // Gather summation of each processor to rank 0, and rearrange it on rank 0
    MPI_Gather(block_c, n * n, MPI_FLOAT, &c_temp[rank * n * n], n * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        rearrange_the_gathered_data(blocks_per_dir, n, N, c_temp, c);
    }

    //You can print if you need
    if (rank == 0 && print_option == 1){
        printf("Matrix a:\n");
        printf_given_2Darray(N, a);
        printf("Matrix b:\n");
        printf_given_2Darray(N, b);
        printf("Product c:\n");
        printf_given_2Darray(N, c);
    }

    //free the heap
    free(a);
    free(b);
    free(c);
    free(c_temp);
    free(block_a_to_send);
    free(block_a_to_recv);
    free(block_b_to_send);
    free(block_b_to_recv);
    free(block_c);
    free(block_c_store_partial);

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


/** funtion storage
Print A
if (step == 1) {
//printf("Sent A, rank:%d, step:%d", rank, step);
//printf_given_2Darray(n, block_a_to_recv);
}

Print B
if (step == 1) {
    //printf("Sent B, rank:%d, step:%d", rank, step);
    //printf_given_2Darray(n, block_b_to_recv);
}
**/
