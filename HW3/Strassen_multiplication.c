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

void copy_matrix(int N, float *original, float *target){
    // copy matrix
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) {
            target[i * N + j] = original[i * N + j];
        }
    }
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
    /**
     *  Add a to b
     */
    for (int row = 0; row < N; row++){
        for (int col = 0; col < N; col++){
            b[row * N + col] = a[row * N + col] + b[row * N + col];
        }
    }
}

void subtract_matrices(int N, float *a, float *b){
    /**
     *  subtract a from b
     */
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


void rearrange_block_matrix(float *a, float *b, float *c, float *d, float *rearranged_matrix, int n, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i < n && j < n){
                rearranged_matrix[i * N + j] = a[i * n + j];
            }
            if (i < n && j >= n){
                //printf("To %d, from %d\n", i * N + j, i * n + (j - n));
                rearranged_matrix[i * N + j] = b[i * n + (j - n)];
            }
            if (i >= n && j < n){
                rearranged_matrix[i * N + j] = c[(i - n) * n + j];
            }
            if (i >= n && j >= n){
                rearranged_matrix[i * N + j] = d[(i - n) * n + (j - n)];
            }
        }
    }
}

void recursive_naive(float *a, float *b, float *c, int n){
    // if n == 2, return the matrix multiplication by using naive method
    if (n == 2){
        multiply_matrices(n, a, b, c);
        return;
    }
    // half the n
    int half_n = n/2;
    // Get block matrix
    float *a11;
    float *a12;
    float *a21;
    float *a22;
    float *b11;
    float *b12;
    float *b21;
    float *b22;
    float *c111;
    float *c121;
    float *c211;
    float *c221;
    float *c112;
    float *c122;
    float *c212;
    float *c222;

    a11 = (float*)malloc(sizeof(float) * half_n * half_n);
    a12 = (float*)malloc(sizeof(float) * half_n * half_n);
    a21 = (float*)malloc(sizeof(float) * half_n * half_n);
    a22 = (float*)malloc(sizeof(float) * half_n * half_n);
    b11 = (float*)malloc(sizeof(float) * half_n * half_n);
    b12 = (float*)malloc(sizeof(float) * half_n * half_n);
    b21 = (float*)malloc(sizeof(float) * half_n * half_n);
    b22 = (float*)malloc(sizeof(float) * half_n * half_n);
    c111 = (float*)malloc(sizeof(float) * half_n * half_n);
    c121 = (float*)malloc(sizeof(float) * half_n * half_n);
    c211 = (float*)malloc(sizeof(float) * half_n * half_n);
    c221 = (float*)malloc(sizeof(float) * half_n * half_n);
    c112 = (float*)malloc(sizeof(float) * half_n * half_n);
    c122 = (float*)malloc(sizeof(float) * half_n * half_n);
    c212 = (float*)malloc(sizeof(float) * half_n * half_n);
    c222 = (float*)malloc(sizeof(float) * half_n * half_n);

    get_block_from_original_matrix(0, 0, n, half_n, a, a11);
    get_block_from_original_matrix(0, half_n, n, half_n, a, a12);
    get_block_from_original_matrix(half_n, 0, n, half_n, a, a21);
    get_block_from_original_matrix(half_n, half_n, n, half_n, a, a22);
    get_block_from_original_matrix(0, 0, n, half_n, b, b11);
    get_block_from_original_matrix(0, half_n, n, half_n, b, b12);
    get_block_from_original_matrix(half_n, 0, n, half_n, b, b21);
    get_block_from_original_matrix(half_n, half_n, n, half_n, b, b22);


    // recursive
    recursive_naive(a11, b11, c111, half_n);
    recursive_naive(a12, b21, c112, half_n);
    recursive_naive(a11, b12, c121, half_n);
    recursive_naive(a12, b22, c122, half_n);
    recursive_naive(a21, b11, c211, half_n);
    recursive_naive(a22, b21, c212, half_n);
    recursive_naive(a21, b12, c221, half_n);
    recursive_naive(a22, b22, c222, half_n);
    sum_matrices(half_n, c112, c111);
    sum_matrices(half_n, c122, c121);
    sum_matrices(half_n, c212, c211);
    sum_matrices(half_n, c222, c221);
    //jprintf_given_2Darray(half_n, c111);
    //printf_given_2Darray(half_n, c121);
    //printf_given_2Darray(half_n, c211);
    //printf_given_2Darray(half_n, c221);

    // reorder the matrix
    rearrange_block_matrix(c111, c121, c211, c221, c, half_n, n);

    // free the matrix
    free(a11);
    free(a12);
    free(a21);
    free(a22);
    free(b11);
    free(b12);
    free(b21);
    free(b22);
}

void do_strassen(float *a, float *b, float *c, int n){
    // if n == 2, return the matrix multiplication by using naive method
    if (n == 2){
        multiply_matrices(n, a, b, c);
        return;
    }
    // half the n
    int half_n = n/2;
    // Get block matrix
    float *a11;
    float *a12;
    float *a21;
    float *a22;
    float *b11;
    float *b12;
    float *b21;
    float *b22;
    float *c11;
    float *c12;
    float *c21;
    float *c22;
    float *m1;
    float *m2;
    float *m3;
    float *m4;
    float *m5;
    float *m6;
    float *m7;
    float *A11_A22;
    float *B11_B22;
    float *A21_A22;
    float *B12_B22;
    float *B21_B11;
    float *A11_A12;
    float *A21_A11;
    float *B11_B12;
    float *A12_A22;
    float *B21_B22;

    a11 = (float*)malloc(sizeof(float) * half_n * half_n);
    a12 = (float*)malloc(sizeof(float) * half_n * half_n);
    a21 = (float*)malloc(sizeof(float) * half_n * half_n);
    a22 = (float*)malloc(sizeof(float) * half_n * half_n);
    b11 = (float*)malloc(sizeof(float) * half_n * half_n);
    b12 = (float*)malloc(sizeof(float) * half_n * half_n);
    b21 = (float*)malloc(sizeof(float) * half_n * half_n);
    b22 = (float*)malloc(sizeof(float) * half_n * half_n);
    c11 = (float*)malloc(sizeof(float) * half_n * half_n);
    c12 = (float*)malloc(sizeof(float) * half_n * half_n);
    c21 = (float*)malloc(sizeof(float) * half_n * half_n);
    c22 = (float*)malloc(sizeof(float) * half_n * half_n);
    m1  = (float*)malloc(sizeof(float) * half_n * half_n);
    m2  = (float*)malloc(sizeof(float) * half_n * half_n);
    m3  = (float*)malloc(sizeof(float) * half_n * half_n);
    m4  = (float*)malloc(sizeof(float) * half_n * half_n);
    m5  = (float*)malloc(sizeof(float) * half_n * half_n);
    m6  = (float*)malloc(sizeof(float) * half_n * half_n);
    m7  = (float*)malloc(sizeof(float) * half_n * half_n);
    A11_A22 = (float*)malloc(sizeof(float) * half_n * half_n);
    B11_B22 = (float*)malloc(sizeof(float) * half_n * half_n);
    A21_A22 = (float*)malloc(sizeof(float) * half_n * half_n);
    B12_B22 = (float*)malloc(sizeof(float) * half_n * half_n);
    B21_B11 = (float*)malloc(sizeof(float) * half_n * half_n);
    A11_A12 = (float*)malloc(sizeof(float) * half_n * half_n);
    A21_A11 = (float*)malloc(sizeof(float) * half_n * half_n);
    B11_B12 = (float*)malloc(sizeof(float) * half_n * half_n);
    A12_A22 = (float*)malloc(sizeof(float) * half_n * half_n);
    B21_B22 = (float*)malloc(sizeof(float) * half_n * half_n);

    get_block_from_original_matrix(0, 0, n, half_n, a, a11);
    get_block_from_original_matrix(0, half_n, n, half_n, a, a12);
    get_block_from_original_matrix(half_n, 0, n, half_n, a, a21);
    get_block_from_original_matrix(half_n, half_n, n, half_n, a, a22);
    get_block_from_original_matrix(0, 0, n, half_n, b, b11);
    get_block_from_original_matrix(0, half_n, n, half_n, b, b12);
    get_block_from_original_matrix(half_n, 0, n, half_n, b, b21);
    get_block_from_original_matrix(half_n, half_n, n, half_n, b, b22);

    // Get sums
    sum_matrices(half_n, a11, A11_A22);
    sum_matrices(half_n, a22, A11_A22);
    sum_matrices(half_n, b11, B11_B22);
    sum_matrices(half_n, b22, B11_B22);
    sum_matrices(half_n, a21, A21_A22);
    sum_matrices(half_n, a22, A21_A22);
    sum_matrices(half_n, b12, B12_B22);
    subtract_matrices(half_n, b22, B12_B22);
    sum_matrices(half_n, b21, B21_B11);
    subtract_matrices(half_n, b11, B21_B11);
    sum_matrices(half_n, a11, A11_A12);
    sum_matrices(half_n, a12, A11_A12);
    sum_matrices(half_n, a21, A21_A11);
    subtract_matrices(half_n, a11, A21_A11);
    sum_matrices(half_n, b11, B11_B12);
    sum_matrices(half_n, b12, B11_B12);
    sum_matrices(half_n, a12, A12_A22);
    subtract_matrices(half_n, a22, A12_A22);
    sum_matrices(half_n, b21, B21_B22);
    sum_matrices(half_n, b22, B21_B22);


    // recursive
    do_strassen(A11_A22, B11_B22, m1, half_n);
    do_strassen(A21_A22, b11, m2, half_n);
    do_strassen(a11, B12_B22, m3, half_n);
    do_strassen(a22, B21_B11, m4, half_n);
    do_strassen(A11_A12, b22, m5, half_n);
    do_strassen(A21_A11, B11_B12, m6, half_n);
    do_strassen(A12_A22, B21_B22, m7, half_n);
    sum_matrices(half_n, m1, c11);
    sum_matrices(half_n, m4, c11);
    subtract_matrices(half_n, m5, c11);
    sum_matrices(half_n, m7, c11);
    sum_matrices(half_n, m3, c12);
    sum_matrices(half_n, m5, c12);
    sum_matrices(half_n, m2, c21);
    sum_matrices(half_n, m4, c21);
    sum_matrices(half_n, m1, c22);
    subtract_matrices(half_n, m2, c22);
    sum_matrices(half_n, m3, c22);
    sum_matrices(half_n, m6, c22);
    //jprintf_given_2Darray(half_n, c111);
    //printf_given_2Darray(half_n, c121);
    //printf_given_2Darray(half_n, c211);
    //printf_given_2Darray(half_n, c221);

    // reorder the matrix
    rearrange_block_matrix(c11, c12, c21, c22, c, half_n, n);

    // free the matrix
    free(a11);
    free(a12);
    free(a21);
    free(a22);
    free(b11);
    free(b12);
    free(b21);
    free(b22);
}

int main(int argc, char **argv) {
    // Variables
    int N = atoi(argv[1]);
    int print_option = atoi(argv[2]);
    int rank, size, i;
    int tag = 0;
    int t;
    double start_time, end_time;
    float *a, *b, *c;

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
    a = (float*)malloc(sizeof(float) * N * N);
    b = (float*)malloc(sizeof(float) * N * N);
    c = (float*)malloc(sizeof(float) * N * N);


    // initialize arrays only when the rank == 0
    if (rank == 0){
        // init
        init_given_2Darray_randomly(N, a);
        init_given_2Darray_randomly(N, b);
    }


    // do matrix multiplication
    recursive_naive(a, b, c, N);
    //do_strassen(a, b, c, N);

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


