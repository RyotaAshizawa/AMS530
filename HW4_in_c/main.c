#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "particle.h"
#include "box.h"

// mpi functions
void mpi_send_n_particles_to_eachrank(int *n_particles_eachrank, const int tag, const int max_rank, MPI_Comm comm, MPI_Request *request){
    for (int dst_rank = 1; dst_rank < max_rank; dst_rank++){
        MPI_Isend(n_particles_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
    }
}
void mpi_send_centerbox_particles(int *n_particles_eachrank, double **coords_each_rank, const int max_rank, MPI_Comm comm, MPI_Request *request, const int tag){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(coords_each_rank[dst_rank], n_particles_eachrank[dst_rank] * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
    }
}
void copy_coords_surrcells(int rank, double *temp_coords_surrcells, double **coords_each_rank, int *n_particles_eachrank, int *map_rank_to_n_surrcells, int **map_rank_to_ranks_of_surrcells){
    int processed_temp_n_particles = 0;
    for (int i_surr_rank = 0; i_surr_rank < map_rank_to_n_surrcells[rank]; i_surr_rank++) {
        int surr_rank = map_rank_to_ranks_of_surrcells[rank][i_surr_rank];
        if (rank == 13){
            printf("Surr rank for Rank %d: %d \n", rank, surr_rank); //OK
            printf("N particles in surr rank %d: %d\n", surr_rank, n_particles_eachrank[surr_rank]);

        }
        for (int j_particle = 0; j_particle < n_particles_eachrank[surr_rank]; j_particle++){
            temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 0] = coords_each_rank[surr_rank][j_particle * 4 + 0];
            temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 1] = coords_each_rank[surr_rank][j_particle * 4 + 1];
            temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 2] = coords_each_rank[surr_rank][j_particle * 4 + 2];
            temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 3] = coords_each_rank[surr_rank][j_particle * 4 + 3];
            if (rank == 13){
                printf("1:processed particle for Rank %d: %d, (x, y, z, id) = (%f, %f, %f, %f)\n", rank, processed_temp_n_particles, coords_each_rank[surr_rank][j_particle * 4 + 0], coords_each_rank[surr_rank][j_particle * 4 + 1], coords_each_rank[surr_rank][j_particle * 4 + 2], coords_each_rank[surr_rank][j_particle * 4 + 3]); //OK
                printf("2:processed particle for Rank %d: %d, (x, y, z, id) = (%f, %f, %f, %f)\n", rank, processed_temp_n_particles, temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 0], temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 1], temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 2], temp_coords_surrcells[processed_temp_n_particles * 4 + j_particle * 4 + 3]); //OK
            }
            processed_temp_n_particles += 1;
        }
    }
}
void mpi_send_surrbox_particles(int *n_particles_eachrank, int *map_rank_to_n_particles_in_surrcells, int *map_rank_to_n_surrcells, int **map_rank_to_ranks_of_surrcells, double **coords_each_rank, int max_rank, int N, MPI_Comm comm, MPI_Request *request, const int tag){
    //double temp_coords_surrcells[4 * N];
    double *temp_coords_surrcells = (double *) malloc(sizeof(double) * 4 * N);
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        copy_coords_surrcells(dst_rank, temp_coords_surrcells, coords_each_rank, n_particles_eachrank, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells);
        //MPI_Isend(temp_coords_surrcells, map_rank_to_n_particles_in_surrcells[dst_rank] * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
        if (dst_rank == 13){
            //printf("N in surr cells for Rank %d: %d \n",dst_rank, map_rank_to_n_particles_in_surrcells[dst_rank]); // map_rank_to_n_particles_in_surrcells is OK
            //printf("N  surr cells for Rank %d: %d \n",dst_rank, map_rank_to_n_surrcells[dst_rank]); // map_rank_to_n_surrcells is OK
            print_particles_in_box(temp_coords_surrcells, map_rank_to_n_particles_in_surrcells[dst_rank]);
        }
        //free(temp_coords_surrcells);
    }
}

int main(int argc, char **argv) {
    int rank, size, i;
    int tag = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request;
    MPI_Status status;

    // box variables
    int N = atoi(argv[1]);
    double box_size = atoi(argv[2]);
    int particles_per_side = ceil(pow(N, 1. / 3.));
    double particle_cellsize = box_size / particles_per_side;

    // mpi variables
    int cpus_per_side = floor(pow(size, 1. / 3.));
    double cell_len_per_cpu = box_size / cpus_per_side;
    int max_rank = cpus_per_side * cpus_per_side * cpus_per_side;

    // memory allocations
    int *n_particles_eachrank = (int *) malloc(sizeof(int) * max_rank);
    int *map_rank_to_cell = (int *) malloc(sizeof(int) * 3 * max_rank);
    int *map_cell_to_rank = (int *) malloc(sizeof(int) * max_rank);
    double **box = (double **) malloc(sizeof(double *) * N); //first assign box memory
    double **coords_each_rank = (double **) malloc(sizeof(double *) * max_rank);
    for (i = 0; i < max_rank; i++) {
        coords_each_rank[i] = (double *) malloc(sizeof(double) * 4 * N);
    }
    double *coords_center_box = (double *) malloc(sizeof(double) * 4 * N);
    int *map_rank_to_n_surrcells = (int *) malloc(sizeof(int) * max_rank);
    int **map_rank_to_ranks_of_surrcells = (int **) malloc(sizeof(int *) * max_rank);
    for (i = 0; i < max_rank; i++) {
        map_rank_to_ranks_of_surrcells[i] = (int *) malloc(sizeof(int) * 26); // 26 is the max possible n of surr cells
    }
    double *particles_surr_box = (double *) malloc(sizeof(double) * 4 * N);
    int *map_rank_to_n_particles_in_surrcells = (int *)malloc(sizeof(int) * max_rank);
    double *temp_coords_surrcells = (double *)malloc(sizeof(double) * 4 * N);

    // Box definition
    for (i = 0; i < N; i++) {
        box[i] = (double *) malloc(sizeof(double) * 8);
    }

    //// 1. initialize box
    if (rank == 0) {
        init_coords_and_forces(box, true, N, particles_per_side, particle_cellsize);
        dump_particles(box, "./test.xyz", N);
        // MPI definition
        init_map_cell_to_rank(cpus_per_side, map_cell_to_rank, map_rank_to_cell);
        assign_rank_to_cell(box, N, n_particles_eachrank, map_cell_to_rank, cpus_per_side, cell_len_per_cpu, max_rank);
        print_particles(box, N); // check it
        mpi_send_n_particles_to_eachrank(n_particles_eachrank, tag, max_rank, MPI_COMM_WORLD, &request);
    }

    //// 2. Send number of particles of each cell
    if (rank < max_rank && rank != 0) {
        MPI_Irecv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        //printf("Rank:%d, Recv:%d\n", rank, n_particles_eachrank[0]);
    }

    //// 3. Copy particle positions of each cell
    if (rank == 0) {
        get_particles_each_rank(box, N, coords_each_rank, max_rank);
        //for (i = 0; i < max_rank; i++){
        //    print_particles_special_rank(coords_each_rank, n_particles_eachrank, i);
        //}
    }

    //// 4. Send and recv the position of particles in the center box
    // receive the data for the centered box
    if (rank == 0) {
        mpi_send_centerbox_particles(n_particles_eachrank, coords_each_rank, max_rank, MPI_COMM_WORLD, &request, tag);
    }
    if (rank < max_rank) {
        MPI_Irecv(coords_center_box, n_particles_eachrank[rank] * 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
    //printf("Received particles for the center box for rank %d:\n", rank);
    //print_particles_in_box(coords_center_box, n_particles_eachrank[rank]);


    //// 5. Send and recv the position of particles in the peripheral boxes
    if (rank == 0) {
        get_n_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_rank_to_n_surrcells);
        get_rank_of_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_cell_to_rank, map_rank_to_ranks_of_surrcells);
        get_tot_particles_in_surrboxes(max_rank, map_rank_to_n_particles_in_surrcells, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells, n_particles_eachrank);
        mpi_send_surrbox_particles(n_particles_eachrank, map_rank_to_n_particles_in_surrcells, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells, coords_each_rank, max_rank, N, MPI_COMM_WORLD, &request, tag);
    }


    // caluclate force



    /**
    //MPI_Recv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    int send_test = 0;
    int recv_test[1];
    std::cout << "Rank of the process..:" << rank;
    **/


    MPI_Finalize();
    return 0;


}
