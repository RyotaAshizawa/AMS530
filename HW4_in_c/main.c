#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "particle.h"
#include "box.h"
#include "force.h"

// mpi functions
void mpi_send_n_particles_to_eachrank(int *n_particles_eachrank, const int tag, const int max_rank, MPI_Comm comm, MPI_Request *request){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(n_particles_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
    }
}
void mpi_send_centerbox_particles(int *n_particles_eachrank, double **coords_each_rank, const int max_rank, MPI_Comm comm, MPI_Request *request, const int tag){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(coords_each_rank[dst_rank], n_particles_eachrank[dst_rank] * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
    }
}
void mpi_send_surrbox_particles(double **map_rank_to_coords_surrbox, int *map_rank_to_n_particles_in_surrcells, int max_rank, MPI_Comm comm, MPI_Request *request, const int tag){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(map_rank_to_coords_surrbox[dst_rank], map_rank_to_n_particles_in_surrcells[dst_rank] * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
    }
}
//mpi_send_surrbox_particles(map_rank_to_coords_surrbox, map_rank_to_n_particles_in_surrcells, max_rank, MPI_COMM_WORLD, &request, tag);


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
    double cell_len_per_cpu = (double)box_size / cpus_per_side;
    int max_rank = cpus_per_side * cpus_per_side * cpus_per_side;

    // memory allocations
    int *n_particles_eachrank_send = (int *) malloc(sizeof(int) * max_rank);
    int *n_particles_eachrank = (int *) malloc(sizeof(int) * max_rank);
    int *map_rank_to_cell = (int *) malloc(sizeof(int) * 3 * max_rank);
    int *map_cell_to_rank = (int *) malloc(sizeof(int) * max_rank);
    //// Box definition
    double **box = (double **) malloc(sizeof(double *) * N); //first assign box memory
    for (i = 0; i < N; i++) {
        box[i] = (double *) malloc(sizeof(double) * 8);
    }
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
    int *map_rank_to_n_particles_in_surrcells_send = (int *) malloc(sizeof(int) * max_rank);
    int *map_rank_to_n_particles_in_surrcells = (int *) malloc(sizeof(int) * max_rank);
    double *force_and_id = (double *) malloc(sizeof(double) * 4 * N);
    double **map_rank_to_coords_surrbox = (double **) malloc(sizeof(double *) * max_rank);
    for (i = 0; i < max_rank; i++) {
        map_rank_to_coords_surrbox[i] = (double *) malloc(sizeof(double) * N * 4); // 26 is the max possible n of surr cells
    }
    double *coords_peripheral_box = (double *) malloc(sizeof(double *) * N * 4);

    // debug option
    int rank_interest = 26;

    // Prepare data on rank 0
    //// 1. initialize box
    if (rank == 0) {
        init_coords_and_forces(box, false, N, particles_per_side, particle_cellsize);
        dump_particles(box, "./test.xyz", N);
        // MPI definition
        init_map_cell_to_rank(cpus_per_side, map_cell_to_rank, map_rank_to_cell, false);
        assign_rank_to_cell(box, N, n_particles_eachrank_send, map_cell_to_rank, cpus_per_side, cell_len_per_cpu, max_rank, false);
        //print_particles(box, N);
    }

    //// 2. Send-Recv number of particles of each cell
    if (rank == 0) {
        mpi_send_n_particles_to_eachrank(n_particles_eachrank_send, tag, max_rank, MPI_COMM_WORLD, &request);
    }
    if (rank < max_rank) {
        MPI_Irecv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
    //printf("Rank:%d, Recv:%d\n", rank, n_particles_eachrank[rank]);

    //// 3. Send-recv particle positions of each cell
    if (rank == 0) {
        get_particles_each_rank(box, N, coords_each_rank, max_rank, false, rank_interest);
        mpi_send_centerbox_particles(n_particles_eachrank, coords_each_rank, max_rank, MPI_COMM_WORLD, &request, tag);
    }
    if (rank < max_rank) {
        MPI_Irecv(coords_center_box, n_particles_eachrank[rank] * 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
    if (rank == rank_interest) {
        //printf("Received particles for the center box for rank %d:\n", rank);
        //print_particles_in_box(coords_center_box, n_particles_eachrank[rank]);
    }

    //// 4. Get info of peripheral cells
    if (rank == 0){
        get_n_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_rank_to_n_surrcells);
        get_rank_of_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_cell_to_rank, map_rank_to_ranks_of_surrcells);
        get_tot_particles_in_surrboxes(max_rank, map_rank_to_n_particles_in_surrcells_send, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells, n_particles_eachrank);
    }

    //// 5. Send-recv n particles of surr cells
    if (rank == 0) {
        mpi_send_n_particles_to_eachrank(map_rank_to_n_particles_in_surrcells_send, tag, max_rank, MPI_COMM_WORLD, &request);
    }
    if (rank < max_rank) {
        MPI_Irecv(map_rank_to_n_particles_in_surrcells, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
    //printf("Rank:%d, N of particles in the surrownding cells:%d\n", rank, map_rank_to_n_particles_in_surrcells[rank]);

    /**
    //// 6. Send-recv cooridinates of surr cells
    if (rank == 0) {
        map_rank_to_coords_surrcells(max_rank, map_rank_to_n_particles_in_surrcells, map_rank_to_coords_surrbox, coords_each_rank, n_particles_eachrank, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells, false, rank_interest);
        mpi_send_surrbox_particles(map_rank_to_coords_surrbox, map_rank_to_n_particles_in_surrcells, max_rank, MPI_COMM_WORLD, &request, tag);
        print_particles_in_box(map_rank_to_coords_surrbox[rank_interest], map_rank_to_n_particles_in_surrcells[rank_interest]);
    }
    if (rank < max_rank) {
        MPI_Irecv(coords_peripheral_box, map_rank_to_n_particles_in_surrcells[rank] * 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        if (rank == rank_interest) {
            printf("%d particles are received for the peripheral box for rank %d:\n", map_rank_to_n_particles_in_surrcells[rank], rank);
            print_particles_in_box(coords_peripheral_box, map_rank_to_n_particles_in_surrcells[rank]);
        }
    }


    //// 7. force debugger
    if (rank == 0) {
        double *coords1 = (double *) malloc(sizeof(double) * 4);
        double *coords2 = (double *) malloc(sizeof(double) * 4);
        double *force_and_id = (double *) malloc(sizeof(double) * 4);
        coords1[0] = 0;
        coords1[1] = 0;
        coords1[2] = 0;
        coords1[3] = 0;
        coords2[0] = 1;
        coords2[1] = 1;
        coords2[2] = 1;
        coords2[3] = 1;
        add_force_between_two_particles_to_vector(force_and_id, coords2, coords1);
        printf("Accumulated (Fx, Fy, Fz) = (%f, %f %f) for id %d\n", force_and_id[0], force_and_id[1], force_and_id[2], (int)force_and_id[3]);
        add_force_between_two_particles_to_vector(force_and_id, coords2, coords1);
        printf("Accumulated (Fx, Fy, Fz) = (%f, %f %f) for id %d\n", force_and_id[0], force_and_id[1], force_and_id[2], (int)force_and_id[3]);
    }

    //// 8. Calculate force inside the main box
    for (int i = 0; i < n_particles_eachrank[rank]; i++){
        for (int j = 0; j < n_particles_eachrank[rank]; j++){
            if (i != j) { // ignore same particle
                add_force_between_two_particles_to_vector(&force_and_id[i * 4],
                                                          &coords_center_box[j * 4],
                                                          &coords_center_box[i * 4],
                                                          rank);
            }
        }
    }

    if (rank == rank_interest){
        for (int i = 0; i < n_particles_eachrank[rank]; i++) {
            printf("Force   (Fx, Fy, Fz) = (%lf, %lf %lf) for id %d for %d:\n", force_and_id[i * 4 + 0], force_and_id[i * 4 + 1], force_and_id[i * 4 + 2], (int) force_and_id[i * 4 + 3], rank);
        }
    }
    **/

    /**
    //// 9. Calculate force betweeen the main box and surrownding boxes
    for (int i_atom_centerbox = 0; i_atom_centerbox < n_particles_eachrank[rank]; i_atom_centerbox++) {
        for (int j_atom_surrboxes = 0; j_atom_surrboxes < map_rank_to_n_particles_in_surrcells[rank]; j_atom_surrboxes++) {
            add_force_between_two_particles_to_vector(&force_and_id[i_atom_centerbox * 4],
                                                      &coords_peripheral_box[j_atom_surrboxes * 4],
                                                      &coords_center_box[i_atom_centerbox * 4],
                                                      rank);
        }
    }
    **/

    // free
    free(n_particles_eachrank);
    free(map_rank_to_cell);
    for (i = 0; i < N; i++) {
        free(box[i]);
    }
    free(box);
    free(map_cell_to_rank);
    for (i = 0; i < max_rank; i++) {
        free(coords_each_rank[i]);
    }
    free(coords_each_rank);
    free(coords_center_box);
    for (i = 0; i < max_rank; i++) {
        free(map_rank_to_ranks_of_surrcells[i]); // 26 is the max possible n of surr cells
    }
    free(map_rank_to_ranks_of_surrcells);
    free(map_rank_to_n_surrcells);
    free(particles_surr_box);
    free(map_rank_to_n_particles_in_surrcells);
    free(force_and_id);
    for (i = 0; i < max_rank; i++) {
        free(map_rank_to_coords_surrbox[i]); // 26 is the max possible n of surr cells
    }
    free(map_rank_to_coords_surrbox);
    free(coords_peripheral_box);

    MPI_Finalize();
    return 0;


}
