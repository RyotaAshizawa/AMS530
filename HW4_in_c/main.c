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
    int *map_rank_to_n_particles_in_surrcells = (int *) malloc(sizeof(int) * max_rank);
    double *force_and_id = (double *) malloc(sizeof(double) * 4 * N);
    double **map_rank_to_coords_surrbox = (double **) malloc(sizeof(double *) * max_rank);
    for (i = 0; i < max_rank; i++) {
        map_rank_to_coords_surrbox[i] = (double *) malloc(sizeof(double) * N * 4); // 26 is the max possible n of surr cells
    }
    double *coords_peripheral_box = (double *) malloc(sizeof(double *) * N * 4);


    // Box definition
    for (i = 0; i < N; i++) {
        box[i] = (double *) malloc(sizeof(double) * 8);
    }

    // Prepare data on rank 0
    //// 1. initialize box
    if (rank == 0) {
        init_coords_and_forces(box, true, N, particles_per_side, particle_cellsize);
        dump_particles(box, "./test.xyz", N);
        // MPI definition
        init_map_cell_to_rank(cpus_per_side, map_cell_to_rank, map_rank_to_cell);
        assign_rank_to_cell(box, N, n_particles_eachrank, map_cell_to_rank, cpus_per_side, cell_len_per_cpu, max_rank);
        print_particles(box, N); // check it
    }

    //// 2. Send-Recv number of particles of each cell
    if (rank == 0) {
        mpi_send_n_particles_to_eachrank(n_particles_eachrank, tag, max_rank, MPI_COMM_WORLD, &request);
    }
    if (rank < max_rank) {
        MPI_Irecv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        //printf("Rank:%d, Recv:%d\n", rank, n_particles_eachrank[0]);
    }

    //// 3. Send-recv particle positions of each cell
    if (rank == 0) {
        get_particles_each_rank(box, N, coords_each_rank, max_rank);
        mpi_send_centerbox_particles(n_particles_eachrank, coords_each_rank, max_rank, MPI_COMM_WORLD, &request, tag);
    }
    if (rank < max_rank) {
        MPI_Irecv(coords_center_box, n_particles_eachrank[rank] * 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
    if (rank == 16) {
        //printf("Received particles for the center box for rank %d:\n", rank);
        //print_particles_in_box(coords_center_box, n_particles_eachrank[rank]);
    }

    //// 4. Get info of peripheral cells
    if (rank == 0){
        get_n_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_rank_to_n_surrcells);
        get_rank_of_surrboxes(max_rank, cpus_per_side, map_rank_to_cell, map_cell_to_rank, map_rank_to_ranks_of_surrcells);
        get_tot_particles_in_surrboxes(max_rank, map_rank_to_n_particles_in_surrcells, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells, n_particles_eachrank);
    }

    //// 5. Send-recv n particles of surr cells
    if (rank == 0) {
        mpi_send_n_particles_to_eachrank(map_rank_to_n_particles_in_surrcells, tag, max_rank, MPI_COMM_WORLD, &request);
    }
    if (rank < max_rank) {
        MPI_Irecv(map_rank_to_n_particles_in_surrcells, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        printf("Rank:%d, N of particles in the surrownding cells:%d\n", rank, map_rank_to_n_particles_in_surrcells[rank]);
    }

    //// 6. Send-recv cooridinates of surr cells
    // allocate memory first
    for (int rank = 0; rank < max_rank; rank++) {
        map_rank_to_coords_surrbox[rank] = (double *) malloc(sizeof(double) * map_rank_to_n_particles_in_surrcells[rank] * 4);
    }
    if (rank == 0) {
        map_rank_to_coords_surrcells(max_rank, map_rank_to_n_particles_in_surrcells, map_rank_to_coords_surrbox, coords_each_rank, n_particles_eachrank, map_rank_to_n_surrcells, map_rank_to_ranks_of_surrcells);
        //mpi_send_surrbox_particles(map_rank_to_coords_surrbox, map_rank_to_n_particles_in_surrcells, max_rank, MPI_COMM_WORLD, &request, tag);
    }
    if (rank == 16) {
        printf("%d particles are received for the peripheral box for rank %d:\n", map_rank_to_n_particles_in_surrcells[rank], rank);
        print_particles_in_box(map_rank_to_coords_surrbox[rank], map_rank_to_n_particles_in_surrcells[rank]);
    }
    /**
    if (rank < max_rank) {
        MPI_Irecv(coords_peripheral_box, map_rank_to_n_particles_in_surrcells[rank] * 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        if (rank == 16) {
            printf("%d particles are received for the peripheral box for rank %d:\n", map_rank_to_n_particles_in_surrcells[rank], rank);
            print_particles_in_box(coords_peripheral_box, map_rank_to_n_particles_in_surrcells[rank]);
        }
    }
    **/

    /**


    //// 7. Calculate the force of each main box
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
    for (int i_atom_centerbox = 0; i_atom_centerbox < n_particles_eachrank[rank]; i_atom_centerbox++){
        for (int j_atom_centerbox = 0; j_atom_centerbox < n_particles_eachrank[rank]; j_atom_centerbox++){
            if (i_atom_centerbox != j_atom_centerbox) {
                add_force_between_two_particles_to_vector(&force_and_id[i_atom_centerbox * 4],
                                                          &coords_center_box[j_atom_centerbox * 4],
                                                          &coords_center_box[i_atom_centerbox * 4],
                                                          rank);
            }
        }
    }

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
    //printf("Accumulated (Fx, Fy, Fz) = (%lf, %lf %lf) for id %d for rank %d:\n", force_and_id[0], force_and_id[1], force_and_id[2], (int)force_and_id[3], rank);


    /**
    //MPI_Recv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    int send_test = 0;
    int recv_test[1];
    std::cout << "Rank of the process..:" << rank;
    **/


    MPI_Finalize();
    return 0;


}
