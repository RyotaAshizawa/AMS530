//
// Created by Ryota Ashizawa on 11/22/22.
//

#include <fstream>
#include <string>
#include <cmath>
#include <mpi.h>
#include <iostream>
#include "particle.hpp"
#include "box.hpp"
#include "mpi_sim.hpp"



void init_map_cell_to_rank(int cpus_per_side, int *map_rank_to_cell, int ***map_cell_to_rank) {
    int rank = 0;
    for (int i = 0; i < cpus_per_side; i++) {
        for (int j = 0; j < cpus_per_side; j++) {
            for (int k = 0; k < cpus_per_side; k++) {
                map_cell_to_rank[i][j][k] = rank;
                map_rank_to_cell[rank * 3 + 0] = i;
                map_rank_to_cell[rank * 3 + 1] = j;
                map_rank_to_cell[rank * 3 + 2] = k;
                rank++;
            }
        }
    }
}

void assign_rank_to_box(Box *box, int max_rank, int *n_particles_eachrank, int ***map_cell_to_rank, double cell_len_per_cpu) {
    // array definition and initialize
    for (int i = 0; i < max_rank; i++){
        n_particles_eachrank[i] = 0;
    }
    for (int i = 0; i < box->get_n_particles(); i++){
        int cellno_in_x = floor(box->get_particle(i).get_x() / cell_len_per_cpu);
        int cellno_in_y = floor(box->get_particle(i).get_y() / cell_len_per_cpu);
        int cellno_in_z = floor(box->get_particle(i).get_z() / cell_len_per_cpu);
        int rank = map_cell_to_rank[cellno_in_x][cellno_in_y][cellno_in_z];
        box->get_particle(i).set_rank(rank);
        n_particles_eachrank[rank]++;
    }
}

// public mpi send functions
void mpi_send_n_particles_eachrank(int *n_particles_eachrank, MPI_Comm comm, MPI_Request *request, int max_rank, const int tag){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(n_particles_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
    }
}


int main(int argc, char *argv[]){
    // initialize mpi
    int rank, size, i;
    int tag = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request;
    MPI_Status status;

    // box variables
    const int N = std::stoi(argv[1]);
    const double box_size = std::stof(argv[2]);

    // mpi variables
    int cpus_per_side = std::floor(std::pow(size, 1./3.));
    double cell_len_per_cpu = box_size / double(cpus_per_side);
    int max_rank = cpus_per_side * cpus_per_side * cpus_per_side;


    //Box *box = new Box(N, box_size, false);
    Box *box = new Box(N, box_size, false);

    // allocate memory
    int *n_particles_eachrank = new int[max_rank];
    int *map_rank_to_cell = new int[3 * max_rank];
    int ***map_cell_to_rank = new int**[cpus_per_side];
    for (int i = 0; i < cpus_per_side; i++) {
        map_cell_to_rank[i] = new int*[cpus_per_side];
        for (int j = 0; j < cpus_per_side; j++) {
            map_cell_to_rank[i][j] = new int[cpus_per_side];
        }
    }

    // initialize box and mpi things
    init_map_cell_to_rank(cpus_per_side, map_rank_to_cell, map_cell_to_rank);
    assign_rank_to_box(box, max_rank, n_particles_eachrank, map_cell_to_rank, cell_len_per_cpu);
    if (rank == 0) {
        box->print_particles();
        for (int i = 0; i< max_rank; i++){
            std::cout << n_particles_eachrank[i] << std::endl;
        }
    }

    // send and recv number of particles in each mpi box
    //if (rank == 0) {
    //   mpi_send_n_particles_eachrank(n_particles_eachrank, MPI_COMM_WORLD, &request, max_rank, tag);
    //}

    //MPI_Recv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    //MPI_Wait(&request, &status);


    MPI_Finalize();
    return 0;
}

