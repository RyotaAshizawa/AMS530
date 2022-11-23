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


int main(int argc, char *argv[]){
    // initialize mpi
    int rank, size, i;
    int max_rank = std::pow(std::floor(std::pow(size, 1./3.)), 3);
    int tag = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request;
    MPI_Status status;

    // set variables
    const int N = std::stoi(argv[1]);
    const double box_size = std::stof(argv[2]);
    const int mpi_tag = 0;

    // initialize box and mpi simulator
    Box box(N, box_size, true);
    box.print_particles();
    MPI_sim mpi_simulator(&box, MPI_COMM_WORLD, &request, &status, size, rank);

    // send and recv number of particles in each mpi box
    if (rank == 0) {
        mpi_simulator.mpi_send_n_particles_eachrank(tag);
    }
    int *n_particles_eachrank = new int[max_rank];
    MPI_Irecv(n_particles_eachrank, max_rank, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    MPI_Finalize();
    return 0;
}

