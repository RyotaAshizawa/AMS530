//
// Created by Ryota Ashizawa on 11/22/22.
//

#include <fstream>
#include <string>
#include "particle.hpp"
#include "box.hpp"
#include "mpi.h"


int main(int argc, char *argv[]){
    // initialize mpi
    int rank, size, i;
    int max_rank = std::pow(std::floor(std::pow(size, 1./3.)), 3);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request *request;
    MPI_Status *status;

    // set variables
    const int N = std::stoi(argv[1]);
    const double box_size = std::stof(argv[2]);
    const int mpi_tag_n_per_proc = 0;
    const int mpi_tag_main_box = 1;
    const int mpi_tag_nsurrbox_per_proc = 2;
    const int mpi_tag_rank_to_surrranks = 3;
    const int mpi_tag_surr_box = 4;

    // create box on rank 0 and send the pointer
    Box box(N, box_size, true);
    if (rank == 0){
        box.mpi_send(MPI_COMM_WORLD, request, 0, size);
    }
    int recv;
    MPI_Irecv(&recv, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, request);
    MPI_Wait(request, status);
    //MPI_simulator.print_rank(rank);

    MPI_Finalize();
    return 0;
}

