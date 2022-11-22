//
// Created by Ryota Ashizawa on 11/22/22.
//

#include "mpi_sim.hpp"

MPI_sim::MPI_sim(const int N, const int box_size, MPI_Comm comm, MPI_Request *request, MPI_Status *status, const int size, const int rank){
    // arg
    this -> size = size;
    this -> rank = rank;
    this -> comm = comm;
    this -> status = status;
    this -> request = request;
    std::cout << "test";
    // others
    //N_box = box->get_n_particles();
    // initalize box
    if (rank == 0) {
        this -> box = Box(N, box_size, true);
    }
}

void MPI_sim::print_rank(const int rank){
    std::cout << rank << "dayo\n";
}