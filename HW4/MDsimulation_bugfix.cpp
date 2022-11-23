//
// Created by Ryota Ashizawa on 11/22/22.
//

#include <fstream>
#include <string>
#include "particle.hpp"
#include "box.hpp"
#include <cmath>
#include <mpi.h>
#include <iostream>


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
    int tag = 0;
    if (rank == 0){
        int item[1] = {10};
        Box box(N, box_size, true);
        box.mpi_send(MPI_COMM_WORLD, request, tag, size);
        //for (int rank = 0; rank < size; rank++) {
        //    MPI_Send(item, 1, MPI_INT, rank, tag, MPI_COMM_WORLD);
        //}
    }
    int recv[1];
    MPI_Recv(recv, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
    std::cout << "Rank:" << rank << ", " << "recv:"  << recv[0] << std::endl;
    MPI_Finalize();
    return 0;
}

