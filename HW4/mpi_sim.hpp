//
// Created by Ryota Ashizawa on 11/22/22.
//
#include <mpi.h>
#include <iostream>
#include "box.hpp"

#ifndef HW4_MPI_SIM_HPP
#define HW4_MPI_SIM_HPP

class MPI_sim {
public:
    MPI_sim(const int N, const int box_size, MPI_Comm comm, MPI_Request *request, MPI_Status *status, const int size, const int rank);
    // member functions
    void print_rank(const int rank);

private:
    // constructor params
    Box box;
    int rank;
    int size;
    int max_rank;
    int N_box;
    MPI_Comm comm;
    MPI_Request *request;
    MPI_Status *status;
    // inner params
    int cpus_per_side;
};

#endif //HW4_MPI_SIM_HPP
