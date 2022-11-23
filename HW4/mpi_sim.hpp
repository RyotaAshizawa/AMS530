//
// Created by Ryota Ashizawa on 11/22/22.
//
#include <mpi.h>
#include <iostream>
#include <cmath>
#include "box.hpp"

#ifndef HW4_MPI_SIM_HPP
#define HW4_MPI_SIM_HPP

class MPI_sim {
public:
    MPI_sim(Box *box, MPI_Comm comm, MPI_Request *request, MPI_Status *status, const int size, const int rank);
    ~MPI_sim();
    // member functions
    void print_rank(const int rank);
    void mpi_send_n_particles_eachrank(const int tag);

private:
    // constructor params
    Box *box;
    int rank;
    int size;
    int max_rank;
    //box
    int N;
    int box_size;
    //mpi
    MPI_Comm comm;
    MPI_Request *request;
    MPI_Status *status;
    // inner params
    int cpus_per_side;
    int cell_len_per_cpu;
    int *n_particles_eachrank;
    int *map_rank_to_cell;
    int ***map_cell_to_rank;


    //private funcs
    void init_map_cell_to_rank();
    void assign_rank_to_box();
};

#endif //HW4_MPI_SIM_HPP
