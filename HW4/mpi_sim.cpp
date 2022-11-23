//
// Created by Ryota Ashizawa on 11/22/22.
//

#include "mpi_sim.hpp"

MPI_sim::MPI_sim(Box *box, MPI_Comm comm, MPI_Request *request, MPI_Status *status, const int size, const int rank){
    // arg
    this -> size = size;
    this -> rank = rank;
    this -> comm = comm;
    this -> status = status;
    this -> request = request;
    std::cout << "test";
    // others
    cpus_per_side = std::floor(std::pow(size, 1./3.));
    max_rank = cpus_per_side * cpus_per_side * cpus_per_side;
    cell_len_per_cpu = box_size / double(cpus_per_side);
    N = box -> get_n_particles();
    box_size = box -> get_box_size();
    // initalize box
    // new
    n_particles_eachrank = new int[max_rank];
    map_rank_to_cell = new int[3 * max_rank];
    map_cell_to_rank = new int**[cpus_per_side];
    for (int i = 0; i < cpus_per_side; i++) {
        map_cell_to_rank[i] = new int*[cpus_per_side];
        for (int j = 0; j < cpus_per_side; j++) {
            map_cell_to_rank[i][j] = new int[cpus_per_side];
        }
    }

    // initialize values
    init_map_cell_to_rank();
    assign_rank_to_box();
}

MPI_sim::~MPI_sim(){
    delete[] n_particles_eachrank;
    delete[] map_rank_to_cell;
}

void MPI_sim::print_rank(const int rank){
    std::cout << rank << "dayo\n";
}

// private mpi functions
void MPI_sim::init_map_cell_to_rank() {
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
void MPI_sim::assign_rank_to_box() {
    // array definition and initialize
    for (int i = 0; i < max_rank; i++){
        n_particles_eachrank[i] = 0;
    }
    for (int i = 0; i < N; i++){
        int cellno_in_x = floor(box->get_particle(i).get_x() / cell_len_per_cpu);
        int cellno_in_y = floor(box->get_particle(i).get_y() / cell_len_per_cpu);
        int cellno_in_z = floor(box->get_particle(i).get_z() / cell_len_per_cpu);
        int rank = map_cell_to_rank[cellno_in_x][cellno_in_y][cellno_in_z];
        box->get_particle(i).set_rank(rank);
        n_particles_eachrank[rank]++;
    }
}

// public mpi send functions
void MPI_sim::mpi_send_n_particles_eachrank(const int tag){
    for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
        MPI_Isend(n_particles_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
    }
}
