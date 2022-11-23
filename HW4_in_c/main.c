#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "particle.h"
#include "box.h"

// mpi functions
void mpi_send_n_particles_to_eachrank(int *n_particles_eachrank, const int tag, const int max_rank, MPI_Comm comm, MPI_Request *request){
    for (int dst_rank = 1; dst_rank < max_rank; dst_rank++){
        MPI_Isend(n_particles_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
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
    int particles_per_side = ceil(pow(N, 1./3.));
    double particle_cellsize = box_size / particles_per_side;

    // mpi variables
    int cpus_per_side = floor(pow(size, 1./3.));
    double cell_len_per_cpu = box_size / cpus_per_side;
    int max_rank = cpus_per_side * cpus_per_side * cpus_per_side;

    // mapping between cell and rank
    int n_particles_eachrank[max_rank];
    int map_rank_to_cell[max_rank];
    int map_cell_to_rank[cpus_per_side][cpus_per_side][cpus_per_side];


    // Box definition
    if (rank == 0) {
        // Box setting
        Box *box = (Box *)malloc(sizeof(Box));
        set_box_size(box, box_size);
        set_n_in_box(box, N);
        init_coords_and_forces(box, true, particles_per_side, particle_cellsize);
        dump_particles(box, "./test.xyz");
        print_particles(box);

        // assign mpi mapping
        init_map_cell_to_rank(box, cpus_per_side, map_cell_to_rank, map_rank_to_cell);
        assign_rank_to_box(box, n_particles_eachrank, cpus_per_side, map_cell_to_rank, cell_len_per_cpu, max_rank);

        /** printing **/
        for (int i = 0; i < max_rank; i++){
            printf("%d\n", n_particles_eachrank[i]);
        }
        mpi_send_n_particles_to_eachrank(n_particles_eachrank, tag, max_rank, MPI_COMM_WORLD, &request);
        //int i = 1;
        //MPI_Isend(n_particles_eachrank, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &request);
    }

    // receive number of particles
    if (rank < max_rank)  {
        MPI_Irecv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        printf("Rank:%d, Recv:%d\n", rank, n_particles_eachrank[0]);
    }


    /**
    //MPI_Recv(n_particles_eachrank, max_rank, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    int send_test = 0;
    int recv_test[1];
    std::cout << "Rank of the process..:" << rank;
    **/


    MPI_Finalize();
    return 0;


}
