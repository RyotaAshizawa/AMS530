//
// Created by Ryota Ashizawa on 11/23/22.
//

#include "box.h"

// member functions
void print_particles(double **box, const int N) {
    for (int i = 0; i < N; i++) {
        print_particle(box[i]);
    }
}
void dump_particles(double **box, char filepath[], const int N){
    FILE *fp = fopen(filepath, "w");
    fprintf(fp, "%d\n", N);
    fprintf(fp, "Initial coordinate\n");
    for (int i = 0; i < N; i++) {
        fprintf(fp, "He\t%f\t%f\t%f\n", box[i][0], box[i][1], box[i][2]);
    }
    fclose(fp);
}

void init_coords_and_forces(double **box, bool use_rand, const int N, const int particles_per_side, const double particle_cellsize) {
    // assume each particle locates at the center of each cell
    bool flag;
    double x, y, z;
    int particle_count = 0;
    for (int i = 0; i < particles_per_side; i++) {
        for (int j = 0; j < particles_per_side; j++) {
            for (int k = 0; k < particles_per_side; k++) {
                // set position
                x = particle_cellsize * (i + 0.5);
                y = particle_cellsize * (j + 0.5);
                z = particle_cellsize * (k + 0.5);
                // add random perturbation
                if (use_rand){
                    x = x + particle_cellsize * ((double)rand() / RAND_MAX * 0.2);
                    y = y + particle_cellsize * ((double)rand() / RAND_MAX * 0.2);
                    z = z + particle_cellsize * ((double)rand() / RAND_MAX * 0.2);
                }
                set_coordinate(box[particle_count], x, y, z, 0, 0, 0, 0, particle_count);
                // Break the loop if the positions of all particles are set.
                particle_count = particle_count + 1;
                if (particle_count == N) {
                    flag = true;
                    break;
                }
            }
            if (flag){break;}
        }
        if (flag){break;}
    }
}


//mpi functions
void init_map_cell_to_rank(int cpus_per_side, int map_cell_to_rank[cpus_per_side][cpus_per_side][cpus_per_side], int *map_rank_to_cell) {
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
void assign_rank_to_box(double **box, const int N, int *n_particles_eachrank, int ***map_cell_to_rank, const int cell_len_per_cpu, const int max_rank) {
    // array definition and initialize
    for (int i = 0; i < max_rank; i++) {
        n_particles_eachrank[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        int cellno_in_x = floor(get_x(box[i]) / cell_len_per_cpu);
        int cellno_in_y = floor(get_y(box[i]) / cell_len_per_cpu);
        int cellno_in_z = floor(get_z(box[i]) / cell_len_per_cpu);
        int rank = map_cell_to_rank[cellno_in_x][cellno_in_y][cellno_in_z];
        set_rank(box[i], rank);
        n_particles_eachrank[rank]++;
    }
}
