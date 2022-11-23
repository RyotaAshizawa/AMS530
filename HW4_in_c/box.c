//
// Created by Ryota Ashizawa on 11/23/22.
//

#include "box.h"

// member functions
Particle get_particle(double *box, const int i){
    return box->particles[i];
}
void print_particles(double *box) {
    for (int i = 0; i < box -> N; i++) {
        print_particle(&(box->particles)[i]);
    }
}
void set_box_size(double *box, const double box_size){
    box -> box_size = box_size;
}
void set_n_in_box(double *box, const int N){
    box -> N = N;
}
void dump_particles(double *box, char filepath[]){
    FILE *fp = fopen(filepath, "w");
    fprintf(fp, "%d\n", box -> N);
    fprintf(fp, "Initial coordinate\n");
    for (int i = 0; i < box -> N; i++) {
        fprintf(fp, "He\t%f\t%f\t%f\n", (box->particles)[i].x, (box->particles)[i].y, (box->particles)[i].z);
    }
    fclose(fp);
}

void init_coords_and_forces(double *box, bool use_rand, const int particles_per_side, const double particle_cellsize) {
    /** assume each particle locates at the center of each cell **/
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
                set_coordinate(&(box -> particles[particle_count]), x, y, z, 0, 0, 0, 0, particle_count);
                // Break the loop if the positions of all particles are set.
                particle_count++;
                if (particle_count == box -> N) {break;}
            }
        }
    }
}


//mpi functions
void init_map_cell_to_rank(double *box, int cpus_per_side, int ***map_cell_to_rank, int *map_rank_to_cell) {
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
void assign_rank_to_box(double *box, int *n_particles_eachrank, int ***map_cell_to_rank, const int cell_len_per_cpu, const int max_rank) {
    // array definition and initialize
    for (int i = 0; i < max_rank; i++) {
        n_particles_eachrank[i] = 0;
    }
    for (int i = 0; i < box->N; i++) {
        int cellno_in_x = floor(get_x(&box->particles[i]) / cell_len_per_cpu);
        int cellno_in_y = floor(get_y(&box->particles[i]) / cell_len_per_cpu);
        int cellno_in_z = floor(get_z(&box->particles[i]) / cell_len_per_cpu);
        int rank = map_cell_to_rank[cellno_in_x][cellno_in_y][cellno_in_z];
        set_rank(&(box->particles)[i], rank);
        n_particles_eachrank[rank]++;
    }
}
