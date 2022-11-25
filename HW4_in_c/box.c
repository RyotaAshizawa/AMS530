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
void init_map_cell_to_rank(int cpus_per_side, int *map_cell_to_rank, int *map_rank_to_cell) {
    int rank = 0;
    for (int i = 0; i < cpus_per_side; i++) {
        for (int j = 0; j < cpus_per_side; j++) {
            for (int k = 0; k < cpus_per_side; k++) {
                map_cell_to_rank[i * cpus_per_side * cpus_per_side + j * cpus_per_side + k] = rank;
                map_rank_to_cell[rank * 3 + 0] = i;
                map_rank_to_cell[rank * 3 + 1] = j;
                map_rank_to_cell[rank * 3 + 2] = k;
                rank++;
            }
        }
    }
}
void assign_rank_to_cell(double **box, const int N, int *n_particles_eachrank, int *map_cell_to_rank, const int cpu_per_side, const int cell_len_per_cpu, const int max_rank) {
    // array definition and initialize
    for (int i = 0; i < max_rank; i++) {
        n_particles_eachrank[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        int cellno_in_x = floor(get_x(box[i]) / cell_len_per_cpu);
        int cellno_in_y = floor(get_y(box[i]) / cell_len_per_cpu);
        int cellno_in_z = floor(get_z(box[i]) / cell_len_per_cpu);
        int rank = map_cell_to_rank[cellno_in_x * cpu_per_side * cpu_per_side + cellno_in_y * cpu_per_side + cellno_in_z];
        set_rank(box[i], rank);
        n_particles_eachrank[rank]++;
    }
}
void get_particles_each_rank(double **box, const int N, double **coords_each_rank, const int max_rank){
    // array definition and initialize
    int n_assigned_p_each_rank[max_rank];
    for (int i = 0; i < max_rank; i++){
        n_assigned_p_each_rank[i] = 0;
    }
    for (int i = 0; i < N; i++){
        int rank = get_rank(box[i]);
        coords_each_rank[rank][n_assigned_p_each_rank[rank] * 4 + 0] = get_x(box[i]);
        coords_each_rank[rank][n_assigned_p_each_rank[rank] * 4 + 1] = get_y(box[i]);
        coords_each_rank[rank][n_assigned_p_each_rank[rank] * 4 + 2] = get_z(box[i]);
        coords_each_rank[rank][n_assigned_p_each_rank[rank] * 4 + 3] = get_id(box[i]);
        n_assigned_p_each_rank[rank]++;
    }
}
void print_particles_special_rank(double **coords_each_rank, int *n_particles_eachrank, const int rank){
    printf("Coords for rank %d:\n", rank);
    for (int i = 0; i < n_particles_eachrank[rank]; i++){
        double *coords_this_rank = coords_each_rank[rank];
        printf("(x, y, z, id, rank) = (%.2f, %.2f, %.2f, %d, %d)\n", coords_this_rank[4 * i + 0], coords_this_rank[4 * i + 1], coords_this_rank[4 * i + 2], (int)coords_this_rank[4 * i + 3], rank);
    }
}
void print_particles_in_box(double *particles_coords, int particles_in_box){
    for (int i = 0; i < particles_in_box; i++){
        printf("%d, (x, y, z, id) = (%.2f, %.2f, %.2f, %d)\n", i, particles_coords[4 * i + 0], particles_coords[4 * i + 1], particles_coords[4 * i + 2], (int)particles_coords[4 * i + 3]);
    }
}
int get_startidx_of_surrboxes(const int centered_box_idx){
    int start_idx = centered_box_idx - 1;
    if (start_idx < 0) {
        start_idx = 0;
    }
    return start_idx;
}
int get_endidx_of_surrboxes(const int cpus_per_side, const int centered_box_idx){
    int end_idx = centered_box_idx + 1;
    if (cpus_per_side <= end_idx) {
        end_idx = cpus_per_side - 1;
    }
    return end_idx;
}
void get_n_surrboxes(const int max_rank, const int cpus_per_side, int *map_rank_to_cell, int *map_rank_to_n_surrcells){
    for (int rank = 0; rank < max_rank; rank++) {
        int i_centered_box = map_rank_to_cell[rank * 3 + 0];
        int j_centered_box = map_rank_to_cell[rank * 3 + 1];
        int k_centered_box = map_rank_to_cell[rank * 3 + 2];
        //std::cout << "Rank:" << rank << ", " << "x, y, z" << cell_main_x << cell_main_y << cell_main_z << std::endl;

        // define search range
        int i_start = get_startidx_of_surrboxes(i_centered_box);
        int j_start = get_startidx_of_surrboxes(j_centered_box);
        int k_start = get_startidx_of_surrboxes(k_centered_box);
        int i_end = get_endidx_of_surrboxes(cpus_per_side, i_centered_box);
        int j_end = get_endidx_of_surrboxes(cpus_per_side, j_centered_box);
        int k_end = get_endidx_of_surrboxes(cpus_per_side, k_centered_box);

        // initialization
        map_rank_to_n_surrcells[rank] = 0;

        //get number of surrownding cells
        for (int i = i_start; i <= i_end; i++) {
            for (int j = j_start; j <= j_end; j++) {
                for (int k = k_start; k <= k_end; k++) {
                    if ((i != i_centered_box) || (j != j_centered_box) || (k != k_centered_box)) {
                        map_rank_to_n_surrcells[rank]++;
                        //printf("Peripheral boxes of (%d, %d, %d): (x, y, z) = (%d, %d, %d)\n", i_centered_box, j_centered_box, k_centered_box, i, j, k);
                    }
                }
            }
        }
        //printf("N of peripheral boxes for (%d, %d, %d): %d\n", i_centered_box, j_centered_box, k_centered_box, map_rank_to_n_surrcells[rank]);
    }
}
void get_rank_of_surrboxes(const int max_rank, const int cpus_per_side, int *map_rank_to_cell, int *map_cell_to_rank, int **map_rank_to_ranks_of_surrcells){
    for (int rank = 0; rank < max_rank; rank++) {
        int i_centered_box = map_rank_to_cell[rank * 3 + 0];
        int j_centered_box = map_rank_to_cell[rank * 3 + 1];
        int k_centered_box = map_rank_to_cell[rank * 3 + 2];
        //std::cout << "Rank:" << rank << ", " << "x, y, z" << cell_main_x << cell_main_y << cell_main_z << std::endl;

        // define search range
        int i_start = get_startidx_of_surrboxes(i_centered_box);
        int j_start = get_startidx_of_surrboxes(j_centered_box);
        int k_start = get_startidx_of_surrboxes(k_centered_box);
        int i_end = get_endidx_of_surrboxes(cpus_per_side, i_centered_box);
        int j_end = get_endidx_of_surrboxes(cpus_per_side, j_centered_box);
        int k_end = get_endidx_of_surrboxes(cpus_per_side, k_centered_box);

        //loop for surrownding cells
        int count = 0;
        for (int i = i_start; i <= i_end; i++) {
            for (int j = j_start; j <= j_end; j++) {
                for (int k = k_start; k <= k_end; k++) {
                    if ((i != i_centered_box) || (j != j_centered_box) || (k != k_centered_box)) {
                        int surr_rank = map_cell_to_rank[i * cpus_per_side * cpus_per_side + j * cpus_per_side + k];
                        map_rank_to_ranks_of_surrcells[rank][count] = surr_rank;
                        count++;
                        //printf("Surr rank for rank %d: %d\n", rank, surr_rank);
                    }
                }
            }
        }
        //printf("N of peripheral boxes for (%d, %d, %d): %d\n", i_centered_box, j_centered_box, k_centered_box, count);
    }
}
void get_tot_particles_in_surrboxes(const int max_rank, int *map_rank_to_n_particles_in_surrcells, int *map_rank_to_n_surrcells, int **map_rank_to_ranks_of_surrcells, int *n_particles_eachrank) {
    for (int rank = 0; rank < max_rank; rank++){
        map_rank_to_n_particles_in_surrcells[rank] = 0; //init
        for (int i_surr = 0; i_surr < map_rank_to_n_surrcells[rank]; i_surr++) {
            int surr_rank = map_rank_to_ranks_of_surrcells[rank][i_surr];
            map_rank_to_n_particles_in_surrcells[rank] += n_particles_eachrank[surr_rank];
        }
        //printf("Total particles in the surrownding box for rank %d: %d\n", rank, map_rank_to_n_particles_in_surrcells[rank]);
    }
}
void map_rank_to_coords_surrcells(int max_rank, int *map_rank_to_n_particles_in_surrcells, double **map_rank_to_coords_surrbox, double **coords_each_rank, int *n_particles_eachrank, int *map_rank_to_n_surrcells, int **map_rank_to_ranks_of_surrcells){
    // get values
    printf("All loop values \n");
    printf("%d, %d\n", max_rank, map_rank_to_n_surrcells[0], n_particles_eachrank[0]);
    for (int rank; rank < max_rank; rank++){
        int processed_temp_n_particles = 0;
        for (int i_surr_rank = 0; i_surr_rank < map_rank_to_n_surrcells[rank]; i_surr_rank++) {
            int surr_rank = map_rank_to_ranks_of_surrcells[rank][i_surr_rank];
            if (rank == 16) {
                //printf("Original coords for rank %d:\n", surr_rank);
                //print_particles_in_box(coords_each_rank[surr_rank], n_particles_eachrank[surr_rank]);
            }
            for (int j_particle = 0; j_particle < n_particles_eachrank[surr_rank]; j_particle++){
                map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 0] = coords_each_rank[surr_rank][j_particle * 4 + 0];
                map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 1] = coords_each_rank[surr_rank][j_particle * 4 + 1];
                map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 2] = coords_each_rank[surr_rank][j_particle * 4 + 2];
                map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 3] = coords_each_rank[surr_rank][j_particle * 4 + 3];
                if (rank == 16) {
                    printf("%d, (x, y, z, id) = (%.2f, %.2f, %.2f, %d)\n", processed_temp_n_particles,
                           coords_each_rank[surr_rank][j_particle * 4 + 0],
                           coords_each_rank[surr_rank][j_particle * 4 + 1],
                           coords_each_rank[surr_rank][j_particle * 4 + 2],
                           (int) coords_each_rank[surr_rank][j_particle * 4 + 3]);
                }
                /**
                if (rank == 16) {
                    printf("%d, (x, y, z, id) = (%.2f, %.2f, %.2f, %d)\n", processed_temp_n_particles,
                           map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 0],
                           map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 1],
                           map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 2],
                           (int)map_rank_to_coords_surrbox[rank][processed_temp_n_particles * 4 + 3]);
                }
                **/
                processed_temp_n_particles = processed_temp_n_particles + 1;
            }
        }
    }
}
