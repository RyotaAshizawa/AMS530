//
// Created by Ryota Ashizawa on 11/23/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "particle.h"

#ifndef HW4_IN_C_BOX_H
#define HW4_IN_C_BOX_H


/** **?
 *
 * @param box is N * 8 vector. Each column has x, y, z, fx, fy, fz, rank, and id
 * */

// member like functions
void print_particles(double **box, const int N);
void init_coords_and_forces(double **box, bool use_rand, const int N, const int particles_per_side, const double particle_cellsize);
void dump_particles(double **box, char filepath[], const int N);
void init_map_cell_to_rank(int cpus_per_side, int *map_cell_to_rank, int *map_rank_to_cell);
void assign_rank_to_cell(double **box, const int N, int *n_particles_eachrank, int *map_cell_to_rank, const int cpu_per_side, const int cell_len_per_cpu, const int max_rank);
void get_particles_each_rank(double **box, const int N, double **coords_each_rank, const int max_rank);
void print_particles_special_rank(double **coords_each_rank, int *n_particles_eachrank, const int rank);
void print_particles_in_box(double *particles_coords, int particles_in_box);
int get_startidx_of_surrboxes(const int centered_box_idx);
int get_endidx_of_surrboxes(const int cpus_per_side, const int centered_box_idx);
void get_n_surrboxes(const int max_rank, const int cpus_per_side, int *map_rank_to_cell, int *map_rank_to_n_surrcells);
void get_rank_of_surrboxes(const int max_rank, const int cpus_per_side, int *map_rank_to_cell, int *map_cell_to_rank, int **map_rank_to_ranks_of_surrcells);
void get_tot_particles_in_surrboxes(const int max_rank, int *map_rank_to_n_particles_in_surrcells, int *map_rank_to_n_surrcells, int **map_rank_to_ranks_of_surrcells, int *n_particles_eachrank);







#endif //HW4_IN_C_BOX_H
