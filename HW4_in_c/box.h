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
void dump_particles(double **box, char filepath[]);
void init_map_cell_to_rank(double **box, int cpus_per_side, int ***map_cell_to_rank, int *map_rank_to_cell);
void assign_rank_to_box(double **box, int *n_particles_eachrank, int ***map_cell_to_rank, const int cell_len_per_cpu, const int max_rank);



#endif //HW4_IN_C_BOX_H
