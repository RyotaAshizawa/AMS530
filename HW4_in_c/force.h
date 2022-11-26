//
// Created by Ryota Ashizawa on 11/24/22.
//

#include <math.h>

#ifndef HW4_IN_C_FORCE_H
#define HW4_IN_C_FORCE_H


double cal_dist_two_particles(double *j_coords, double *i_coords);
void add_force_between_two_particles_to_vector(double *forces_and_id, double *j_coords, double *i_coords, const int rank, double cutoff);

#endif //HW4_IN_C_FORCE_H
