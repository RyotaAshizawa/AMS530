//
// Created by Ryota Ashizawa on 11/24/22.
//

#include "force.h"
#include <stdio.h>

double cal_dist_two_particles(double *j_coords, double *i_coords){
    double x_j = j_coords[0];
    double y_j = j_coords[1];
    double z_j = j_coords[2];
    double x_i = i_coords[0];
    double y_i = i_coords[1];
    double z_i = i_coords[2];

    double dist = sqrt((pow(x_j - x_i, 2) + pow(y_j - y_i, 2) + pow(z_j - z_i, 2)));
    return dist;
}
void add_force_between_two_particles_to_vector(double *forces_and_id, double *j_coords, double *i_coords, int rank) {
    double dist = cal_dist_two_particles(i_coords, j_coords);
    double force_term1_x = (j_coords[0] - i_coords[0]) / pow(dist, 8);
    double force_term1_y = (j_coords[1] - i_coords[1]) / pow(dist, 8);
    double force_term1_z = (j_coords[2] - i_coords[2]) / pow(dist, 8);
    double force_term2_x = (j_coords[0] - i_coords[0]) / pow(dist, 14);
    double force_term2_y = (j_coords[1] - i_coords[1]) / pow(dist, 14);
    double force_term2_z = (j_coords[2] - i_coords[2]) / pow(dist, 14);

    double force_x = force_term1_x + force_term2_x;
    double force_y = force_term1_y + force_term2_y;
    double force_z = force_term1_z + force_term2_z;
    if (isnan(force_x) || isnan(force_y) || isnan(force_z)){
        printf("Nan happend for ((x1, y1, z1, id), (x2, y2, z2, id), dist) = ((%f, %f, %f, %d), (%f, %f, %f, %d), %f) in the rank %d\n",
               i_coords[0], i_coords[1], i_coords[2], (int)i_coords[3],
               j_coords[0], j_coords[1], j_coords[2], (int)j_coords[3],
               dist, rank);
    }

    forces_and_id[0] += force_x;
    forces_and_id[1] += force_y;
    forces_and_id[2] += force_z;
    forces_and_id[3] = i_coords[3];
}
