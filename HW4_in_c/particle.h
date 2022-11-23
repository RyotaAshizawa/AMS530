//
// Created by Ryota Ashizawa on 11/23/22.
//

#ifndef HW4_IN_C_PARTICLE_H
#define HW4_IN_C_PARTICLE_H

extern void set_coordinate(double *particle, const double x, const double y, const double z, const double fx, const double fy, const double fz, const double rank, const double id);
extern void set_rank(double *particle, const int rank);

// get
int get_rank(double *particle);
double get_x(double *particle);
double get_y(double *particle);
double get_z(double *particle);
int get_id(double *particle);
void print_particle(double *particle);


#endif //HW4_IN_C_PARTICLE_H
