//
// Created by Ryota Ashizawa on 11/23/22.
//

#ifndef HW4_IN_C_PARTICLE_H
#define HW4_IN_C_PARTICLE_H

struct Particle{
    double x;
    double y;
    double z;
    double fx;
    double fy;
    double fz;
    int rank;
    int id;
};
typedef struct Particle Particle;

extern void set_coordinate(Particle *particle, const double x, const double y, const double z, const double fx, const double fy, const double fz, const int rank, const int id);
extern void set_rank(struct Particle *particle, const int rank);

// get
int get_rank(struct Particle *particle);
double get_x(struct Particle *particle);
double get_y(struct Particle *particle);
double get_z(struct Particle *particle);
int get_id(struct Particle *particle);
void print_particle(struct Particle *particle);


#endif //HW4_IN_C_PARTICLE_H
