//
// Created by Ryota Ashizawa on 11/23/22.
//

#include "stdio.h"
#include "particle.h"

void set_coordinate(Particle *particle, const double x, const double y, const double z, const double fx, const double fy, const double fz, const int rank, const int id){
    particle -> x = x;
    particle -> y = y;
    particle -> z = z;
    particle -> fx = fx;
    particle -> fy = fy;
    particle -> fz = fz;
    particle -> id = id;
}
void set_rank(struct Particle *particle, const int rank){
    particle -> rank = rank;
}

// get
int get_rank(struct Particle *particle){
    return particle -> rank;
}
double get_x(struct Particle *particle){
    return particle -> x;
}
double get_y(struct Particle *particle){
    return particle -> y;
}
double get_z(struct Particle *particle){
    return particle -> z;
}
int get_id(struct Particle *particle){
    return particle -> id;
}

// file printing
void print_particle(struct Particle *particle){
    printf("(x, y, z, fx, fy, fz, rank, id) = (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %d, %d)\n",
           particle->x, particle->y, particle->z, particle->fx, particle->fy, particle->fz, particle->rank, particle->id);
}
/**
std::string dump_coordinate(){
    return std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(z);
}
**/
