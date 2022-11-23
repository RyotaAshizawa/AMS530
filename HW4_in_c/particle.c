//
// Created by Ryota Ashizawa on 11/23/22.
//

#include "stdio.h"
#include "particle.h"

/** **?
 *
 * @param box is N * 7 vector. Each column has x, y, z, fx, fy, fz, rank, and id
 * */

void set_coordinate(double *particle, const double x, const double y, const double z, const double fx, const double fy, const double fz, const double rank, const double id){
    particle[0] = x;
    particle[1] = y;
    particle[2] = z;
    particle[3] = fx;
    particle[4] = fy;
    particle[5] = fz;
    particle[6] = rank;
    particle[7] = id;
}
void set_rank(double *particle, const int rank){
    particle[6] = rank;
}

// get
int get_rank(double *particle){
    return particle[6];
}
double get_x(double *particle){
    return particle[0];
}
double get_y(double *particle){
    return particle[1];
}
double get_z(double *particle){
    return particle[2];
}
int get_id(double *particle){
    return particle[7];
}

// file printing
void print_particle(double *particle){
    printf("(x, y, z, fx, fy, fz, rank, id) = (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %d, %d)\n",
           particle[0], particle[1], particle[2], particle[3], particle[4], particle[5], (int)particle[6], (int)particle[7]);
}
/**
std::string dump_coordinate(){
    return std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(z);
}
**/
