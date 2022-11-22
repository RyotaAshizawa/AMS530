//
// Created by Ryota Ashizawa on 11/22/22.
//
#include <iostream>

#ifndef HW4_PARTICLE_HPP
#define HW4_PARTICLE_HPP

class Particle{
public:
    Particle();
    // set
    void set_coordinate(const double x, const double y, const double z, const double fx, const double fy, const double fz, const int id);
    void set_rank(const int rank);
    // get
    int get_rank(const int rank);
    double get_x();
    double get_y();
    double get_z();
    int get_id();
    // stream
    void print();
    // io
    std::string dump_coordinate();

private:
    double x;
    double y;
    double z;
    double fx;
    double fy;
    double fz;
    int rank;
    int id;
};

#endif //HW4_PARTICLE_HPP
