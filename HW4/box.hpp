//
// Created by Ryota Ashizawa on 11/22/22.
//
#include <fstream>
#include <random>
#include "particle.hpp"

#ifndef HW4_BOX_HPP
#define HW4_BOX_HPP

class Box {
public:
    Box(const int N = 0, const int box_size = 0, const bool rand_coordinates = false);
    ~Box();

    // member functions
    void print_particles();
    void init_coords_and_forces(bool use_rand);
    void dump_particles(std::string &filepath);
    int get_n_particles();

private:
    int N;
    int box_size;
    int particles_per_side;
    double particle_cellsize;
    bool rand_coordinates;
    Particle *particles;
};

#endif //HW4_BOX_HPP
