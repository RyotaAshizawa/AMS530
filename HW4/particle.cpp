//
// Created by Ryota Ashizawa on 11/22/22.
//
#include "particle.hpp"

Particle::Particle(){
    x = 0;
    y = 0;
    z = 0;
    fx = 0;
    fy = 0;
    fz = 0;
    rank = 0;
    id = 0;
}

void Particle::set_coordinate(const double x, const double y, const double z, const double fx, const double fy, const double fz, const int id){
    this -> x = x;
    this -> y = y;
    this -> z = z;
    this -> fx = fx;
    this -> fy = fy;
    this -> fz = fz;
    this -> id = id;
}
void Particle::set_rank(const int rank){
    this -> rank = rank;
}

// get
int Particle::get_rank(const int rank){
    return rank;
}
double Particle::get_x(){
    return x;
}
double Particle::get_y(){
    return y;
}
double Particle::get_z(){
    return z;
}
int Particle::get_id(){
    return id;
}

// file printing
void Particle::print(){
    std::cout << "(x, y, z, fx, fy, fz, rank, id) = (" << x << ", " << y << ", " << z << ", " << fx << ", " << fy << ", " << fz << ", " << rank << ", " << id << ")" << std::endl;
}
std::string Particle::dump_coordinate(){
    return std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(z);
}
