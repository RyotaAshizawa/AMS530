//
// Created by Ryota Ashizawa on 11/22/22.
//
#include "box.hpp"

Box::Box(const int N, const int box_size, const bool rand_coordinates) {
    // arguments
    this -> N = N;
    this -> box_size = box_size;
    this -> rand_coordinates = rand_coordinates;
    // initialization
    particles_per_side = std::ceil(std::pow(N, 1./3.));
    particle_cellsize = box_size / double(particles_per_side);
    // new
    Particle particles[N];
    // initialization
    init_coords_and_forces(rand_coordinates);
}
// member functions
int Box::get_n_particles(){
    return N;
}
int Box::get_box_size(){
    return box_size;
}
Particle Box::get_particle(const int i){
    return particles[i];
}
void Box::print_particles() {
    for (int i = 0; i < N; i++) {
        particles[i].print();
    }
}
void Box::dump_particles(std::string &filepath){
    std::ofstream ofs(filepath);
    ofs << N << std::endl;
    ofs << "MD example" << std::endl;
    for (int i = 0; i < N; i++) {
        ofs << "He\t" << particles[i].dump_coordinate() << std::endl;
    }
}
void Box::init_coords_and_forces(bool use_rand) {
    // assume each particle locates at the center of each cell
    double x, y, z;
    // rand generator
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_int_distribution<> rand(-20000, 20000);
    // assign values
    int particle_count = 0;
    for (int i = 0; i < particles_per_side; i++) {
        for (int j = 0; j < particles_per_side; j++) {
            for (int k = 0; k < particles_per_side; k++) {
                // set position
                x = particle_cellsize * (i + 0.5);
                y = particle_cellsize * (j + 0.5);
                z = particle_cellsize * (k + 0.5);
                // add random perturbation
                if (use_rand){
                    x = x + particle_cellsize * (rand(mt) / 100000.0);
                    y = y + particle_cellsize * (rand(mt) / 100000.0);
                    z = z + particle_cellsize * (rand(mt) / 100000.0);
                }
                particles[particle_count].set_coordinate(x, y, z, 0, 0, 0, particle_count);
                // Break the loop if the positions of all particles are set.
                particle_count++;
                if (particle_count == N) {break;}
            }
        }
    }
}
void Box::set_coords(double *coords){
    for (int i = 0; i < N; i++){
        coords[i * 4 + 0] = particles[i].get_x();
        coords[i * 4 + 1] = particles[i].get_y();
        coords[i * 4 + 2] = particles[i].get_z();
        coords[i * 4 + 3] = particles[i].get_id();
    }
}

