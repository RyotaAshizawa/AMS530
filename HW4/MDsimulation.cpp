//
// Created by Ryota Ashizawa on 11/20/22.
//

#include <iostream>
#include <mpi.h>
#include <cmath>
#include <random>
#include <fstream>

class Particle{
public:
    Particle(){
        x = 0;
        y = 0;
        z = 0;
        fx = 0;
        fy = 0;
        fz = 0;
        rank = 0;
        id = 0;
    }
    // set
    void set_coordinate(const double x, const double y, const double z, const double fx, const double fy, const double fz, const int id){
        this -> x = x;
        this -> y = y;
        this -> z = z;
        this -> fx = fx;
        this -> fy = fy;
        this -> fz = fz;
        this -> id = id;
    }
    void set_rank(const int rank){
        this -> rank = rank;
    }
    // get
    int get_rank(){
        return rank;
    }
    double get_x(){
        return x;
    }
    double get_y(){
        return y;
    }
    double get_z(){
        return z;
    }
    int get_id(){
        return id;
    }
    // file printing
    void print(){
        std::cout << "(x, y, z, fx, fy, fz, rank, id) = (" << x << ", " << y << ", " << z << ", " << fx << ", " << fy << ", " << fz << ", " << rank << ", " << id << ")" << std::endl;
    }
    std::string dump_coordinate(){
        return std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(z);
    }
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

class Particles{
public:
    Particles(const int N, const int box_size, const int size, const bool rand_coordinates){
        // arguments
        this -> N = N;
        this -> box_size = box_size;
        this -> size = size;
        this -> cpus_per_side = std::floor(std::pow(size, 1./3.));
        this -> pcells_per_side = std::ceil(std::pow(N, 1./3.));
        this -> pcell_size = box_size / double(pcells_per_side);
        this -> ccell_size = box_size / double(cpus_per_side);
        this -> max_rank = cpus_per_side * cpus_per_side * cpus_per_side;
        // new
        particles = new Particle[N];
        n_eachrank = new int[max_rank];
        map_rank_to_cell = new int[3 * max_rank];
        map_rank_to_n_surrownd_cells = new int[max_rank];
        for (int i = 0; i < max_rank; i++) {
            map_rank_to_n_surrownd_cells[i] = 0;
        }
        map_rank_to_surrownd_ranks = new int*[max_rank];
        particles_coord_rank = new double*[max_rank];
        map_cell_to_rank = new int**[cpus_per_side];
        for (int i = 0; i < cpus_per_side; i++) {
            map_cell_to_rank[i] = new int*[cpus_per_side];
            for (int j = 0; j < cpus_per_side; j++) {
                map_cell_to_rank[i][j] = new int[cpus_per_side];
            }
        }
        // initialization
        init_coords_and_forces(rand_coordinates);
        init_map_cell_to_rank();
        init_rank();
        init_mpi_particles_in_mainbox();
        init_map_rank_to_n_peripheral_boxes();
    }
    ~Particles(){
        delete[] particles;
        delete[] n_eachrank;
        delete[] map_rank_to_cell;
        delete[] map_rank_to_n_surrownd_cells;
        //delete[] particles_each_rank;
    }
    // member functions
    void print(){
        for (int i = 0; i < N; i++){
            particles[i].print();
        }
    }
    void print_rank(const int rank){
        std::cout << "All particles on rank:" << rank << std::endl;
        for (int i = 0; i < n_eachrank[rank]; i++){
            double x = particles_coord_rank[rank][i * 4 + 0];
            double y = particles_coord_rank[rank][i * 4 + 1];
            double z = particles_coord_rank[rank][i * 4 + 2];
            std::cout << "(x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
        }
    }
    // get
    int get_n_of_rank(const int rank){
        return n_eachrank[rank];
    }
    int get_start_idx_peripheral_box(const int main_box_idx){
        int start_idx = main_box_idx - 1;
        if (start_idx < 0) {
            start_idx = 0;
        }
        return start_idx;
    }
    int get_end_idx_peripheral_box(const int main_box_idx){
        int end_idx = main_box_idx + 1;
        if (cpus_per_side <= end_idx) {
            end_idx = cpus_per_side - 1;
        }
        return end_idx;
    }
    void mpi_send_n_eachrank(MPI_Comm comm, MPI_Request *request, const int tag){
        for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
            MPI_Isend(n_eachrank, max_rank, MPI_INT, dst_rank, tag, comm, request);
        }
    }
    void mpi_send_n_surrboxes(MPI_Comm comm, MPI_Request *request, const int tag){
        for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
            MPI_Isend(map_rank_to_n_surrownd_cells, max_rank, MPI_INT, dst_rank, tag, comm, request);
            //std::cout << map_rank_to_n_surrownd_cells[dst_rank] << std::endl;
        }
    }
    void mpi_send_map_rank_to_surrranks(MPI_Comm comm, MPI_Request *request, const int tag) {
        for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
            MPI_Isend(map_rank_to_surrownd_ranks[dst_rank], map_rank_to_n_surrownd_cells[dst_rank], MPI_INT, dst_rank, tag, comm, request);
        }
    }
    void mpi_send_mainbox_particles(MPI_Comm comm, MPI_Request *request, const int tag){
        for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
            MPI_Isend(particles_coord_rank[dst_rank], n_eachrank[dst_rank] * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
        }
    }
    void mpi_send_surrboxes_particles(MPI_Comm comm, MPI_Request *request, const int tag){
        for (int dst_rank = 0; dst_rank < max_rank; dst_rank++){
            // get particles sum in the surr box
            int n_in_surrbox = 0;
            for (int i = 0; i < map_rank_to_n_surrownd_cells[dst_rank]; i++){ // loop for surr box times
                int surr_rank = map_rank_to_surrownd_ranks[dst_rank][i];
                n_in_surrbox += n_eachrank[surr_rank];
            }
            //std::cout << "dst_rank:" << dst_rank << "," << n_in_surrbox << "GAAAAAAAAA\n";

            // explictly create passing array
            int count = 0;
            double *surrbox_coords = new double[n_in_surrbox * 4];
            for (int surr_i = 0; surr_i < map_rank_to_n_surrownd_cells[dst_rank]; surr_i++) {
                int surr_rank = map_rank_to_surrownd_ranks[dst_rank][surr_i];
                //std::cout << "(surr_rank, n)" << surr_rank << "," << n_eachrank[surr_rank] << std::endl;
                for (int atom_n = 0; atom_n < n_eachrank[surr_rank]; atom_n++){
                    surrbox_coords[(count + atom_n) * 4 + 0] = particles_coord_rank[surr_rank][atom_n * 4 + 0];
                    surrbox_coords[(count + atom_n) * 4 + 1] = particles_coord_rank[surr_rank][atom_n * 4 + 1];
                    surrbox_coords[(count + atom_n) * 4 + 2] = particles_coord_rank[surr_rank][atom_n * 4 + 2];
                    surrbox_coords[(count + atom_n) * 4 + 3] = particles_coord_rank[surr_rank][atom_n * 4 + 3];
                    //std::cout << "Rank:" << dst_rank <<  "(x, y, z) = (" << particles_coord_rank[surr_rank][atom_n * 4 + 0] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 1] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 2] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 3] << ")" << std::endl;
                    //std::cout << "(x, y, z) = (" << surrbox_coords[surr_i * n_eachrank[surr_rank] + atom_n * 4 + 0] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 1] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 2] << ", " << particles_coord_rank[surr_rank][atom_n * 4 + 3] << ")" << std::endl;
                    count++;
                }
            }
            // mpi send
            MPI_Isend(surrbox_coords, n_in_surrbox * 4, MPI_DOUBLE, dst_rank, tag, comm, request);
        }
    }
    // file
    void dump_particles(std::string &filepath){
        std::ofstream ofs(filepath);
        ofs << N << std::endl;
        ofs << "MD example" << std::endl;
        for (int i = 0; i < N; i++) {
            ofs << "He\t" << particles[i].dump_coordinate() << std::endl;
        }
    }

private:
    // variables
    int N;
    int size;
    int max_rank;
    int cpus_per_side;
    int pcells_per_side;
    int *n_eachrank;
    int *map_rank_to_n_surrownd_cells;
    int *map_rank_to_cell;
    int **map_rank_to_surrownd_ranks;
    int ***map_cell_to_rank;
    double box_size;
    double pcell_size; // particle cell size
    double ccell_size; // particle cell size
    double **particles_coord_rank; // t
    Particle *particles;


    // functions
    void init_coords_and_forces(bool use_rand) {
        // assume each particle locates at the center of each cell
        double x, y, z;
        // rand generator
        std::random_device rnd;
        std::mt19937 mt(rnd());
        std::uniform_int_distribution<> rand(-20000, 20000);
        // assign values
        int particle_count = 0;
        for (int i = 0; i < pcells_per_side; i++) {
            for (int j = 0; j < pcells_per_side; j++) {
                for (int k = 0; k < pcells_per_side; k++) {
                    // set position
                    x = pcell_size * (i + 0.5);
                    y = pcell_size * (j + 0.5);
                    z = pcell_size * (k + 0.5);
                    // add random perturbation
                    if (use_rand){
                        x = x + pcell_size * (rand(mt) / 100000.0);
                        y = y + pcell_size * (rand(mt) / 100000.0);
                        z = z + pcell_size * (rand(mt) / 100000.0);
                    }
                    particles[particle_count].set_coordinate(x, y, z, 0, 0, 0, particle_count);
                    // Break the loop if the positions of all particles are set.
                    particle_count++;
                    if (particle_count == N) {break;}
                }
            }
        }
    }
    void init_map_cell_to_rank() {
        int rank = 0;
        for (int i = 0; i < cpus_per_side; i++) {
            for (int j = 0; j < cpus_per_side; j++) {
                for (int k = 0; k < cpus_per_side; k++) {
                    map_cell_to_rank[i][j][k] = rank;
                    map_rank_to_cell[rank * 3 + 0] = i;
                    map_rank_to_cell[rank * 3 + 1] = j;
                    map_rank_to_cell[rank * 3 + 2] = k;
                    rank++;
                }
            }
        }
    }
    void init_rank() {
        // array definition and initialize
        for (int i = 0; i < max_rank; i++){
            n_eachrank[i] = 0;
        }
        for (int i = 0; i < N; i++){
            int x_cell = floor(particles[i].get_x() / ccell_size);
            int y_cell = floor(particles[i].get_y() / ccell_size);
            int z_cell = floor(particles[i].get_z() / ccell_size);
            int rank = map_cell_to_rank[x_cell][y_cell][z_cell];
            particles[i].set_rank(rank);
            n_eachrank[rank]++;
        }
    }
    // mpi
    void init_mpi_particles_in_mainbox(){
        // Allocate memory for each rank
        for (int i = 0; i < max_rank; i++) {
            particles_coord_rank[i] = new double[n_eachrank[i] * 4]; // *3 to store x, y, z, id
        }
        // array definition and initialize
        int n_assigned_p_each_rank[max_rank];
        for (int i = 0; i < max_rank; i++){
            n_assigned_p_each_rank[i] = 0;
        }
        for (int i = 0; i < N; i++){
            int rank = particles[i].get_rank();
            particles_coord_rank[rank][n_assigned_p_each_rank[rank] * 4 + 0] = particles[i].get_x();
            particles_coord_rank[rank][n_assigned_p_each_rank[rank] * 4 + 1] = particles[i].get_y();
            particles_coord_rank[rank][n_assigned_p_each_rank[rank] * 4 + 2] = particles[i].get_z();
            particles_coord_rank[rank][n_assigned_p_each_rank[rank] * 4 + 3] = particles[i].get_id();
            n_assigned_p_each_rank[rank]++;
        }
    }
    // allocate
    void init_map_rank_to_n_peripheral_boxes(){
        for (int rank = 0; rank < max_rank; rank++) {
            int cell_main_x = map_rank_to_cell[rank * 3 + 0];
            int cell_main_y = map_rank_to_cell[rank * 3 + 1];
            int cell_main_z = map_rank_to_cell[rank * 3 + 2];
            //std::cout << "Rank:" << rank << ", " << "x, y, z" << cell_main_x << cell_main_y << cell_main_z << std::endl;

            // define search range
            int start_x = get_start_idx_peripheral_box(cell_main_x);
            int start_y = get_start_idx_peripheral_box(cell_main_y);
            int start_z = get_start_idx_peripheral_box(cell_main_z);
            int end_x = get_end_idx_peripheral_box(cell_main_x);
            int end_y = get_end_idx_peripheral_box(cell_main_y);
            int end_z = get_end_idx_peripheral_box(cell_main_z);

            //get number of surrownding cells
            for (int i = start_x; i <= end_x; i++) {
                for (int j = start_y; j <= end_y; j++) {
                    for (int k = start_z; k <= end_z; k++) {
                        if ((i != cell_main_x) || (j != cell_main_y) || (k != cell_main_z)) {
                            map_rank_to_n_surrownd_cells[rank]++;
                            //std::cout << "Peripheral boxes for rank:" << rank << ", " << "x, y, z" << i << j << k << "rank:" << map_cell_to_rank[i][j][k]  << std::endl;
                        }
                    }
                }
            }

            // assign peripheral ranks
            map_rank_to_surrownd_ranks[rank] = new int[map_rank_to_n_surrownd_cells[rank]];
            int count = 0;
            for (int i = start_x; i <= end_x; i++) {
                for (int j = start_y; j <= end_y; j++) {
                    for (int k = start_z; k <= end_z; k++) {
                        if ((i != cell_main_x) || (j != cell_main_y) || (k != cell_main_z)) {
                            map_rank_to_surrownd_ranks[rank][count] = map_cell_to_rank[i][j][k];
                            count ++;
                        }
                    }
                }
            }
            /**
            std::cout << "Rank:" << rank << std::endl;
            for (int i = 0; i < map_rank_to_n_surrownd_cells[rank]; i++) {
                std::cout << map_rank_to_surrownd_ranks[rank][i] << std::endl;
            }
            **/
        }
    }
};

int main(int argc, char *argv[]){
    // initialize mpi
    int rank, size, i;
    int max_rank = std::pow(std::floor(std::pow(size, 1./3.)), 3);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request request;

    // set variables
    const int N = std::stoi(argv[1]);
    const double box_size = std::stof(argv[2]);
    const int mpi_tag_n_per_proc = 0;
    const int mpi_tag_main_box = 1;
    const int mpi_tag_nsurrbox_per_proc = 2;
    const int mpi_tag_rank_to_surrranks = 3;
    const int mpi_tag_surr_box = 4;

    // set particles
    if (rank == 0) {
        // make particles
        Particles particles(N, box_size, size, true);
        particles.print();
        // dump particles
        std::string file_path = "./test.xyz";
        particles.dump_particles(file_path);
        // Send info
        particles.mpi_send_n_eachrank(MPI_COMM_WORLD, &request, mpi_tag_n_per_proc);
        //particles.mpi_send_mainbox_particles(MPI_COMM_WORLD, &request, mpi_tag_main_box);
        //particles.mpi_send_n_surrboxes(MPI_COMM_WORLD, &request, mpi_tag_nsurrbox_per_proc);
        //particles.mpi_send_map_rank_to_surrranks(MPI_COMM_WORLD, &request, mpi_tag_rank_to_surrranks);
        //particles.mpi_send_surrboxes_particles(MPI_COMM_WORLD, &request, mpi_tag_surr_box);
        particles.print_rank(3);
    }

    // recv number of particles in main box
    int *n_eachrank = new int[max_rank];
    MPI_Irecv(n_eachrank, max_rank, MPI_INT, 0, mpi_tag_n_per_proc, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);


    /**
    // recv particles in the main box for this processor
    double *particles_main_box = new double [n_eachrank[rank] * 4];
    MPI_Irecv(particles_main_box, n_eachrank[rank] * 4, MPI_DOUBLE, 0, mpi_tag_main_box, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
    for (int i = 0; i < n_eachrank[rank]; i++){
        if (rank == 2) {
            std::cout << "(rank, x, y, z, i) = (" << rank << ", " << particles_main_box[i * 4] << ", " << particles_main_box[i * 4 + 1] << ", " << particles_main_box[i * 4 + 2] << ", " << particles_main_box[i * 4 + 3]  << ")"  << std::endl;
        }
    }

    // recv n of surrownding boxes
    int *n_surrbox_eachrank = new int[max_rank];
    MPI_Irecv(n_surrbox_eachrank, max_rank, MPI_INT, 0, mpi_tag_nsurrbox_per_proc, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    //recv map i to surr ranks
    int *map_i_to_surrranks = new int[n_surrbox_eachrank[rank]];
    MPI_Irecv(map_i_to_surrranks, n_surrbox_eachrank[rank], MPI_INT, 0, mpi_tag_rank_to_surrranks, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    // get particles in the surr box
    int n_particles_in_surrbox = 0;
    for (int i = 0; i < n_surrbox_eachrank[rank]; i++){ // loop for surr box times
        n_particles_in_surrbox += n_eachrank[map_i_to_surrranks[i]];
    }

    //distribute peripheral particles on each proc
    double *surr_box_coords = new double[n_particles_in_surrbox * 4];
    MPI_Irecv(surr_box_coords, n_particles_in_surrbox * 4, MPI_DOUBLE, 0, mpi_tag_surr_box, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    if (rank == 3) {
        for (int i = 0; i < n_particles_in_surrbox; i++) {
            std::cout << "Passed data: (rank, x, y, z, i) = (" << rank << ", " << surr_box_coords[i * 4] << ", "
                      << surr_box_coords[i * 4 + 1] << ", " << surr_box_coords[i * 4 + 2] << ", "
                      << surr_box_coords[i * 4 + 3] << ")" << std::endl;
        }
    }
    **/

    // Calc force
    // return




    /**
    // delete heap
    delete[] particles_main_box;
    delete[] n_mainbox_eachrank;
    **/

    MPI_Finalize();
    return 0;
}