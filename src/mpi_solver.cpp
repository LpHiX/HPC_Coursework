/**
 * @file mpi_solver.cpp
 * @author Martin Leung
 */

#include "mpi_solver.h"
#include <iostream>
#include <cmath>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <fstream>

MPISolver::MPISolver(int Nx, int Ny, int Nz, int Px, int Py, int Pz, double epsilon):
Nx(Nx),
Ny(Ny),
Nz(Nz),
nx(Nx - 2),
ny(Ny - 2),
nz(Nz - 2),
Px(Px),
Py(Py),
Pz(Pz),
epsilon(epsilon)
{
    const int dims = 3;
    int sizes[dims] = {Px, Py, Pz};
    int periods[dims] = {0, 0, 0};
    int reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &gridcomm);

    // gridrank;
    MPI_Comm_rank(gridcomm, &gridrank);

    MPI_Cart_shift(gridcomm, 0, 1, &neighbours_ranks[NEG_X], &neighbours_ranks[POS_X]);
    MPI_Cart_shift(gridcomm, 1, 1, &neighbours_ranks[NEG_Y], &neighbours_ranks[POS_Y]);
    MPI_Cart_shift(gridcomm, 2, 1, &neighbours_ranks[NEG_Z], &neighbours_ranks[POS_Z]);

    // Handle remainders

    int quotient_x = nx / Px;
    int quotient_y = ny / Py;
    int quotient_z = nz / Pz;
    int remainder_x = nx % Px;
    int remainder_y = ny % Py;
    int remainder_z = nz % Pz;

    int coords[3];
    MPI_Cart_coords(gridcomm, gridrank, 3, coords);
    int Pi = coords[0];
    int Pj = coords[1];
    int Pk = coords[2];

    lnx = (Pi < remainder_x) ? quotient_x + 1 : quotient_x;
    lny = (Pj < remainder_y) ? quotient_y + 1 : quotient_y;
    lnz = (Pk < remainder_z) ? quotient_z + 1 : quotient_z;

    lNx = lnx + 2;
    lNy = lny + 2;
    lNz = lnz + 2;



    // N distribution: (13)
    // 0 1 2 3 4 5 6 7 8 9 10 11 12
    // 0 1 2 3 4 5 
    //         4 5 6 7 8 9
    //                 8 9 10 11 12
    // Pi * (quotient + 1) (quotient = 3)
    // 0 1 2 3
    //         4 5 6 7
    //                 8 9 10
    start_i = (Pi < remainder_x) ? Pi * (quotient_x + 1) : Pi * quotient_x + remainder_x;
    start_j = (Pj < remainder_y) ? Pj * (quotient_y + 1) : Pj * quotient_y + remainder_y;
    start_k = (Pk < remainder_z) ? Pk * (quotient_z + 1) : Pk * quotient_z + remainder_z;

    localstate = std::make_unique<JacobiLocalState>(Nx, Ny, Nz, lNx, lNy, lNz);

    for (int i = 0; i < 6; i++){
        if(neighbours_ranks[i] != MPI_PROC_NULL){
            int plane_size = localstate->get_planesize(static_cast<Direction>(i));
            send_buffers[i] = new double[plane_size]();
            recv_buffers[i] = new double[plane_size]();
        } 
    }
    
}

void MPISolver::initialize(int test){
    localstate->apply_test_conditions(test, start_i, start_j, start_k);
}

void MPISolver::initialize(std::string filename){
    localstate->read_set_forcing(filename, start_i, start_j, start_k);
}

void MPISolver::solve(){
    boost::timer::cpu_timer wall_timer;
    double last_report_time = -2;

    int iterations = 0;
    double residual = get_residual();

    while (residual > epsilon){
        localstate->jacobi_step();
        
        MPI_Request requests[12];
        int req_count = 0;

        for (int i = 0; i < 6; i++){
            if (neighbours_ranks[i] != MPI_PROC_NULL){
                int buf_size = localstate->get_planesize(static_cast<Direction>(i));
                MPI_Irecv(recv_buffers[i], buf_size, MPI_DOUBLE, neighbours_ranks[i], 0, gridcomm, &requests[req_count]);
                req_count++;
            }
        }
        for (int i = 0; i < 6; i++){
            if (neighbours_ranks[i] != MPI_PROC_NULL){
                Direction dir = static_cast<Direction>(i);
                int buf_size = localstate->get_planesize(dir);
                localstate->get_u_boundary(send_buffers[i], dir);
                MPI_Isend(send_buffers[i], buf_size, MPI_DOUBLE, neighbours_ranks[i], 0, gridcomm, &requests[req_count]);
                req_count++;
            }
        }
        
        MPI_Waitall(req_count, requests, MPI_STATUS_IGNORE);
        for (int i = 0; i < 6; i++){
            if (neighbours_ranks[i] != MPI_PROC_NULL){
                Direction dir = static_cast<Direction>(i);
                localstate->set_u_boundary(recv_buffers[i], dir);
            }
        }
        if (iterations % 100 == 0) {
            residual = get_residual();
            double current_time = wall_timer.elapsed().wall / 1e9; // Convert nanoseconds to seconds        
            if (current_time - last_report_time >= 2.0) {
                if (gridrank == 0) {
                    std::cout 
                    << std::setw(15) << "Step: " << std::setw(10) << iterations 
                    << std::setw(15) << "Time: " << std::setw(10) << current_time << "s"
                    << std::setw(15) << "Residual: " << std::setw(15) << residual << std::endl;
                    last_report_time = current_time;
                }
                last_report_time = current_time;
            }
        }
        iterations++;
    }
    if (gridrank == 0) {
        
        std::cout << "Converged in " << iterations << " iterations, residual: " << residual << std::endl;
    }
}

MPISolver::~MPISolver(){
    for (int i = 0; i < 6; i++){
        if(neighbours_ranks[i] != MPI_PROC_NULL){
            delete[] send_buffers[i];
            delete[] recv_buffers[i];
        }
    }
}

double MPISolver::get_residual() const{
    double local_residual = localstate->get_residualsquared();
    double residual_square_sum = 0;
    MPI_Allreduce(&local_residual, &residual_square_sum, 1, MPI_DOUBLE, MPI_SUM, gridcomm);
    return sqrt(residual_square_sum);
}

void MPISolver::write_solution(std::string filename, int test){ // This destroys u2 to save memory
    localstate->pack_solution_u2();

    double *recv_interior_u = nullptr;
    int* displs = nullptr;
    int* recvcounts = nullptr;
    int total_ranks = Px*Py*Pz;

    if (gridrank == 0){
        recv_interior_u = new double[nx * ny * nz];
        displs = new int[total_ranks]();
        recvcounts = new int[total_ranks];
    }

    int sendcount = lnx * lny * lnz;
    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, gridcomm);
    if (gridrank == 0){
        int displcount = 0;
        for (int i = 0; i < total_ranks; i++){
            displs[i] = displcount;
            displcount += recvcounts[i];
        }
    }

    MPI_Gatherv(localstate->lu2, sendcount, MPI_DOUBLE, recv_interior_u, recvcounts, displs, MPI_DOUBLE, 0, gridcomm);

    if (gridrank == 0){
        std::ofstream fileOutput(filename, std::ios::out | std::ios::trunc);
        fileOutput.precision(std::numeric_limits<double>::max_digits10);
        fileOutput << Nx << " " << Ny << " " << Nz << std::endl;
        if (!fileOutput.good()){
            throw std::runtime_error("Could not open output file");
        }
        int quotient_x = nx / Px;
        int quotient_y = ny / Py;
        int quotient_z = nz / Pz;
        int remainder_x = nx % Px;
        int remainder_y = ny % Py;
        int remainder_z = nz % Pz;
        double hx = 1.0 / (Nx - 1);
        double hy = 1.0 / (Ny - 1);
        double hz = 1.0 / (Nz - 1);

        int idx = 0;
        for (int p = 0; p < total_ranks; p++){
            int coords[3];
            MPI_Cart_coords(gridcomm, p, 3, coords);
            int Pi_p = coords[0];
            int Pj_p = coords[1];
            int Pk_p = coords[2];

            int lnx_p = (Pi_p < remainder_x) ? quotient_x + 1 : quotient_x;
            int lny_p = (Pj_p < remainder_y) ? quotient_y + 1 : quotient_y;
            int lnz_p = (Pk_p < remainder_z) ? quotient_z + 1 : quotient_z;

            int start_i_p = (Pi_p < remainder_x) ? Pi_p * (quotient_x + 1) : Pi_p * quotient_x + remainder_x;
            int start_j_p = (Pj_p < remainder_y) ? Pj_p * (quotient_y + 1) : Pj_p * quotient_y + remainder_y;
            int start_k_p = (Pk_p < remainder_z) ? Pk_p * (quotient_z + 1) : Pk_p * quotient_z + remainder_z;
            
            for (int k = 0; k < lnz_p; k++){
                double z = (start_k_p + k + 1) * hz;
                for (int j = 0; j < lny_p; j++){
                    double y = (start_j_p + j + 1) * hy;
                    for (int i = 0; i < lnx_p; i++){
                        double x = (start_i_p + i + 1) * hx;
                        fileOutput << x << " " << y << " " << z << " " << recv_interior_u[idx] << "\n"; 
                        idx ++;
                    }
                }
            }
        }
        for (int k = 0; k < Nz; k++){
            double z = k * hz;
            for (int j = 0; j < Ny; j++){
                double y = j * hy;
                for (int i = 0; i < Nx; i++){
                    if (i == 0 || i == Nx-1 || 
                        j == 0 || j == Ny-1 || 
                        k == 0 || k == Nz-1) {
                        double x = i * hx;
                        double val = (test == 1) ? x*x + y*y + z*z : 0.0;
                        fileOutput << x << " " << y << " " << z << " " << val << '\n'; 
                    }
                }
            }
        }

        fileOutput.close();
        delete[] recvcounts;
        delete[] displs;
        delete[] recv_interior_u;
    }
}