/**
 * @file mpi_solver.cpp
 * @author Martin Leung
 */

#include "mpi_solver.h"
#include <iostream>
#include <cmath>
#include <boost/timer/timer.hpp>
#include <iomanip>

MPISolver::MPISolver(int Nx, int Ny, int Nz, int Px, int Py, int Pz, double epsilon):
Nx(Nx),
Ny(Ny),
Nz(Nz),
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
    int nx = Nx - 2;
    int ny = Ny - 2;
    int nz = Nz - 2;
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

    int lNx = (Pi < remainder_x) ? quotient_x + 1 : quotient_x;
    int lNy = (Pj < remainder_y) ? quotient_y + 1 : quotient_y;
    int lNz = (Pk < remainder_z) ? quotient_z + 1 : quotient_z;

    lNx += 2;
    lNy += 2;
    lNz += 2;

    // int Start_i = (Pi < remainder_x) ? Pi * (quotient_x + 1) : Pi * quotient_x + remainder_x;
    // int Start_j = (Pj < remainder_y) ? Pj * (quotient_y + 1) : Pj * quotient_y + remainder_y;
    // int Start_k = (Pk < remainder_z) ? Pk * (quotient_z + 1) : Pk * quotient_z + remainder_z;

    double* lf = new double[lNx*lNy*lNz];
    for (int i = 0; i < lNx*lNy*lNz; i++){lf[i]=6;}

    localstate = std::make_unique<JacobiLocalState>(Nx, Ny, Nz, lNx, lNy, lNz, lf);

    for (int i = 0; i < 6; i++){
        if(neighbours_ranks[i] != MPI_PROC_NULL){
            int plane_size = localstate->get_planesize(static_cast<Direction>(i));
            send_buffers[i] = new double[plane_size]();
            recv_buffers[i] = new double[plane_size]();
        } 
    }
    
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
    if (gridrank == 0) std::cout << "Converged in " << iterations << " iterations, residual: " << residual << std::endl;
}

MPISolver::~MPISolver(){
    for (int i = 0; i < 6; i++){
        if(neighbours_ranks[i] != MPI_PROC_NULL){
            delete[] send_buffers[i];
            delete[] recv_buffers[i];
        }
    }
}

double MPISolver::get_residual(){
    double local_residual = localstate->get_residualsquared();
    double residual_square_sum = 0;
    MPI_Allreduce(&local_residual, &residual_square_sum, 1, MPI_DOUBLE, MPI_SUM, gridcomm);
    return sqrt(residual_square_sum);
}