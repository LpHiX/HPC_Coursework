/**
 * @file mpi_solver.cpp
 * @author Martin Leung
 */

#include "mpi_solver.h"
#include <mpi.h>
#include <iostream>
#include <algorithm>

MPISolver::MPISolver(int Nx, int Ny, int Nz, int Px, int Py, int Pz, double epsilon):
Nx(Nx),
Ny(Ny),
Nz(Nz),
Px(Px),
Py(Py),
Pz(Pz),
epsilon(epsilon)
{
    MPI_Comm gridcomm;
    const int dims = 3;
    int sizes[dims] = {Px, Py, Pz};
    int periods[dims] = {0, 0, 0};
    int reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &gridcomm);

    int rank;
    MPI_Comm_rank(gridcomm, &rank);

    MPI_Cart_shift(gridcomm, 0, 1, &neighbours_ranks[NEG_X], &neighbours_ranks[POS_X]);
    MPI_Cart_shift(gridcomm, 1, 1, &neighbours_ranks[NEG_Y], &neighbours_ranks[POS_Y]);
    MPI_Cart_shift(gridcomm, 2, 1, &neighbours_ranks[NEG_Z], &neighbours_ranks[POS_Z]);

    // Handle remainders

    int quotient_x = Nx / Px;
    int quotient_y = Ny / Py;
    int quotient_z = Nz / Pz;
    int remainder_x = Nx % Px;
    int remainder_y = Ny % Py;
    int remainder_z = Nz % Pz;

    int coords[3];
    MPI_Cart_coords(gridcomm, rank, 3, coords);
    int Pi = coords[0];
    int Pj = coords[1];
    int Pk = coords[2];

    int lNx = (Pi < remainder_x) ? quotient_x + 1 : quotient_x;
    int lNy = (Pj < remainder_y) ? quotient_y + 1 : quotient_y;
    int lNz = (Pk < remainder_z) ? quotient_z + 1 : quotient_z;

    int Start_i = (Pi < remainder_x) ? Pi * (quotient_x + 1) : Pi * quotient_x + remainder_x;
    int Start_j = (Pj < remainder_y) ? Pj * (quotient_y + 1) : Pj * quotient_y + remainder_y;
    int Start_k = (Pk < remainder_z) ? Pk * (quotient_z + 1) : Pk * quotient_z + remainder_z;

    double* lf = new double[lNx*lNy*lNz];
    for (int i = 0; i < lNx*lNy*lNz; i++){lf[i]=6;}

    localstate = std::make_unique<JacobiLocalState>(Nx, Ny, Nz, lNx, lNy, lNz, lf);

    
    // transfer_array = new double[lNy*lNz];
}

void MPISolver::solve(){
    int iterations = 0;
    double residual = localstate->get_residual();

    while (residual > epsilon){
        localstate->jacobi_step();
        
        for (int i = 0; i < 6; i++){
            
        }
    }
}